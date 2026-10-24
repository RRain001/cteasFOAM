#!/usr/bin/env python3
"""M1C-5 independent source-time integration and species BE checker."""
from pathlib import Path
import math
import re
import sys

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
INV = re.compile(r"CTEAS_M1C activity_HT_Bq=(" + NUM + r") activity_HTO_Bq=(" + NUM + r")")
EM = re.compile(r"CTEAS_M1C emitted_HT_Bq=(" + NUM + r") emitted_HTO_Bq=(" + NUM + r") source_overlap_s=(" + NUM + r")")
BUD = re.compile(r"CTEAS_M1C closed_HT_residual_Bq=(" + NUM + r") closed_HTO_residual_Bq=(" + NUM + r") closed_total_residual_Bq=(" + NUM + r")")

def fail(message):
    raise SystemExit("FAIL: M1C-5 " + message)

def native(path, size, dims):
    if not path.is_file():
        fail(f"missing native field {path}")
    data = path.read_text()
    if not re.search(r"dimensions\s+\[" + re.escape(dims) + r"\]", data):
        fail(f"incorrect native field dimensions {path}")
    u = re.search(r"internalField\s+uniform\s+(" + NUM + r")\s*;", data)
    if u:
        vals = [float(u.group(1))] * size
    else:
        n = re.search(
            r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
            data, re.S)
        if not n or int(n.group(1)) != size:
            fail(f"incomplete native {size}-cell field {path}")
        vals = [float(v) for v in n.group(2).split()]
    if len(vals) != size or not all(math.isfinite(v) for v in vals):
        fail(f"nonfinite or wrong count {path}")
    return vals

def audit(case, log_path, mode):
    if mode not in ("pulse", "ht_only", "zero_rate"):
        fail("unsupported test mode")
    case, log_path = Path(case), Path(log_path)
    if not (case / "constant/polyMesh/owner").is_file():
        fail("no actual native OpenFOAM mesh")
    quality = (case / "mesh-check.log").read_text()
    if "Mesh OK." not in quality or not re.search(r"(?m)^\s*cells:\s*24\s*$", quality):
        fail("24 real OpenFOAM cells/checkMesh missing")
    volume = native(case / "0/V", 24, "0 3 0 0 0 0 0")
    if min(volume) <= 0 or abs(sum(volume)-1) > 1e-12:
        fail("invalid actual native 3D cell volumes")
    text = log_path.read_text()
    times = [float(x) for x in re.findall(r"(?m)^Time = (" + NUM + r")\s*$", text)]
    inventories = [tuple(map(float, m)) for m in INV.findall(text)]
    emissions = [tuple(map(float, m)) for m in EM.findall(text)]
    budgets = [tuple(map(float, m)) for m in BUD.findall(text)]
    if any(len(arr) != 20 for arr in (times, inventories, emissions, budgets)):
        fail("not twenty real time/inventory/emission/budget records")
    if (len(re.findall(r"Solving for HT,", text)) != 20 or
        len(re.findall(r"Solving for HTO,", text)) != 20):
        fail("not twenty HT/HTO native matrix solves")
    if not all(math.isfinite(x) for r in inventories+emissions+budgets for x in r):
        fail("nonfinite activity history")
    props=(case / "constant/transportProperties").read_text()
    q_ht = 60000
    q_hto = 12000 if mode == "pulse" else 0
    if mode == "zero_rate":
        q_ht = 0
    for key, dims, value in (
        ("sourceHT_Bq_s", "0 0 -2 0 0 0 0", q_ht),
        ("sourceHTO_Bq_s", "0 0 -2 0 0 0 0", q_hto),
        ("sourceStart", "0 0 1 0 0 0 0", .075),
        ("sourceEnd", "0 0 1 0 0 0 0", .225),
        ("kOx", "0 0 -1 0 0 0 0", .4),
        ("lambdaT", "0 0 -1 0 0 0 0", .2),
    ):
        m=re.search(r"(?m)^" + key + r"\s+\[([^]]+)\]\s+(" + NUM + r")\s*;", props)
        if not m or m.group(1) != dims or abs(float(m.group(2))-value)>1e-12:
            fail(f"wrong physical source/decay setting {key}")
    if "auditClosedActivity true;" not in props:
        fail("closed balance audit not enabled")
    for name in ("HT","HTO"):
        if max(map(abs, native(case / "0" / name, 24, "0 -3 -1 0 0 0 0"))) > 1e-14:
            fail(f"nonzero initial {name}")
    predicted_ht=predicted_hto=0.
    previous_ht=previous_hto=0.
    max_pred=max_budget=max_total=max_overlap=0.
    after_end=[]
    for i, (t, (a,z), (e_ht,e_hto,e_seconds), (r_ht,r_hto,r_total)) in enumerate(
        zip(times,inventories,emissions,budgets),1
    ):
        if abs(t-i*.05) > 2e-12:
            fail(f"wrong actual solver clock step {i}: {t}")
        seconds=max(0.,min(t,.225)-max(t-.05,.075))
        j_ht=q_ht*seconds
        j_hto=q_hto*seconds
        max_overlap=max(max_overlap, abs(e_seconds-seconds), abs(e_ht-j_ht), abs(e_hto-j_hto))
        predicted_ht=(predicted_ht+j_ht)/(1+.05*(.4+.2))
        predicted_hto=(predicted_hto+j_hto+.05*.4*predicted_ht)/(1+.05*.2)
        max_pred=max(max_pred,abs(a-predicted_ht),abs(z-predicted_hto))
        calc_ht=a-previous_ht+.05*(.4+.2)*a-j_ht
        calc_hto=z-previous_hto+.05*.2*z-.05*.4*a-j_hto
        calc_total=a+z-previous_ht-previous_hto+.05*.2*(a+z)-j_ht-j_hto
        max_budget=max(max_budget,abs(calc_ht),abs(calc_hto),
                       abs(r_ht),abs(r_hto),abs(calc_ht-r_ht),abs(calc_hto-r_hto))
        max_total=max(max_total,abs(calc_total),abs(r_total),abs(calc_total-r_total))
        if min(a,z) < -1e-12:
            fail(f"negative activity inventory step {i}")
        if t >= .3:
            after_end.append(a+z)
        previous_ht, previous_hto=a,z
    for name, expected, logged in (
        ("HT", predicted_ht, inventories[-1][0]),
        ("HTO", predicted_hto, inventories[-1][1])
    ):
        field=native(case / "1" / name, 24, "0 -3 -1 0 0 0 0")
        max_pred=max(max_pred,max(abs(v-expected) for v in field))
        if abs(sum(v*w for v,w in zip(field,volume))-logged)>1e-6:
            fail(f"native volume-weighted field inventory mismatch {name}")
    if max_overlap > 1e-8:
        fail(f"actual time-window emitted activity mismatch {max_overlap:.12g}")
    if max_pred > 1e-6 or max_budget > 1e-6 or max_total > 2e-6:
        fail(f"analytic/closed budget max BE={max_pred:.12g} budget={max_budget:.12g} total={max_total:.12g}")
    if mode == "zero_rate":
        if max(abs(v) for record in inventories for v in record)>1e-12:
            fail("zero-rate injection made activity")
    else:
        if not all(after_end[i+1]<after_end[i] for i in range(len(after_end)-1)):
            fail("post-injection radioactive total failed to decay")
        if inventories[1][0]<=0:
            fail("fractional-window initial injection missing")
    print(f"PASS: M1C-5 mode={mode} cells=24 steps=20 max_BE_Bq={max_pred:.12g} "
          f"max_budget_Bq={max_budget:.12g} max_total_Bq={max_total:.12g} "
          f"max_source_Bq={max_overlap:.12g} finalHT_Bq={predicted_ht:.12g} "
          f"finalHTO_Bq={predicted_hto:.12g}")

if len(sys.argv)!=4:
    raise SystemExit("usage: check_timed_emission.py CASE LOG pulse|ht_only|zero_rate")
audit(*sys.argv[1:])
