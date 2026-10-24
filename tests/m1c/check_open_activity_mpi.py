#!/usr/bin/env python3
"""M1C-4 independent serial vs two-rank native open-pressure species parity."""
from pathlib import Path
import math
import re
import sys

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
FIELD = re.compile(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s+(\d+)\s*\((.*?)\)\s*;", re.S)
VEC = re.compile(r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)")
INV = re.compile(r"CTEAS_M1C activity_HT_Bq=(" + NUM + r") activity_HTO_Bq=(" + NUM + r")")
HT = re.compile(r"CTEAS_M1C flux_HT_xMin_Bq_s=(" + NUM + r") flux_HT_xMax_Bq_s=(" + NUM + r") flux_HT_net_Bq_s=(" + NUM + r")")
HTO = re.compile(r"CTEAS_M1C flux_HTO_xMin_Bq_s=(" + NUM + r") flux_HTO_xMax_Bq_s=(" + NUM + r") flux_HTO_net_Bq_s=(" + NUM + r")")
FLOW = re.compile(r"CTEAS_M1B1 pressure_inlet_flux=(" + NUM + r") pressure_outlet_flux=(" + NUM + r")")

def fail(message):
    raise SystemExit("FAIL: M1C-4 " + message)

def field(path, name):
    if not path.is_file():
        fail(f"missing genuine native field {path}")
    text = path.read_text()
    if name in ("HT", "HTO") and not re.search(
        r"dimensions\s+\[0 -3 -1 0 0 0 0\]", text
    ):
        fail(f"invalid Bq/m3 dimensions {path}")
    m = FIELD.search(text)
    if m:
        if int(m.group(2)) != 960:
            fail(f"incorrect native 960-cell field size {path}")
        if name == "U":
            if m.group(1) != "vector":
                fail(f"incorrect vector field {path}")
            values = [tuple(map(float, item)) for item in VEC.findall(m.group(3))]
        else:
            if m.group(1) != "scalar":
                fail(f"incorrect scalar field {path}")
            values = [(float(x),) for x in m.group(3).split()]
    else:
        # OpenFOAM legitimately compresses identically zero final fields
        # to native 'internalField uniform 0' instead of List<scalar>.
        # The independent M1C-3 check and source-mesh check already verify
        # all 960 actual cells; only permit the exact uniform syntax here.
        if name == "U":
            u = re.search(r"internalField\s+uniform\s+\(([^()]*)\)\s*;", text)
            if not u:
                fail(f"no valid native uniform vector field {path}")
            values = [tuple(map(float, u.group(1).split()))] * 960
        else:
            u = re.search(r"internalField\s+uniform\s+(" + NUM + r")\s*;", text)
            if not u:
                fail(f"no valid native uniform scalar field {path}")
            values = [(float(u.group(1)),)] * 960
    if len(values) != 960 or not all(math.isfinite(x) for row in values for x in row):
        fail(f"nonfinite / incomplete native {name} field")
    return values

def history(log, tag):
    if not log.is_file():
        fail(f"missing native OpenFOAM log {log}")
    text = log.read_text()
    if tag == "parallel" and not re.search(r"nProcs\s*:\s*2", text):
        fail(f"not actual two-rank MPI {log}")
    if (len(re.findall(r"(?m)^Time = ", text)) != 20
        or len(re.findall(r"Solving for p,", text)) != 180
        or len(re.findall(r"Solving for HT,", text)) != 20
        or len(re.findall(r"Solving for HTO,", text)) != 20):
        fail(f"incomplete true PISO species steps in {log}")
    out = {}
    for name, pattern in (("inventory", INV), ("htflux", HT),
                          ("htoflux", HTO), ("flow", FLOW)):
        rows = [tuple(map(float, row)) for row in pattern.findall(text)]
        if len(rows) != 20 or any(not math.isfinite(x) for r in rows for x in r):
            fail(f"missing or nonfinite 20-step {name} in {log}")
        out[name] = rows
    return out

def compare(mode, work):
    serial = work / f"serial-{mode}"
    parallel = work / f"parallel-{mode}"
    if (serial / "constant/polyMesh/points").read_bytes() != (
        parallel / "constant/polyMesh/points"
    ).read_bytes():
        fail(f"{mode} serial/MPI mesh source changed")
    for i in range(2):
        processor = parallel / f"processor{i}"
        if not (processor / "constant/polyMesh/owner").is_file():
            fail(f"{mode} processor{i} decomposed mesh missing")
        for name in ("U", "HT", "HTO"):
            if not (processor / "0.4" / name).is_file():
                fail(f"{mode} processor{i} missing final {name}")
    s_hist = history(work / f"serial-{mode}.log", "serial")
    p_hist = history(work / f"parallel-{mode}.log", "parallel")
    for key, max_delta in (("inventory", 20.0), ("htflux", 25.0),
                           ("htoflux", 25.0), ("flow", 1e-6)):
        delta = max(abs(a - b) for sr, pr in zip(s_hist[key], p_hist[key])
                    for a, b in zip(sr, pr))
        if delta > max_delta:
            fail(f"{mode} {key} serial/MPI delta={delta:.12g} limit={max_delta}")
        print(f"M1C-4 {mode} max_{key}_delta={delta:.12g}")
    for name, linf_limit, rms_limit in (("HT", 100.0, 30.0),
                                        ("HTO", 20.0, 6.0),
                                        ("U", 0.003, 0.003)):
        a = field(serial / "0.4" / name, name)
        b = field(parallel / "0.4" / name, name)
        diffs = [x - y for ra, rb in zip(a, b) for x, y in zip(ra, rb)]
        linf = max(map(abs, diffs))
        rms = math.sqrt(sum(x*x for x in diffs) / len(diffs))
        print(f"M1C-4 {mode} {name} max_Linf={linf:.12g} RMS={rms:.12g}")
        if linf > linf_limit or rms > rms_limit:
            fail(f"{mode} {name} native field Linf={linf:.12g} RMS={rms:.12g}")
    if mode == "zero":
        for name in ("HT", "HTO", "U"):
            for component in field(parallel / "0.4" / name, name):
                if max(map(abs, component)) > 1e-8:
                    fail(f"zero-pressure MPI {name} generated nonphysical field")

if len(sys.argv) != 2:
    raise SystemExit("usage: check_open_activity_mpi.py WORK_DIRECTORY")
work = Path(sys.argv[1])
for mode in ("warped", "reverse", "zero"):
    compare(mode, work)
print("PASS: M1C-4 x-decomposed official v2606 pressure-open HT/HTO MPI parity")
