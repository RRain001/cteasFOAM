#!/usr/bin/env python3
"""M1C-2 independent native serial / genuine two-rank HT/HTO parity checker."""
import math
import re
import sys
from pathlib import Path

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
SCALAR = re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", re.S)
INVENTORY = re.compile(r"CTEAS_M1C activity_HT_Bq=(" + NUM + r") activity_HTO_Bq=(" + NUM + r")")
HT_FLUX = re.compile(r"CTEAS_M1C flux_HT_xMin_Bq_s=(" + NUM + r") flux_HT_xMax_Bq_s=(" + NUM + r") flux_HT_net_Bq_s=(" + NUM + r") flux_HT_budget_residual_Bq=(" + NUM + r")")
HTO_FLUX = re.compile(r"CTEAS_M1C flux_HTO_xMin_Bq_s=(" + NUM + r") flux_HTO_xMax_Bq_s=(" + NUM + r") flux_HTO_net_Bq_s=(" + NUM + r") flux_HTO_budget_residual_Bq=(" + NUM + r") flux_total_budget_residual_Bq=(" + NUM + r")")

def fail(msg):
    raise SystemExit("FAIL: M1C-2 " + msg)

def native_field(path,n):
    if not path.is_file():
        fail(f"missing native field {path}")
    s=path.read_text()
    if not re.search(r"dimensions\s+\[0 -3 -1 0 0 0 0\]",s):
        fail(f"incorrect Bq/m3 field dimensions {path}")
    m=SCALAR.search(s)
    if m:
        if int(m.group(1))!=n:
            fail(f"native field count {path}")
        v=[float(x) for x in m.group(2).split()]
    else:
        m=re.search(r"internalField\s+uniform\s+("+NUM+r")\s*;",s)
        if not m: fail(f"unreadable native field {path}")
        v=[float(m.group(1))]*n
    if len(v)!=n or not all(math.isfinite(x) and x>=0 for x in v):
        fail(f"nonfinite/negative native field {path}")
    return v

def logs(path):
    s=path.read_text()
    if len(re.findall(r"(?m)^Time = ",s))!=40 or len(re.findall(r"Solving for HT,",s))!=40 or len(re.findall(r"Solving for HTO,",s))!=40:
        fail(f"not 40 actual species steps {path}")
    inv=[tuple(map(float,v)) for v in INVENTORY.findall(s)]
    h=[tuple(map(float,v)) for v in HT_FLUX.findall(s)]
    z=[tuple(map(float,v)) for v in HTO_FLUX.findall(s)]
    if len(inv)!=40 or len(h)!=40 or len(z)!=40:
        fail(f"incomplete native activity history {path}")
    if not all(math.isfinite(x) for row in [*inv,*h,*z] for x in row):
        fail(f"nonfinite activity log {path}")
    maxbudget=0.
    for hi,zi in zip(h,z):
        maxbudget=max(maxbudget,abs(hi[3]),abs(zi[3]),abs(zi[4]))
    if maxbudget>1e-5:
        fail(f"species/total activity budget {maxbudget:.12g} Bq")
    return inv,h,z,maxbudget

def audit(serial,parallel,slog,plog,n):
    if (serial/"constant/polyMesh/points").read_bytes() != (parallel/"constant/polyMesh/points").read_bytes():
        fail("serial/parallel source mesh not identical")
    for rank in (0,1):
        processor=parallel/f"processor{rank}"
        if not (processor/"constant/polyMesh/owner").is_file():
            fail(f"processor{rank} actual decomposed mesh missing")
        for name in ("HT","HTO"):
            if not (processor/"0.4"/name).is_file():
                fail(f"processor{rank} final {name} field missing")
    plogtext=plog.read_text()
    if "nProcs : 2" not in plogtext and "nProcs: 2" not in plogtext:
        fail("MPI solver did not report two processes")
    s_inv,s_h,s_z,s_budget=logs(slog)
    p_inv,p_h,p_z,p_budget=logs(plog)
    max_field=0.; max_rms=0.
    for name in ("HT","HTO"):
        sv=native_field(serial/"0.4"/name,n)
        pv=native_field(parallel/"0.4"/name,n)
        diffs=[abs(a-b) for a,b in zip(sv,pv)]
        linf=max(diffs)
        rms=math.sqrt(sum(v*v for v in diffs)/n)
        if linf>1e-3 or rms>5e-4:
            fail(f"{name} serial/MPI native field difference Linf={linf:.12g} RMS={rms:.12g}")
        max_field=max(max_field,linf)
        max_rms=max(max_rms,rms)
    inv_diff=max(abs(a-b) for si,pi in zip(s_inv,p_inv) for a,b in zip(si,pi))
    flux_diff=max(abs(a-b) for srows,prows in ((s_h,p_h),(s_z,p_z))
                  for sr,pr in zip(srows,prows) for a,b in zip(sr[:3],pr[:3]))
    if inv_diff>1e-4 or flux_diff>1e-3:
        fail(f"serial/MPI activity history inventory={inv_diff:.12g} Bq flux={flux_diff:.12g} Bq/s")
    print(f"M1C-2 cells={n} max_HT_HTO_Linf={max_field:.12g} Bq/m3 "
          f"max_RMS={max_rms:.12g} Bq/m3 max_inventory_delta={inv_diff:.12g} Bq "
          f"max_flux_delta={flux_diff:.12g} Bq/s "
          f"max_serial_budget={s_budget:.12g} Bq max_MPI_budget={p_budget:.12g} Bq")

if len(sys.argv)!=9:
    raise SystemExit("usage: check_mpi_activity.py COARSE_SERIAL COARSE_PARALLEL COARSE_SERIAL_LOG COARSE_MPI_LOG FINE_SERIAL FINE_PARALLEL FINE_SERIAL_LOG FINE_MPI_LOG")
audit(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),Path(sys.argv[4]),960)
audit(Path(sys.argv[5]),Path(sys.argv[6]),Path(sys.argv[7]),Path(sys.argv[8]),3240)
print("PASS: M1C-2 official v2606 two-rank native HT/HTO serial-MPI parity")
