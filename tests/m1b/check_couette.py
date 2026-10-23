#!/usr/bin/env python3
"""Check true OpenFOAM transient planar Couette against Fourier-series solution."""
import math
import re
from pathlib import Path
import sys

if len(sys.argv)!=3:
    raise SystemExit("usage: check_couette.py CASE SOLVER_LOG")
case=Path(sys.argv[1])
log=Path(sys.argv[2]).read_text()
field=case/"1"/"U"
if not field.is_file():
    raise SystemExit(f"FAIL: solver did not write U at time=1: {field}")
data=field.read_text()
m=re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",data,re.S)
if not m:
    raise SystemExit("FAIL: no nonuniform OpenFOAM U vector field")
raw_vectors=re.findall(r"\(([^()]*)\)",m.group(2))
vectors=[tuple(float(w) for w in v.split()) for v in raw_vectors]
if int(m.group(1))!=40 or len(vectors)!=40 or any(len(v)!=3 for v in vectors):
    raise SystemExit("FAIL: wrong dimensionality or cell count")
if not all(math.isfinite(w) for vec in vectors for w in vec):
    raise SystemExit("FAIL: nonfinite velocity")
def analytic(y,t=1.,nu=.05):
    steady=y
    transient=sum(
        2*((-1)**n)/(n*math.pi)*math.sin(n*math.pi*y)
        *math.exp(-nu*(n*math.pi)**2*t)
        for n in range(1,101)
    )
    return steady+transient

ref=[analytic((i+.5)/40) for i in range(40)]
max_error=max(abs(v[0]-ref[i]) for i,v in enumerate(vectors))
max_transverse=max(max(abs(v[1]),abs(v[2])) for v in vectors)
all_ux=[v[0] for v in vectors]
print(
    f"M1B Couette max_velocity_error={max_error:.12g} "
    f"max_transverse_velocity={max_transverse:.12g} "
    f"u_min={min(all_ux):.12g} u_max={max(all_ux):.12g}"
)
if "mode=laminarShear" not in log or log.count("CTEAS_M1B momentum_integral=")<200:
    raise SystemExit("FAIL: laminar momentum mode or time-step log absent")
if max_error>0.003 or max_transverse>1e-10:
    raise SystemExit("FAIL: transient Couette analytical profile or transverse field")
if min(all_ux)<-1e-8 or max(all_ux)>1+1e-8:
    raise SystemExit("FAIL: shear velocity outside boundary-constrained range")
if any(a>b+1e-8 for a,b in zip(all_ux,all_ux[1:])):
    raise SystemExit("FAIL: wrong Couette shear profile ordering")
print("PASS: solved viscous transient momentum and independent Couette eigenfunction")
