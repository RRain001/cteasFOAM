#!/usr/bin/env python3
"""Independent serial/MPI reconstructed field and inventory comparison."""
from pathlib import Path
import math
import re
import sys

FIELD_RE = re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", re.S)
INTEGRAL_RE = re.compile(r"CTEAS_M1A (?:initial_integral|integral)=([-+0-9.eE]+)")

def field(case):
    path=Path(case)/"2"/"c"
    if not path.exists():
        raise SystemExit(f"FAIL: missing solver final field {path}")
    m=FIELD_RE.search(path.read_text())
    if not m:
        raise SystemExit(f"FAIL: invalid solver field {path}")
    values=[float(s) for s in m.group(2).split()]
    if int(m.group(1))!=40 or len(values)!=40 or not all(map(math.isfinite,values)):
        raise SystemExit("FAIL: invalid reconstructed 40-cell result")
    return values

def integrals(filename):
    values=[float(s) for s in INTEGRAL_RE.findall(Path(filename).read_text())]
    if len(values)<201 or not all(map(math.isfinite, values)):
        raise SystemExit(f"FAIL: missing 200-step solver budget from {filename}")
    return values

if len(sys.argv)!=5:
    raise SystemExit("usage: check_parallel_diffusion.py SERIAL_CASE PARALLEL_CASE SERIAL_LOG PARALLEL_LOG")
s,p=field(sys.argv[1]),field(sys.argv[2])
si,pi=integrals(sys.argv[3]),integrals(sys.argv[4])
max_field=max(abs(a-b) for a,b in zip(s,p))
max_inventory=max(max(abs(x-si[0]) for x in si),max(abs(x-pi[0]) for x in pi))
final_budget=abs(si[-1]-pi[-1])
analytic=[
    1+0.2*math.exp(-0.05*math.pi**2*2)*math.cos(math.pi*(i+0.5)/40)
    for i in range(40)
]
max_analytic=max(abs(a-b) for a,b in zip(p,analytic))
print(
    f"M1A MPI max_cell_difference={max_field:.12g} "
    f"max_inventory_drift={max_inventory:.12g} "
    f"final_budget_difference={final_budget:.12g} "
    f"parallel_analytic_error={max_analytic:.12g}"
)
if max_field>1e-9 or final_budget>1e-9:
    raise SystemExit("FAIL: MPI and serial fields or budgets differ")
if max_inventory>2e-8 or max_analytic>1.5e-3:
    raise SystemExit("FAIL: MPI inventory or analytic benchmark")
print("PASS: actual two-rank MPI PDE and serial finite-volume equivalence")
