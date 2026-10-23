#!/usr/bin/env python3
"""Independent no-slip unsteady Poiseuille Fourier reference and budgets."""
from pathlib import Path
import math
import re
import sys
if len(sys.argv)!=3:
    raise SystemExit("usage: check_poiseuille.py CASE_PATH SOLVER_LOG")
case=Path(sys.argv[1])
log=Path(sys.argv[2]).read_text()
path=case/"1"/"U"
if not path.is_file():
    raise SystemExit(f"FAIL: missing real OpenFOAM output {path}")
m=re.search(
    r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",
    path.read_text(),re.S
)
if not m or int(m.group(1))!=160:
    raise SystemExit("FAIL: wrong OpenFOAM cell list for 4x40 channel")
actual=[tuple(map(float,s.split())) for s in re.findall(r"\(([^()]*)\)",m.group(2))]
if len(actual)!=160 or any(len(v)!=3 for v in actual):
    raise SystemExit("FAIL: invalid vector field")
if not all(math.isfinite(x) for v in actual for x in v):
    raise SystemExit("FAIL: nonfinite velocity field")
nu=.05
G=.1
def expected(y,t=1):
    steady=G/(2*nu)*y*(1-y)
    transient=sum(
        4*G/(nu*math.pi**3*n**3)
        *math.sin(n*math.pi*y)
        *math.exp(-nu*n*n*math.pi**2*t)
        for n in range(1,101,2)
    )
    return steady-transient
# blockMesh indexing: i advances before j in this structured 4x40x1 case.
max_error=max(
    abs(actual[4*j+i][0]-expected((j+.5)/40))
    for j in range(40) for i in range(4)
)
max_transverse=max(max(abs(v[1]),abs(v[2])) for v in actual)
homogeneity=max(
    abs(actual[4*j+i][0]-actual[4*j][0])
    for j in range(40) for i in range(4)
)
ux=[v[0] for v in actual]
continuity=re.findall(
    r"CTEAS_M1B1 max_div_phi=([-+0-9.eE]+) "
    r"boundary_net_flux=([-+0-9.eE]+) "
    r"kinetic_energy=([-+0-9.eE]+)",log
)
inventory=[float(x) for x in re.findall(r"CTEAS_M1A integral=([-+0-9.eE]+)",log)]
if len(continuity)!=100 or len(inventory)!=100:
    raise SystemExit("FAIL: incomplete 100-step physical and scalar budgets")
if not all(math.isfinite(x) for row in continuity for x in map(float,row)):
    raise SystemExit("FAIL: nonfinite pressure corrected continuity audit")
max_div=max(abs(float(row[0])) for row in continuity)
max_net=max(abs(float(row[1])) for row in continuity)
max_scalar=max(abs(z-1.0) for z in inventory)
print(f"M1B-1c Poiseuille max_U_error={max_error:.12g} "
      f"max_transverse={max_transverse:.12g} "
      f"max_streamwise_nonuniformity={homogeneity:.12g} "
      f"max_div_phi={max_div:.12g} "
      f"max_boundary_net_flux={max_net:.12g} "
      f"max_passive_scalar_drift={max_scalar:.12g} "
      f"min_U={min(ux):.12g} max_U={max(ux):.12g}")
if max_error>0.003 or max_transverse>1e-8 or homogeneity>1e-7:
    raise SystemExit("FAIL: transient no-slip Poiseuille analytic profile")
if max_div>1e-7 or max_net>1e-7 or max_scalar>1e-8:
    raise SystemExit("FAIL: pressure corrected continuity or scalar conservation")
if not 0 < min(ux) < max(ux) < 0.25:
    raise SystemExit("FAIL: wrong driven flow direction or bounds")
print("PASS: actual pressure-driven-equivalent no-slip transient Poiseuille")
