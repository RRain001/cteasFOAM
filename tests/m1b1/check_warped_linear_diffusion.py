#!/usr/bin/env python3
"""M1B-1g-8: external affine Laplace solution and signed FVM face budgets.

Exact continuum reference: u=0, p=0, c=1-x in [0,1]^3, D=0.01.
Fixed c=1/0 on x ends, insulating y/z walls. This is an invariance
test of nonorthogonal finite-volume diffusion, not external field validation.
"""
from pathlib import Path
import math
import re
import sys
from check_warped_open_pressure import (
    read_field, audit_physical_points, quality, FLOW, BUDGET
)
FLUX = re.compile(
    r"CTEAS_M1A scalar_flux_xMin=([-+0-9.eE]+) "
    r"scalar_flux_xMax=([-+0-9.eE]+) "
    r"scalar_flux_net=([-+0-9.eE]+) "
    r"scalar_budget_residual=([-+0-9.eE]+)"
)
MASS = re.compile(r"CTEAS_M1A integral=([-+0-9.eE]+)")
INITIAL = re.compile(r"CTEAS_M1A initial_integral=([-+0-9.eE]+)")

if len(sys.argv) != 4:
    raise SystemExit("usage: check_warped_linear_diffusion.py CASE CHECKMESH_LOG SOLVER_LOG")
case, mesh_log, solver_log = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
if audit_physical_points(case, True) != 665:
    raise SystemExit("FAIL: not the prescribed warped physical mesh")
nonorth, skew, minvol = quality(mesh_log, True)
properties = (case / "constant/transportProperties").read_text()
if len(re.findall(r"(?m)^D \[0 2 -1 0 0 0 0\] 0\.01;$", properties)) != 1:
    raise SystemExit("FAIL: analytic D=0.01 is not the actual transport coefficient")
if properties.count("auditPassiveScalarFlux true;") != 1:
    raise SystemExit("FAIL: native matrix-flux audit disabled")
boundary = (case / "0/c").read_text()
for patch, value in (("xMin", 1), ("xMax", 0)):
    pat = rf"\b{patch}\s*\{{\s*type\s+fixedValue;\s*value\s+uniform\s+{value};"
    if not re.search(pat, boundary):
        raise SystemExit(f"FAIL: wrong analytical Dirichlet patch {patch}")
for patch in ("yMin", "yMax", "zMin", "zMax"):
    if not re.search(rf"\b{patch}\s*\{{\s*type\s+zeroGradient;", boundary):
        raise SystemExit(f"FAIL: noninsulating analytical wall {patch}")
C = read_field(case, "0", "C", True)
seed = read_field(case, "0", "c", False)
final = read_field(case, "0.4", "c", False)
U = read_field(case, "0.4", "U", True)
p = read_field(case, "0.4", "p", False)
reference = [1-x for x,y,z in C]
initial_err = max(abs(a-b) for a,b in zip(seed, reference))
max_err = max(abs(a-b) for a,b in zip(final, reference))
rms_err = math.sqrt(sum((a-b)**2 for a,b in zip(final, reference))/len(final))
if initial_err > 1e-12:
    raise SystemExit(f"FAIL: initial native-centroid affine reference corrupted {initial_err}")
if max_err > 1e-5 or rms_err > 1e-6:
    raise SystemExit(f"FAIL: analytical c=1-x mismatch Linf={max_err} RMS={rms_err}")
if max(abs(v) for vel in U for v in vel) > 1e-8 or max(abs(v) for v in p) > 1e-8:
    raise SystemExit("FAIL: spurious momentum/pressure under zero drive")
data = solver_log.read_text()
if "mode=pisoLaminar" not in data:
    raise SystemExit("FAIL: wrong physical solver mode")
if len(re.findall(r"(?m)^CTEAS_M1B1 nonorthogonal_correctors=2$", data)) != 20:
    raise SystemExit("FAIL: missing genuine 20 PISO steps")
if len(re.findall(r"Solving for p,", data)) != 180:
    raise SystemExit("FAIL: missing genuine 180 p solves")
flux = [[float(x) for x in row] for row in FLUX.findall(data)]
mass = [float(x) for x in MASS.findall(data)]
start = [float(x) for x in INITIAL.findall(data)]
flow = [[float(x) for x in row] for row in FLOW.findall(data)]
p_budget = [[float(x) for x in row] for row in BUDGET.findall(data)]
if any(len(a) != 20 for a in (flux, mass, flow, p_budget)) or len(start) != 1:
    raise SystemExit("FAIL: incomplete actual 20-step finite-volume budgets")
if not all(math.isfinite(x) for x in [*mass,*start,*(v for rows in (flux,flow,p_budget) for row in rows for v in row)]):
    raise SystemExit("FAIL: nonfinite logged numerical result")
if abs(start[0]-0.5) > 1e-6:
    raise SystemExit(f"FAIL: wrong analytical initial inventory {start[0]}")
previous = start[0]
worst_balance = 0.0
worst_flux = 0.0
for step, (m, (left,right,net,resid), q, pb) in enumerate(zip(mass,flux,flow,p_budget),1):
    err = abs(m-previous+0.02*net)
    worst_balance = max(worst_balance,err,abs(resid))
    worst_flux = max(worst_flux,abs(left+0.01),abs(right-0.01))
    if abs(left+0.01)>2e-5 or abs(right-0.01)>2e-5:
        raise SystemExit(f"FAIL: analytical F=-D grad(c)·n flux mismatch at step {step}: {left}, {right}")
    if abs(net)>2e-5 or abs(net-left-right)>1e-6:
        raise SystemExit(f"FAIL: nonconservative signed boundary flux at step {step}")
    if err>2e-6 or abs(resid)>1e-7:
        raise SystemExit(f"FAIL: time-discrete mass budget at step {step}: {err}, {resid}")
    if abs(m-0.5)>1e-5:
        raise SystemExit(f"FAIL: affine steady inventory drift at step {step}: {m}")
    if max(abs(x) for x in q)>1e-8 or abs(pb[0])>1e-7 or abs(pb[1])>1e-7:
        raise SystemExit(f"FAIL: zero-drive volume flux/divergence at step {step}")
    previous=m
print(f"M1B-1g-8 cells=960 nonorth_deg={nonorth:.9g} skew={skew:.9g} "
      f"initial_Linf={initial_err:.12g} final_Linf={max_err:.12g} "
      f"final_RMS={rms_err:.12g} max_flux_error={worst_flux:.12g} "
      f"max_mass_balance={worst_balance:.12g} final_inventory={mass[-1]:.12g}")
print("PASS: real warped 3-D corrected FVM positive diffusion vs independent affine solution")
