#!/usr/bin/env python3
"""Independent parser of M1B-1g-6 real OpenFOAM scalar finite-volume budgets."""
import math
from pathlib import Path
import re
import sys
from check_warped_open_pressure import (
    audit_physical_points, quality, read_field, FLOW, BUDGET
)

FLUX = re.compile(
    r"CTEAS_M1A scalar_flux_xMin=([-+0-9.eE]+) "
    r"scalar_flux_xMax=([-+0-9.eE]+) "
    r"scalar_flux_net=([-+0-9.eE]+) "
    r"scalar_budget_residual=([-+0-9.eE]+)"
)
MASS = re.compile(r"CTEAS_M1A integral=([-+0-9.eE]+)")
INITIAL = re.compile(r"CTEAS_M1A initial_integral=([-+0-9.eE]+)")
DT = 0.02


def finite(*values):
    return all(math.isfinite(x) for x in values)


def audit(case, log, mode):
    case, log = Path(case), Path(log)
    data = log.read_text()
    if "mode=pisoLaminar" not in data:
        raise SystemExit(f"FAIL: no genuine laminar PISO {mode}")
    if len(re.findall(r"(?m)^CTEAS_M1B1 nonorthogonal_correctors=2$", data)) != 20:
        raise SystemExit(f"FAIL: missing native PISO steps {mode}")
    if len(re.findall(r"Solving for p,", data)) != 180:
        raise SystemExit(f"FAIL: missing real p matrix solves {mode}")
    rows = [[float(x) for x in row] for row in FLUX.findall(data)]
    masses = [float(x) for x in MASS.findall(data)]
    pressure = [[float(x) for x in row] for row in FLOW.findall(data)]
    p_budgets = [[float(x) for x in row] for row in BUDGET.findall(data)]
    starts = [float(x) for x in INITIAL.findall(data)]
    if any(len(x) != 20 for x in (rows, masses, pressure, p_budgets)) or len(starts) != 1:
        raise SystemExit(f"FAIL: missing 20 signed-flux/mass/pressure records {mode}")
    if not finite(starts[0], *masses, *(v for row in rows+pressure+p_budgets for v in row)):
        raise SystemExit(f"FAIL: nonfinite budget {mode}")
    if abs(starts[0]) > 1e-12:
        raise SystemExit(f"FAIL: nonzero tracer at start {mode}")
    if audit_physical_points(case, True) != 665:
        raise SystemExit(f"FAIL: changed native warped polyMesh {mode}")

    C = read_field(case, "0", "C", True)
    seed = read_field(case, "0", "c", False)
    output = read_field(case, "0.4", "c", False)
    U = read_field(case, "0.4", "U", True)
    p = read_field(case, "0.4", "p", False)
    if len(C) != 960 or max(abs(v) for v in seed) > 1e-12:
        raise SystemExit(f"FAIL: native c initial zero field absent {mode}")
    minc, maxc = min(output), max(output)
    if minc < -1e-9 or maxc > 1 + 1e-9:
        raise SystemExit(f"FAIL: negative/overshooting scalar {mode}")
    patch = (case/"0"/"c").read_text()
    for name, expected in (("xMin", 0 if mode == "reverse" else 1),
                           ("xMax", 1 if mode == "reverse" else 0)):
        m = re.search(rf"\b{name}\s*\{{\s*type\s+inletOutlet;\s*"
                      r"inletValue\s+uniform\s+([-+0-9.eE]+)", patch)
        if not m or abs(float(m.group(1))-expected)>1e-12:
            raise SystemExit(f"FAIL: wrong physical inletOutlet tracer patch {mode}/{name}")
    for wall in ("yMin", "yMax", "zMin", "zMax"):
        if not re.search(rf"\b{wall}\s*\{{\s*type\s+zeroGradient", patch):
            raise SystemExit(f"FAIL: scalar wall condition lost {mode}/{wall}")
    dp = 0 if mode == "zero" else (-.1 if mode == "reverse" else .1)
    pressure_text = (case/"0.4"/"p").read_text()
    for end, value in (("xMin",dp),("xMax",0)):
        m = re.search(rf"\b{end}\s*\{{\s*type\s+fixedValue;\s*"
                      r"value\s+uniform\s+([-+0-9.eE]+)", pressure_text)
        if not m or abs(float(m.group(1))-value)>1e-12:
            raise SystemExit(f"FAIL: missing real fixed pressure patch {mode}/{end}")
    if max(abs(c[0]) for c in C) > 1:
        raise SystemExit("FAIL: native centroid outside domain")

    prev = starts[0]
    max_balance, max_recomputed, max_pair_flux, max_div = 0, 0, 0, 0
    for i, (m, f, pr, pb) in enumerate(zip(masses, rows, pressure, p_budgets)):
        x_min, x_max, net, residual = f
        max_balance = max(max_balance, abs(residual))
        max_recomputed = max(max_recomputed, abs(m-prev+DT*net))
        max_pair_flux = max(max_pair_flux, abs(pr[0]+pr[1]))
        max_div = max(max_div, abs(pb[0]))
        if abs(net-(x_min+x_max)) > 1e-6:
            raise SystemExit(f"FAIL: unexpected non-wall scalar face flux {mode}/{i}")
        if abs(residual) > 1e-7:
            raise SystemExit(f"FAIL: scalar budget residual {mode}/{i}: {residual}")
        if abs(m-prev+DT*net)>2e-6:
            raise SystemExit(f"FAIL: independent time-discrete scalar balance {mode}/{i}")
        if m < prev - 1e-8 or m < -1e-8 or m > 1+1e-8:
            raise SystemExit(f"FAIL: nonphysical/inconsistent tracer inventory {mode}/{i}")
        prev = m
    if max_pair_flux > 1e-7 or max_div > 1e-7:
        raise SystemExit(f"FAIL: real 3D pressure-corrected volume imbalance {mode}")
    x_min, x_max, net, _ = rows[-1]
    q_min,q_max = pressure[-1]
    if mode == "zero":
        if max(abs(v) for row in rows for v in row)>1e-8 or (
                max(abs(v) for v in masses)>1e-8 or
                max(abs(v) for v in p)>1e-8 or
                max(abs(v) for vel in U for v in vel)>1e-8 or
                max(abs(v) for row in pressure for v in row)>1e-8):
            raise SystemExit("FAIL: tracer spontaneously generated with zero flow")
    else:
        if not masses[-1] > .001 or not net < -.001:
            raise SystemExit(f"FAIL: genuine inward tracer flux absent {mode}")
        if mode == "warped" and not (x_min < -.001 and q_min < -.003 and q_max > .003):
            raise SystemExit("FAIL: forward tracer/volume flow sign")
        if mode == "reverse" and not (x_max < -.001 and q_min > .003 and q_max < -.003):
            raise SystemExit("FAIL: reversed tracer/volume flow sign")
    print(f"M1B-1g-6 mode={mode} native_cells=960 inventory_final={masses[-1]:.12g} "
          f"tracer_min={minc:.12g} tracer_max={maxc:.12g} "
          f"scalar_flux_xMin={x_min:.12g} scalar_flux_xMax={x_max:.12g} "
          f"scalar_flux_net={net:.12g} q_in={q_min:.12g} q_out={q_max:.12g} "
          f"budget_max={max_balance:.12g} recomputed_max={max_recomputed:.12g} "
          f"max_div_phi={max_div:.12g}")
    return masses[-1]


if len(sys.argv)!=8:
    raise SystemExit("usage: check_open_pressure_tracer.py CHECK_LOG FWD_CASE FWD_LOG REV_CASE REV_LOG ZERO_CASE ZERO_LOG")
qlog, fwd, flog, reverse, rlog, zero, zlog = sys.argv[1:]
quality(qlog, True)
for a,b in ((fwd,reverse),(fwd,zero)):
    for rel in ("constant/polyMesh/points","constant/polyMesh/faces",
                "constant/polyMesh/owner","constant/polyMesh/neighbour",
                "constant/polyMesh/boundary","0/C","0/U"):
        if (Path(a)/rel).read_bytes() != (Path(b)/rel).read_bytes():
            raise SystemExit(f"FAIL: unequal real geometry/seed U {rel}")
forward_mass = audit(fwd,flog,"warped")
reverse_mass = audit(reverse,rlog,"reverse")
zero_mass = audit(zero,zlog,"zero")
if abs(zero_mass) > 1e-8:
    raise SystemExit("FAIL: no-flow negative control tracer inventory")
print("PASS: native 3D nonuniform scalar inlet advection, signed mass budget, reverse/zero controls")
