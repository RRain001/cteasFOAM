#!/usr/bin/env python3
"""M1B-1g-2: true finite-volume time-step self-convergence audit.

Identical OpenFOAM warped physical mesh and initial U/p/c; only deltaT differs.
The independent transient Stokes duct expression tests absolute physical
accuracy, while successive numerical differences test time-discretization
consistency. Neither establishes spatial convergence or field validation.
"""
import math
from pathlib import Path
import re
import sys

from check_warped_open_pressure import (
    BUDGET, FLOW, SCALAR, audit_physical_points, duct_velocity, quality, read_field,
)

if len(sys.argv) != 8:
    raise SystemExit(
        "usage: check_warped_open_pressure_timestep.py "
        "BASE_CHECK CASE20 LOG20 CASE10 LOG10 CASE05 LOG05"
    )
mesh_log = sys.argv[1]
cases = [sys.argv[2], sys.argv[4], sys.argv[6]]
logs = [sys.argv[3], sys.argv[5], sys.argv[7]]
dts = [0.02, 0.01, 0.005]
steps = [20, 40, 80]

# Baseline mesh quality and all three exact physical inputs.
nonorth, skew, minvol = quality(mesh_log, warped=True)
moved = audit_physical_points(cases[0], warped=True)
if moved != 665:
    raise SystemExit("FAIL: wrong number of warped vertices")
for case in cases[1:]:
    if audit_physical_points(case, warped=True) != moved:
        raise SystemExit("FAIL: different warped mesh vertex count")
    for relative in ("constant/polyMesh/points", "constant/polyMesh/faces",
                     "constant/polyMesh/owner", "constant/polyMesh/neighbour",
                     "constant/polyMesh/boundary", "0/C", "0/U", "0/p", "0/c"):
        if (Path(case)/relative).read_bytes() != (Path(cases[0])/relative).read_bytes():
            raise SystemExit(f"FAIL: candidate timestep changed physical input {relative}")

C = read_field(cases[0], "0", "C", True)
results = []
for case, log, dt, expected_steps in zip(cases, logs, dts, steps):
    U = read_field(case, "0.4", "U", True)
    p = read_field(case, "0.4", "p", False)
    c = read_field(case, "0.4", "c", False)
    txt = Path(log).read_text()
    history = [tuple(map(float, z)) for z in BUDGET.findall(txt)]
    signed = [tuple(map(float, z)) for z in FLOW.findall(txt)]
    scalar = [float(z) for z in SCALAR.findall(txt)]
    nsolves = len(re.findall(r"Solving for p,", txt))
    ncorr = len(re.findall(r"^CTEAS_M1B1 nonorthogonal_correctors=2$", txt, re.M))
    if (len(history) != expected_steps or len(signed) != expected_steps
            or len(scalar) != expected_steps or ncorr != expected_steps
            or nsolves != 9*expected_steps):
        raise SystemExit(f"FAIL: truncated real OpenFOAM PDE solve at dt={dt}")
    if "mode=pisoLaminar" not in txt or "Completed mode=pisoLaminar" not in txt:
        raise SystemExit(f"FAIL: no completed laminar pressure solve dt={dt}")
    if not all(math.isfinite(v) for z in history + signed for v in z):
        raise SystemExit(f"FAIL: nonfinite actual pressure/flux budgets dt={dt}")
    if not all(math.isfinite(v) for v in scalar):
        raise SystemExit(f"FAIL: nonfinite scalar budget dt={dt}")

    # Physical fixedValue pressure and no-slip boundary conditions must survive.
    ptext = (Path(case)/"0.4"/"p").read_text()
    utext = (Path(case)/"0.4"/"U").read_text()
    for patch, expected in (("xMin", 0.1), ("xMax", 0.0)):
        pattern = rf"\b{patch}\s*\{{\s*type\s+fixedValue;\s*value\s+uniform\s+([-+0-9.eE]+)\s*;"
        match = re.search(pattern, ptext)
        if not match or abs(float(match.group(1))-expected) > 1e-12:
            raise SystemExit(f"FAIL: physical fixedValue pressure patch {patch} dt={dt}")
    for patch in ("yMin", "yMax", "zMin", "zMax"):
        pattern = rf"\b{patch}\s*\{{\s*type\s+fixedValue;\s*value\s+uniform\s+\(0 0 0\)"
        if not re.search(pattern, utext):
            raise SystemExit(f"FAIL: missing no-slip patch {patch} dt={dt}")

    maxdiv = max(abs(row[0]) for row in history)
    net = max(abs(row[1]) for row in history)
    mismatch = max(abs(qin+qout) for qin,qout in signed)
    drift = max(abs(v-1) for v in scalar)
    uerr = max(abs(v[0]-duct_velocity(y,z))
               for v, (x,y,z) in zip(U,C))
    perror = max(abs(pv-0.1*(1-C[i][0])) for i,pv in enumerate(p))
    transverse = max(max(abs(v[1]),abs(v[2])) for v in U)
    qin,qout = signed[-1]
    print(f"M1B-1g-2 dt={dt:.6g} steps={expected_steps} psolves={nsolves} "
          f"q_in={qin:.12g} q_out={qout:.12g} "
          f"max_div_phi={maxdiv:.12g} max_net_boundary_flux={net:.12g} "
          f"flux_pair_mismatch={mismatch:.12g} inventory_drift={drift:.12g} "
          f"analytic_U_Linf={uerr:.12g} linear_p_Linf={perror:.12g} "
          f"transverse_U_Linf={transverse:.12g} "
          f"nonorth_deg={nonorth:.12g} skew={skew:.12g} minvol={minvol:.12g}")
    if (maxdiv > 1e-7 or net > 1e-7 or mismatch > 1e-7
            or drift > 1e-8 or max(abs(z-1) for z in c) > 1e-8):
        raise SystemExit(f"FAIL: broken conservative finite-volume flux at dt={dt}")
    if qin >= -0.003 or qout <= 0.003 or uerr > 0.025 or perror > 0.03 or transverse > 0.02:
        raise SystemExit(f"FAIL: lost pressure-driven analytical duct reference dt={dt}")
    results.append((U, p, qout))

def vdiff(a,b):
    return max(abs(x-y) for av,bv in zip(a,b) for x,y in zip(av,bv))

def sdiff(a,b):
    return max(abs(x-y) for x,y in zip(a,b))

coarse_u = vdiff(results[0][0], results[1][0])
fine_u = vdiff(results[1][0], results[2][0])
coarse_p = sdiff(results[0][1], results[1][1])
fine_p = sdiff(results[1][1], results[2][1])
coarse_q = abs(results[0][2]-results[1][2])
fine_q = abs(results[1][2]-results[2][2])

def ratio(coarse, fine, kind):
    if not all(math.isfinite(v) for v in (coarse,fine)):
        raise SystemExit(f"FAIL: nonfinite temporal differences for {kind}")
    if coarse <= 1e-11:
        raise SystemExit(f"FAIL: insufficient measurable dt sensitivity in {kind}")
    return fine/coarse

u_ratio = ratio(coarse_u, fine_u, "U")
q_ratio = ratio(coarse_q, fine_q, "outlet Q")
print(f"M1B-1g-2 temporal_U_delta_20_10={coarse_u:.12g} "
      f"temporal_U_delta_10_05={fine_u:.12g} U_ratio={u_ratio:.12g} "
      f"temporal_p_delta_20_10={coarse_p:.12g} "
      f"temporal_p_delta_10_05={fine_p:.12g} "
      f"temporal_outlet_Q_delta_20_10={coarse_q:.12g} "
      f"temporal_outlet_Q_delta_10_05={fine_q:.12g} "
      f"outlet_Q_ratio={q_ratio:.12g}")
# Temporal differences should reduce substantially under dt halving.
# Ratios under 0.85 are intentionally weaker than the asymptotic
# first-order Euler expectation of ~0.5, allowing coupled PISO errors.
if u_ratio > 0.85 or q_ratio > 0.85:
    raise SystemExit("FAIL: actual warped 3D pressure duct has no temporal U/Q self-convergence")
print("PASS: same-mesh genuine 3D warped fixed-pressure duct time-step self-convergence")
