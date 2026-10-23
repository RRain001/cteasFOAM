#!/usr/bin/env python3
"""M1B-1g-4: official OpenFOAM icoFoam cross-implementation reference.

No shared finite-volume C++ solver implementation. Both genuinely
integrate the same physical fixed-pressure 3D duct. This compares
reconstructed native fields; it does not assert external scientific
validation or formal independence of research reviewers.
"""
import math
from pathlib import Path
import re
import sys

from check_warped_open_pressure import (
    audit, duct_velocity, read_field, audit_physical_points, quality,
)

if len(sys.argv) != 10:
    raise SystemExit("usage: check_icofoam_cross_solver.py CT_CASE IC_CASE MESH_LOG CT_LOG IC_LOG CT_ZERO IC_ZERO CT_ZERO_LOG IC_ZERO_LOG")
ct, official, meshlog, ctlog, icolog, zero_ct, zero_ic, zctlog, zicolog = sys.argv[1:]
quality(meshlog, True)
assert audit_physical_points(ct, True) == 665

for a, b in ((ct, official), (zero_ct, zero_ic)):
    for relative in ("constant/polyMesh/points",
                     "constant/polyMesh/faces",
                     "constant/polyMesh/owner",
                     "constant/polyMesh/neighbour",
                     "constant/polyMesh/boundary",
                     "0/C", "0/U", "0/p", "0/c",
                     "constant/transportProperties", "system/fvSchemes",
                     "system/fvSolution"):
        pa,pb = Path(a)/relative, Path(b)/relative
        if pa.read_bytes() != pb.read_bytes():
            raise SystemExit(f"FAIL: cross-solver physical input mismatch {relative}")
    control_a = (Path(a)/"system/controlDict").read_text()
    control_b = (Path(b)/"system/controlDict").read_text()
    expected = control_a.replace("application cteasAtmosphereFoam;", "application icoFoam;")
    if expected == control_a or control_b != expected:
        raise SystemExit("FAIL: cross-solver controlDict differs beyond application")
    if "PISO" not in (Path(a)/"system/fvSolution").read_text():
        raise SystemExit("FAIL: no explicit reference PISO configuration")

for label, ref_log in (("driven", icolog), ("zero", zicolog)):
    text = Path(ref_log).read_text()
    steps = re.findall(r"(?m)^Time = ([0-9.]+)\s*$", text)
    pressure_solves = len(re.findall(r"Solving for p,", text))
    if len(steps) != 20 or abs(float(steps[-1])-0.4)>1e-12:
        raise SystemExit(f"FAIL: incomplete upstream icoFoam {label} history")
    if pressure_solves != 180:
        raise SystemExit(f"FAIL: upstream icoFoam not running 3x3 PISO pressure solves ({label}, count={pressure_solves})")
    if "End" not in text:
        raise SystemExit(f"FAIL: upstream icoFoam not completed {label}")
    print(f"M1B-1g-4 upstream {label} steps={len(steps)} p_solves={pressure_solves}")

# The existing independent pressure/flux/analytic checker audits the
# CTEAS solver, including true signed inlet/outlet and zero-drive.
audit(ct, meshlog, ctlog, "warped")
audit(zero_ct, meshlog, zctlog, "zero")

C = read_field(ct,"0","C",True)
U = read_field(ct,"0.4","U",True)
Ur = read_field(official,"0.4","U",True)
p = read_field(ct,"0.4","p",False)
pr = read_field(official,"0.4","p",False)
u_err = max(max(abs(a-b) for a,b in zip(v,w)) for v,w in zip(U,Ur))
p_err = max(abs(a-b) for a,b in zip(p,pr))
u_rms = math.sqrt(math.fsum(sum((a-b)**2 for a,b in zip(v,w)) for v,w in zip(U,Ur))/len(U))
max_u = max(max(abs(x) for x in u) for u in U)
analytic = max(abs(u[0]-duct_velocity(c[1],c[2])) for u,c in zip(Ur,C))
ref_transverse = max(max(abs(v[1]),abs(v[2])) for v in Ur)
print(f"M1B-1g-4 driven cross_solver_U_Linf={u_err:.12g} "
      f"cross_solver_U_RMS={u_rms:.12g} "
      f"cross_solver_absolute_p_Linf={p_err:.12g} "
      f"CTEAS_U_max={max_u:.12g} icoFoam_stokes_U_Linf={analytic:.12g} "
      f"icoFoam_transverse_U_Linf={ref_transverse:.12g}")
# A priori broad parity limits for separate legitimate PISO formulations;
# tighter-than-model analytic comparison is not assumed or retrofitted.
if (u_err > 0.007 or p_err > 0.007
        or u_rms > 0.003 or analytic > 0.025 or ref_transverse > 0.02):
    raise SystemExit("FAIL: official icoFoam and experimental solver diverge beyond predeclared bounds")

for name, candidate, reference in (("U",zero_ct,zero_ic),("p",zero_ct,zero_ic)):
    vector = name=="U"
    a = read_field(candidate,"0.4",name,vector)
    b = read_field(reference,"0.4",name,vector)
    max_a = max(abs(x) for v in a for x in (v if vector else (v,)))
    max_b = max(abs(x) for v in b for x in (v if vector else (v,)))
    if max_a > 1e-8 or max_b > 1e-8:
        raise SystemExit(f"FAIL: nonzero spontaneous {name} in independent upstream zero drive")
    print(f"M1B-1g-4 zero_drive_{name}_CTEAS_Linf={max_a:.12g} "
          f"zero_drive_{name}_icoFoam_Linf={max_b:.12g}")

print("PASS: physical-input-identical native v2606 upstream icoFoam pressure-duct cross-solver comparison")
