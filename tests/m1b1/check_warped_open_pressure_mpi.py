#!/usr/bin/env python3
"""M1B-1g-1: compare *reconstructed* real 3-D open-pressure duct fields.

All cases use identical OpenFOAM-generated warped geometry, initial fields and
fixed pressure BCs; pressure has a physical Dirichlet gauge (compare absolute p).
"""
import math
from pathlib import Path
import re
import sys

from check_warped_open_pressure import audit, read_field, FLOW, BUDGET, SCALAR

if len(sys.argv) != 11:
    raise SystemExit(
        "usage: check_warped_open_pressure_mpi.py "
        "SERIAL_CASE SERIAL_CHECK SERIAL_LOG X_CASE X_LOG "
        "Z_CASE Z_LOG ZERO_CASE ZERO_CHECK ZERO_LOG"
    )

serial_case, serial_check, serial_log, x_case, x_log, z_case, z_log, zero_case, zero_check, zero_log = sys.argv[1:]
cases = {"serial": serial_case, "x": x_case, "z": z_case, "zero-z": zero_case}
logs = {"serial": serial_log, "x": x_log, "z": z_log, "zero-z": zero_log}

# Confirm that reconstructed cases have the same *real* centroids and initial
# conditions, rather than merely comparing three implementations of a generator.
for name in ("C", "U", "p", "c"):
    vector = name in ("C", "U")
    ref = read_field(serial_case, "0", name, vector)
    # The zero-pressure negative control intentionally has p=0, not
    # p=0.1*(1-x). All other input fields, including its geometry, match.
    subset = cases if name != "p" else {k: cases[k] for k in ("serial", "x", "z")}
    for tag, case in subset.items():
        if read_field(case, "0", name, vector) != ref:
            raise SystemExit(f"FAIL: {tag} initial {name} differs from real serial case")

# Same independent analytic, mesh, signed flux, zero-drive and convergence
# checks as the M1B-1g serial gate, applied to each reconstructed MPI result.
serial_u, serial_q = audit(serial_case, serial_check, serial_log, "warped")
xu, xq = audit(x_case, serial_check, x_log, "warped")
zu, zq = audit(z_case, serial_check, z_log, "warped")
audit(zero_case, zero_check, zero_log, "zero")

serial_p = read_field(serial_case, "0.4", "p", False)
serial_c = read_field(serial_case, "0.4", "c", False)

def history(path):
    txt = Path(path).read_text()
    flux = [tuple(map(float, row)) for row in FLOW.findall(txt)]
    cont = [tuple(map(float, row)) for row in BUDGET.findall(txt)]
    scalar = [float(x) for x in SCALAR.findall(txt)]
    if len(flux) != 20 or len(cont) != 20 or len(scalar) != 20:
        raise SystemExit(f"FAIL: truncated MPI histories: {path}")
    if not all(math.isfinite(v) for row in flux + cont for v in row):
        raise SystemExit(f"FAIL: nonfinite flux/continuity history: {path}")
    if not all(math.isfinite(v) for v in scalar):
        raise SystemExit(f"FAIL: nonfinite scalar inventory: {path}")
    return flux, cont, scalar

ref_flux, ref_cont, ref_scalar = history(serial_log)
for tag, case, vel, q in (("x", x_case, xu, xq), ("z", z_case, zu, zq)):
    pp = read_field(case, "0.4", "p", False)
    cc = read_field(case, "0.4", "c", False)
    flux, cont, scalar = history(logs[tag])
    uerr = max(abs(a-b) for ua, ub in zip(serial_u, vel) for a,b in zip(ua,ub))
    perr = max(abs(a-b) for a,b in zip(serial_p, pp))  # fixed Dirichlet gauge
    cerr = max(abs(a-b) for a,b in zip(serial_c, cc))
    qerr = max(abs(a-b) for row_a,row_b in zip(ref_flux,flux)
               for a,b in zip(row_a,row_b))
    qrel = max(abs(q[1]-serial_q[1])/max(abs(serial_q[1]),1e-12),
               abs(q[0]-serial_q[0])/max(abs(serial_q[0]),1e-12))
    kinetic = abs(cont[-1][2]-ref_cont[-1][2])
    print(f"M1B-1g-1 MPI {tag} exact_mesh=1 nCells=960 "
          f"U_Linf={uerr:.12g} absolute_p_Linf={perr:.12g} "
          f"c_Linf={cerr:.12g} pressure_patch_flux_history_Linf={qerr:.12g} "
          f"final_relative_inout_flux_change={qrel:.12g} "
          f"final_kinetic_energy_difference={kinetic:.12g}")
    # Predeclared bounds, initially adapted from historical M1B-1f3
    # (0.005 m/s U, 0.02 m2/s2 p) plus a stronger physical-gauge check.
    if uerr > 0.005 or perr > 0.01 or cerr > 1e-7:
        raise SystemExit(f"FAIL: real warped pressure-duct MPI field equivalence {tag}")
    if qerr > 0.001 or qrel > 0.05 or kinetic > 0.003:
        raise SystemExit(f"FAIL: real warped pressure-duct MPI flux/energy equivalence {tag}")

xp, zp = (read_field(cases[n], "0.4", "p", False) for n in ("x", "z"))
cross_u = max(abs(a-b) for ua,ub in zip(xu,zu) for a,b in zip(ua,ub))
cross_p = max(abs(a-b) for a,b in zip(xp,zp))
print(f"M1B-1g-1 x_vs_z_U_Linf={cross_u:.12g} x_vs_z_absolute_p_Linf={cross_p:.12g}")
if cross_u > 0.005 or cross_p > 0.01:
    raise SystemExit("FAIL: nonperiodic warped 3D pressure solver decomposition orientation")
print("PASS: warped 3-D open-pressure duct serial/x-MPI/z-MPI/zero-drive pressure and flux consistency")
