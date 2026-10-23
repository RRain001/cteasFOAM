#!/usr/bin/env python3
"""Independent D=0 versus D=0.01 comparison on real v2606 warped-duct fields.

The separate check_open_pressure_tracer.py validates each complete 3-case
set's physical signed scalar/pressure budgets before this differential test.
"""
from pathlib import Path
import math
import re
import sys
from check_warped_open_pressure import read_field, FLOW

FLUX = re.compile(
    r"CTEAS_M1A scalar_flux_xMin=([-+0-9.eE]+) "
    r"scalar_flux_xMax=([-+0-9.eE]+) "
    r"scalar_flux_net=([-+0-9.eE]+) "
    r"scalar_budget_residual=([-+0-9.eE]+)"
)
MASS = re.compile(r"CTEAS_M1A integral=([-+0-9.eE]+)")
PHYSICAL_FILES = (
    "constant/polyMesh/points", "constant/polyMesh/faces",
    "constant/polyMesh/owner", "constant/polyMesh/neighbour",
    "constant/polyMesh/boundary", "0/C", "0/U", "0/p", "0/c",
)


def check_pair(mode, zero_case, zero_log, diff_case, diff_log):
    z, d = Path(zero_case), Path(diff_case)
    for rel in PHYSICAL_FILES:
        if (z / rel).read_bytes() != (d / rel).read_bytes():
            raise SystemExit(f"FAIL: changed physical mesh or initial field {mode}/{rel}")
    for case, expected in ((z, "0"), (d, "0.01")):
        text = (case / "constant/transportProperties").read_text()
        found = re.findall(r"(?m)^D \[0 2 -1 0 0 0 0\] ([-+0-9.eE]+);$", text)
        if len(found) != 1 or float(found[0]) != float(expected):
            raise SystemExit(f"FAIL: wrong declared diffusivity {mode}/{case.name}")
    zero_data, diff_data = Path(zero_log).read_text(), Path(diff_log).read_text()
    z_flux, d_flux = FLUX.findall(zero_data), FLUX.findall(diff_data)
    z_mass, d_mass = MASS.findall(zero_data), MASS.findall(diff_data)
    z_flow, d_flow = FLOW.findall(zero_data), FLOW.findall(diff_data)
    if not all(len(a) == 20 for a in (z_flux, d_flux, z_mass, d_mass, z_flow, d_flow)):
        raise SystemExit(f"FAIL: incomplete 20-step paired flux history {mode}")
    z_m, d_m = float(z_mass[-1]), float(d_mass[-1])
    z_f, d_f = [float(x) for x in z_flux[-1]], [float(x) for x in d_flux[-1]]
    z_q, d_q = [float(x) for x in z_flow[-1]], [float(x) for x in d_flow[-1]]
    z_c, d_c = read_field(z, "0.4", "c", False), read_field(d, "0.4", "c", False)
    z_u, d_u = read_field(z, "0.4", "U", True), read_field(d, "0.4", "U", True)
    z_p, d_p = read_field(z, "0.4", "p", False), read_field(d, "0.4", "p", False)
    if any(len(a) != 960 for a in (z_c, d_c, z_u, d_u, z_p, d_p)):
        raise SystemExit(f"FAIL: incomplete 960-cell paired field {mode}")
    values = [z_m, d_m, *z_f, *d_f, *z_q, *d_q, *z_c, *d_c, *z_p, *d_p]
    values.extend(x for u in z_u + d_u for x in u)
    if not all(math.isfinite(v) for v in values):
        raise SystemExit(f"FAIL: nonfinite paired field {mode}")
    c_delta = max(abs(a-b) for a,b in zip(z_c,d_c))
    u_delta = max(abs(a-b) for u,v in zip(z_u,d_u) for a,b in zip(u,v))
    p_delta = max(abs(a-b) for a,b in zip(z_p,d_p))
    q_delta = max(abs(a-b) for a,b in zip(z_q,d_q))
    if u_delta > 1e-10 or p_delta > 1e-10 or q_delta > 1e-10:
        raise SystemExit(f"FAIL: passive tracer changed solved flow {mode}")
    if mode == "zero":
        if max(abs(x) for x in d_c) > 1e-8 or abs(d_m)>1e-8 or max(abs(x) for x in d_f)>1e-8:
            raise SystemExit("FAIL: diffusion generated tracer with zero normal flux")
    else:
        active = 1 if mode == "reverse" else 0
        gain = d_m - z_m
        influx_gain = abs(d_f[active]) - abs(z_f[active])
        if gain < 0.005 or c_delta < 0.02 or d_f[active] >= -0.001 or influx_gain < 0.005:
            raise SystemExit(
                f"FAIL: positive D sensitivity {mode} mass_gain={gain} "
                f"field_delta={c_delta} influx_gain={influx_gain}"
            )
    print(f"M1B-1g-7 mode={mode} M_D0={z_m:.12g} M_D001={d_m:.12g} "
          f"mass_gain={d_m-z_m:.12g} c_Linf_diff={c_delta:.12g} "
          f"inlet_flux_D0={z_f[1 if mode=='reverse' else 0]:.12g} "
          f"inlet_flux_D001={d_f[1 if mode=='reverse' else 0]:.12g} "
          f"U_Linf_diff={u_delta:.12g} p_Linf_diff={p_delta:.12g} "
          f"q_Linf_diff={q_delta:.12g}")


if len(sys.argv) != 13:
    raise SystemExit("usage: checker Z_FWD Z_FWD_LOG D_FWD D_FWD_LOG Z_REV Z_REV_LOG D_REV D_REV_LOG Z_ZERO Z_ZERO_LOG D_ZERO D_ZERO_LOG")
for i, mode in enumerate(("warped", "reverse", "zero")):
    check_pair(mode, *sys.argv[1+4*i:5+4*i])
print("PASS: positive diffusivity adds conservative scalar transport without flow feedback")
