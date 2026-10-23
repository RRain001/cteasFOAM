#!/usr/bin/env python3
"""Independent real three-component 3D ABC serial-vs-two-axis MPI audit."""
from pathlib import Path
import math
import re
import sys

VEC_RE=re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
SC_RE=re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
PISO_RE=re.compile(
    r"CTEAS_M1B1 max_div_phi=([-+0-9.eE]+) "
    r"boundary_net_flux=([-+0-9.eE]+) kinetic_energy=([-+0-9.eE]+)"
)
SCALAR_RE=re.compile(r"CTEAS_M1A integral=([-+0-9.eE]+)")

def field(case,name,kind):
    f=Path(case)/"0.05"/name
    if not f.is_file():
        raise SystemExit(f"FAIL: missing actual reconstructed 3D {name} field {f}")
    data=f.read_text()
    m=(VEC_RE if kind=="vector" else SC_RE).search(data)
    if not m or int(m.group(1))!=512:
        raise SystemExit(f"FAIL: expected 512 real finite-volume {name} cells in {f}")
    if kind=="vector":
        vals=[tuple(map(float,v.split())) for v in re.findall(r"\(([^()]*)\)",m.group(2))]
        if any(len(v)!=3 for v in vals):
            raise SystemExit("FAIL: three velocity components required")
        flat=[v for xyz in vals for v in xyz]
    else:
        vals=[float(v) for v in m.group(2).split()]
        flat=vals
    if len(vals)!=512 or not all(math.isfinite(z) for z in flat):
        raise SystemExit("FAIL: invalid/nonfinite field")
    return vals

def budget(log):
    s=Path(log).read_text()
    history=[tuple(map(float,z)) for z in PISO_RE.findall(s)]
    inventory=[float(v) for v in SCALAR_RE.findall(s)]
    if len(history)!=20 or len(inventory)!=20:
        raise SystemExit(f"FAIL: incomplete real 20-step PISO and scalar budgets {log}")
    if not all(math.isfinite(v) for row in history for v in row) or not all(map(math.isfinite,inventory)):
        raise SystemExit("FAIL: nonfinite pressure or scalar budget")
    return history,inventory

def pressure_error(a,b):
    ma=sum(a)/len(a)
    mb=sum(b)/len(b)
    return max(abs((u-ma)-(v-mb)) for u,v in zip(a,b))

def velocity_error(a,b):
    return max(abs(u-v) for x,y in zip(a,b) for u,v in zip(x,y))

if len(sys.argv)!=7:
    raise SystemExit("usage: check_abc_mpi.py SERIAL_CASE SERIAL_LOG X_CASE X_LOG Z_CASE Z_LOG")
serialCase,serialLog,xCase,xLog,zCase,zLog=sys.argv[1:]
serialU=field(serialCase,"U","vector")
serialP=field(serialCase,"p","scalar")
serialHistory,serialScalar=budget(serialLog)
xU=field(xCase,"U","vector")
xp=field(xCase,"p","scalar")
zU=field(zCase,"U","vector")
zp=field(zCase,"p","scalar")

for name,U,p,log in (("x-split",xU,xp,xLog),("z-split",zU,zp,zLog)):
    hist,scalar=budget(log)
    du=velocity_error(U,serialU)
    dp=pressure_error(p,serialP)
    div=max(abs(h[0]) for h in hist)
    net=max(abs(h[1]) for h in hist)
    mass=max(abs(x-1) for x in scalar)
    de=abs(hist[-1][2]-serialHistory[-1][2])
    three_comp=[max(abs(U[i][d]-serialU[i][d]) for i in range(512)) for d in range(3)]
    print(f"M1B-1d-2 3D MPI {name} max_U_difference={du:.12g} "
          f"Ux_difference={three_comp[0]:.12g} "
          f"Uy_difference={three_comp[1]:.12g} "
          f"Uz_difference={three_comp[2]:.12g} "
          f"gauge_free_p_difference={dp:.12g} "
          f"max_div_phi={div:.12g} max_net_boundary_flux={net:.12g} "
          f"max_scalar_drift={mass:.12g} final_kinetic_energy_difference={de:.12g}")
    if div>1e-7 or net>1e-7 or mass>1e-8:
        raise SystemExit(f"FAIL: 3D pressure-flux or scalar conservation {name}")
    if du>0.015 or dp>0.05 or de>0.01:
        raise SystemExit(f"FAIL: 3D MPI pressure-velocity differs from serial {name}")

crossU=velocity_error(xU,zU)
crossP=pressure_error(xp,zp)
print(f"M1B-1d-2 3D orientation_U_difference={crossU:.12g} "
      f"orientation_gauge_free_p_difference={crossP:.12g}")
if crossU>0.015 or crossP>0.05:
    raise SystemExit("FAIL: x vs z decompositions have inconsistent three-component PDE")
print("PASS: actual 3D ABC pressure-corrected two-rank x/z decomposition equivalence")
