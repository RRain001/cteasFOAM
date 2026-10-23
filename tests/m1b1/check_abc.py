#!/usr/bin/env python3
"""Independent 3D ABC/Navier-Stokes exact flow, pressure, continuity checker."""
from pathlib import Path
import math
import re
import sys

VEC_RE=re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
SC_RE=re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
LOG_RE=re.compile(
    r"CTEAS_M1B1 max_div_phi=([-+0-9.eE]+) "
    r"boundary_net_flux=([-+0-9.eE]+) kinetic_energy=([-+0-9.eE]+)"
)
SCALAR_BUDGET_RE=re.compile(r"CTEAS_M1A integral=([-+0-9.eE]+)")

def field(case,n,name,kind):
    path=Path(case)/"0.05"/name
    if not path.is_file():
        raise SystemExit(f"FAIL: actual OpenFOAM output missing {path}")
    m=(VEC_RE if kind=="vector" else SC_RE).search(path.read_text())
    if not m or int(m.group(1))!=n**3:
        raise SystemExit(f"FAIL: malformed full {name} OpenFOAM field for n={n}")
    if kind=="vector":
        data=[tuple(float(y) for y in z.split())
              for z in re.findall(r"\(([^()]*)\)",m.group(2))]
        if any(len(v)!=3 for v in data):
            raise SystemExit("FAIL: vector field not 3-dimensional")
        flat=[v for xyz in data for v in xyz]
    else:
        data=[float(x) for x in m.group(2).split()]
        flat=data
    if len(data)!=n**3 or not all(map(math.isfinite,flat)):
        raise SystemExit(f"FAIL: invalid/nonfinite {name} field n={n}")
    return data

def check(case,logfile,n):
    U=field(case,n,"U","vector")
    p=field(case,n,"p","scalar")
    log=Path(logfile).read_text()
    history=[tuple(map(float,values)) for values in LOG_RE.findall(log)]
    scalar=[float(x) for x in SCALAR_BUDGET_RE.findall(log)]
    if len(history)!=20 or len(scalar)!=20 or "mode=pisoLaminar" not in log:
        raise SystemExit(f"FAIL: incomplete actual 20-step 3D PISO/scalar logs n={n}")
    if not all(math.isfinite(v) for row in history for v in row):
        raise SystemExit("FAIL: invalid 3D continuum budget")
    maxdiv=max(abs(row[0]) for row in history)
    maxnet=max(abs(row[1]) for row in history)
    maxmass=max(abs(v-1) for v in scalar)
    if not all(math.isfinite(v) for v in scalar):
        raise SystemExit("FAIL: nonfinite scalar activity-independent inventory")

    time=.05
    viscosity=.01
    k=2*math.pi
    amplitude=math.exp(-viscosity*k*k*time)
    pressure_amplitude=amplitude*amplitude
    uref=[]
    pref=[]
    for iz in range(n):
        z=(iz+.5)/n
        for iy in range(n):
            y=(iy+.5)/n
            for ix in range(n):
                x=(ix+.5)/n
                # Independent closed-form Beltrami field, all 3 components.
                base_u=math.sin(k*z)+math.cos(k*y)
                base_v=math.sin(k*x)+math.cos(k*z)
                base_w=math.sin(k*y)+math.cos(k*x)
                uref.append((amplitude*base_u,
                             amplitude*base_v,
                             amplitude*base_w))
                pref.append(-.5*pressure_amplitude*
                            (base_u**2+base_v**2+base_w**2))
    max_component_errors=[
        max(abs(row[d]-truth[d]) for row,truth in zip(U,uref))
        for d in range(3)
    ]
    max_u=max(max_component_errors)
    mean_p=sum(p)/len(p)
    mean_ref=sum(pref)/len(pref)
    max_p=max(abs((x-mean_p)-(y-mean_ref)) for x,y in zip(p,pref))
    kinetic_energy=sum(.5*sum(v*v for v in row) for row in U)/(n**3)
    ideal_energy=1.5*pressure_amplitude
    # Ensure all three fields are nontrivial, not a disguised 2D benchmark.
    component_ranges=[max(row[d] for row in U)-min(row[d] for row in U)
                      for d in range(3)]
    print(
        f"M1B-1d ABC n={n} max_U_error={max_u:.12g} "
        f"Ux_error={max_component_errors[0]:.12g} "
        f"Uy_error={max_component_errors[1]:.12g} "
        f"Uz_error={max_component_errors[2]:.12g} "
        f"gauge_free_p_error={max_p:.12g} "
        f"max_div_phi={maxdiv:.12g} "
        f"max_boundary_net_flux={maxnet:.12g} "
        f"max_scalar_drift={maxmass:.12g} "
        f"kinetic_energy={kinetic_energy:.12g} "
        f"exact_kinetic_energy={ideal_energy:.12g}"
    )
    if min(component_ranges)<2.0:
        raise SystemExit("FAIL: ABC field must have 3 dynamic velocity components")
    if maxdiv>1e-7 or maxnet>1e-7 or maxmass>1e-8:
        raise SystemExit("FAIL: 3D face continuity or conserved scalar")
    if not 1.0<kinetic_energy<1.55 or not kinetic_energy<1.5:
        raise SystemExit("FAIL: inconsistent 3D kinetic energy")
    if n==8 and (max_u>0.28 or max_p>0.60):
        raise SystemExit("FAIL: coarse 3D ABC exact U/p error")
    if n==16 and (max_u>0.18 or max_p>0.45):
        raise SystemExit("FAIL: refined 3D ABC exact U/p error")
    return max_u,max_p

if len(sys.argv)!=5:
    raise SystemExit("usage: check_abc.py COARSE_CASE COARSE_LOG FINE_CASE FINE_LOG")
e8,p8=check(sys.argv[1],sys.argv[2],8)
e16,p16=check(sys.argv[3],sys.argv[4],16)
print(f"M1B-1d ABC velocity_refinement_ratio={e16/e8:.12g} "
      f"pressure_refinement_ratio={p16/p8:.12g}")
if e16>=.85*e8 or p16>=.9*p8:
    raise SystemExit("FAIL: true 3D ABC exact-solution mesh refinement absent")
print("PASS: all 3 momentum components, 3D pressure, continuity and mesh refinement")
