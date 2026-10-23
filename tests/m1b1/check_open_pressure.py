#!/usr/bin/env python3
"""Independent physical two-pressure-patch channel audit and zero-pressure control."""
from pathlib import Path
import re
import math
import sys

vec=re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
sca=re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
vu=re.compile(r"internalField\s+uniform\s+\(([^()]*)\)\s*;")
su=re.compile(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;")
fluxpat=re.compile(r"CTEAS_M1B1 pressure_inlet_flux=([-+0-9.eE]+) pressure_outlet_flux=([-+0-9.eE]+)")
contpat=re.compile(r"CTEAS_M1B1 max_div_phi=([-+0-9.eE]+) boundary_net_flux=([-+0-9.eE]+) kinetic_energy=([-+0-9.eE]+)")
scalarpat=re.compile(r"CTEAS_M1A integral=([-+0-9.eE]+)")

def load(case,name,vector):
    f=Path(case)/"1"/name
    if not f.is_file():
        raise SystemExit(f"FAIL: missing actual OpenFOAM pressure-channel field {f}")
    data=f.read_text()
    m=(vec if vector else sca).search(data)
    if m:
        if int(m.group(1))!=320:
            raise SystemExit("FAIL: missing 320 finite-volume values")
        out=[tuple(map(float,v.split())) for v in re.findall(r"\(([^()]*)\)",m.group(2))] if vector else [float(z) for z in m.group(2).split()]
    else:
        m=(vu if vector else su).search(data)
        if not m:
            raise SystemExit(f"FAIL: no internalField in {f}")
        base=tuple(map(float,m.group(1).split())) if vector else float(m.group(1))
        out=[base]*320
    if len(out)!=320 or (vector and any(len(x)!=3 for x in out)):
        raise SystemExit("FAIL: incomplete vector/scalar field")
    if not all(math.isfinite(v) for z in out for v in (z if vector else [z])):
        raise SystemExit("FAIL: nonfinite field")
    return out

def audit(case,logfile,driven):
    U=load(case,"U",True)
    p=load(case,"p",False)
    log=Path(logfile).read_text()
    flow=[tuple(map(float,x)) for x in fluxpat.findall(log)]
    cont=[tuple(map(float,x)) for x in contpat.findall(log)]
    scalars=[float(x) for x in scalarpat.findall(log)]
    if len(flow)!=100 or len(cont)!=100 or len(scalars)!=100:
        raise SystemExit("FAIL: missing full 100-step open-boundary, continuity and scalar records")
    if not all(math.isfinite(x) for row in flow+cont for x in row):
        raise SystemExit("FAIL: nonfinite open pressure boundary flow")
    if "mode=pisoLaminar" not in log:
        raise SystemExit("FAIL: not running pressure-corrected laminar solver")
    dp=.1 if driven else 0.
    # fixedValue pressure boundary assertions in real final field, not only the generator.
    pressureTxt=(Path(case)/"1"/"p").read_text()
    for patch,expected in (("xMin",dp),("xMax",0.0)):
        pattern=rf"\b{patch}\s*\{{\s*type\s+fixedValue;\s*value\s+uniform\s+([-+0-9.eE]+)\s*;"
        m=re.search(pattern,pressureTxt)
        if not m or abs(float(m.group(1))-expected)>1e-12:
            raise SystemExit(f"FAIL: physical Dirichlet pressure lost on {patch}")
    div=max(abs(v[0]) for v in cont)
    net=max(abs(v[1]) for v in cont)
    inlet,outlet=flow[-1]
    mass=max(abs(v-1) for v in scalars)
    transverse=max(max(abs(v[1]),abs(v[2])) for v in U)
    inout=abs(inlet+outlet)
    ux=[v[0] for v in U]
    print(
        f"M1B-1e fixedPressure kind={'driven' if driven else 'zero'} "
        f"q_inlet={inlet:.12g} q_outlet={outlet:.12g} flux_mismatch={inout:.12g} "
        f"max_div_phi={div:.12g} max_net_boundary_flux={net:.12g} "
        f"max_scalar_drift={mass:.12g} max_transverse={transverse:.12g} "
        f"Ux_min={min(ux):.12g} Ux_max={max(ux):.12g}"
    )
    if div>1e-7 or net>1e-7 or inout>1e-7 or mass>1e-8:
        raise SystemExit("FAIL: fixed-pressure inlet/outlet flux or scalar budget")
    if driven:
        if inlet>=-0.005 or outlet<=0.005:
            raise SystemExit("FAIL: imposed pressure drop did not drive real throughflow")
        nu=.05
        G=.1
        def reference(y):
            return G*y*(1-y)/(2*nu) - sum(
                4*G/(nu*math.pi**3*n**3)
                *math.sin(n*math.pi*y)
                *math.exp(-nu*n*n*math.pi**2)
                for n in range(1,101,2)
            )
        uxErr=max(abs(U[8*j+i][0]-reference((j+.5)/40)) for j in range(40) for i in range(8))
        xnon=max(abs(U[8*j+i][0]-U[8*j][0]) for j in range(40) for i in range(8))
        pErr=max(abs(p[8*j+i]-dp*(1-(i+.5)/8)) for j in range(40) for i in range(8))
        print(f"M1B-1e fixedPressure driven axial_analytic_error={uxErr:.12g} "
              f"streamwise_velocity_nonuniformity={xnon:.12g} "
              f"pressure_linear_error={pErr:.12g}")
        if uxErr>0.02 or transverse>0.005 or xnon>0.015 or pErr>0.025:
            raise SystemExit("FAIL: physical pressure boundary flow solution/analytical profile")
    else:
        if max(abs(x) for x in ux)>1e-8 or transverse>1e-8:
            raise SystemExit("FAIL: zero pressure difference still drives flow")
        if abs(inlet)>1e-8 or abs(outlet)>1e-8:
            raise SystemExit("FAIL: zero pressure control has nonzero flow")

if len(sys.argv)!=5:
    raise SystemExit("usage: check_open_pressure.py DRIVEN_CASE DRIVEN_LOG ZERO_CASE ZERO_LOG")
audit(sys.argv[1],sys.argv[2],True)
audit(sys.argv[3],sys.argv[4],False)
print("PASS: actual pressure fixedValue open-inlet/outlet, pressure-driven flow, negative control")
