#!/usr/bin/env python3
"""Independent same-input original icoFoam parity and scalar-budget audit."""
from pathlib import Path
import math
import re
import sys

NUM=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VEC=re.compile(r"\(\s*("+NUM+r")\s+("+NUM+r")\s+("+NUM+r")\s*\)")
FIELDS={
 "U":re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S),
 "p":re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",re.S),
 "c":re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",re.S),
 "C":re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S),
}
FLUX=re.compile(r"CTEAS_M1A scalar_flux_xMin=("+NUM+r") scalar_flux_xMax=("+NUM+r") scalar_flux_net=("+NUM+r") scalar_budget_residual=("+NUM+r")")
PISO=re.compile(r"CTEAS_M1B1 max_div_phi=("+NUM+r") boundary_net_flux=("+NUM+r") kinetic_energy=("+NUM+r")")

def field(case,t,name,n):
    path=case/t/name
    s=path.read_text()
    m=FIELDS[name].search(s)
    vector=name in ("U","C")
    if m:
        if int(m.group(1))!=n: raise SystemExit(f"FAIL: wrong field count {path}")
        vals=[tuple(map(float,v)) for v in VEC.findall(m.group(2))] if vector else [float(v) for v in m.group(2).split()]
    else:
        pattern=(r"internalField\s+uniform\s+\(([^()]*)\)\s*;" if vector else
                 r"internalField\s+uniform\s+("+NUM+r")\s*;")
        m=re.search(pattern,s)
        if not m: raise SystemExit(f"FAIL: missing field {path}")
        v=tuple(map(float,m.group(1).split())) if vector else float(m.group(1))
        vals=[v]*n
    if len(vals)!=n or not all(math.isfinite(v) for row in vals for v in (row if vector else (row,))):
        raise SystemExit(f"FAIL: invalid/nonfinite field {path}")
    return vals

def audit(ct,ref,meshlog,ctlog,reflog,n,lims):
    files=("constant/polyMesh/points","constant/polyMesh/faces","constant/polyMesh/owner",
           "constant/polyMesh/neighbour","constant/polyMesh/boundary","0/U","0/p","0/c",
           "constant/transportProperties","system/fvSchemes","system/fvSolution")
    for f in files:
        if (ct/f).read_bytes()!=(ref/f).read_bytes():
            raise SystemExit(f"FAIL: original input mismatch {f}")
    cctl=(ct/"system/controlDict").read_text()
    rctl=(ref/"system/controlDict").read_text()
    if rctl!=cctl.replace("application cteasAtmosphereFoam;","application icoFoam;") or rctl==cctl:
        raise SystemExit("FAIL: reference controlDict differs beyond executable")
    meshtext=meshlog.read_text()
    if "Mesh OK." not in meshtext: raise SystemExit("FAIL: official checkMesh")
    match=re.search(r"Mesh non-orthogonality Max:\s*("+NUM+r")",meshtext)
    if not match or not (1<float(match.group(1))<65):
        raise SystemExit("FAIL: physical warped mesh non-orthogonality")
    props=(ct/"constant/transportProperties").read_text()
    for term in ("nu [0 2 -1 0 0 0 0] 0.05;","D [0 2 -1 0 0 0 0] 0.01;",
                 "drivingAcceleration [0 1 -2 0 0 0 0] (0 0 0);",
                 "momentumMode pisoLaminar;","nPressureCorrectors 3;",
                 "nNonOrthogonalCorrectors 2;","auditPassiveScalarFlux true;"):
        if term not in props: raise SystemExit(f"FAIL: coefficient mismatch {term}")
    for log,tag in ((ctlog,"native"),(reflog,"original")):
        s=log.read_text()
        if len(re.findall(r"(?m)^Time = ",s))!=40 or len(re.findall(r"Solving for p,",s))!=360:
            raise SystemExit(f"FAIL: incomplete {tag} 40-step 3x3 PISO")
        if tag=="native" and (len(re.findall(r"Solving for c,",s))!=40 or
                              len(re.findall(r"Solving for Ux,",s))!=40):
            raise SystemExit("FAIL: incomplete native momentum/scalar")
        if tag=="original" and "End" not in s:
            raise SystemExit("FAIL: upstream icoFoam incomplete")
    uc=field(ct,"0.4","U",n); ur=field(ref,"0.4","U",n)
    pc=field(ct,"0.4","p",n); pr=field(ref,"0.4","p",n)
    for name,arr in (("native",uc),("original",ur)):
        err=max(abs(v-t) for row in arr for v,t in zip(row,(0.25,0.,0.)))
        if err>2e-5: raise SystemExit(f"FAIL: {name} U analytical uniformity {err}")
    for name,p in (("native",pc),("original",pr)):
        if max(p)-min(p)>2e-5: raise SystemExit(f"FAIL: {name} p not spatially uniform")
    du=max(abs(v-w) for a,b in zip(uc,ur) for v,w in zip(a,b))
    dp=max(abs((v-pc[0])-(w-pr[0])) for v,w in zip(pc,pr))
    if du>2e-5 or dp>2e-5:
        raise SystemExit(f"FAIL: original icoFoam periodic PISO U/p parity U={du:.12g} p={dp:.12g}")
    piso=[tuple(map(float,x)) for x in PISO.findall(ctlog.read_text())]
    if len(piso)!=40 or max(abs(x) for row in piso for x in row[:2])>1e-8:
        raise SystemExit("FAIL: native PISO divergence/boundary budget")
    centres=field(ct,"0","C",n)
    c0=field(ct,"0","c",n); cf=field(ct,"0.4","c",n)
    init=[0.5+0.25*math.cos(2*math.pi*x) for x,y,z in centres]
    exact=[0.5+0.25*math.exp(-0.01*(2*math.pi)**2*0.4)*
           math.cos(2*math.pi*(x-0.25*0.4)) for x,y,z in centres]
    if max(abs(a-b) for a,b in zip(c0,init))>1e-12:
        raise SystemExit("FAIL: original parity scalar initial centre field")
    err=[abs(a-b) for a,b in zip(cf,exact)]
    linf=max(err); rms=math.sqrt(sum(e*e for e in err)/n)
    if linf>lims[0] or rms>lims[1] or not all(0.23<=v<=0.77 for v in cf):
        raise SystemExit(f"FAIL: periodic PISO scalar analytical error Linf={linf:.12g} RMS={rms:.12g}")
    s=ctlog.read_text()
    flux=[tuple(map(float,v)) for v in FLUX.findall(s)]
    mass=[float(v) for v in re.findall(r"CTEAS_M1A integral=("+NUM+r")",s)]
    initial=[float(v) for v in re.findall(r"CTEAS_M1A initial_integral=("+NUM+r")",s)]
    if len(flux)!=40 or len(mass)!=40 or len(initial)!=1:
        raise SystemExit("FAIL: incomplete scalar flux budget")
    old=initial[0]; worst=0.
    for step,(new,(left,right,net,resid)) in enumerate(zip(mass,flux),1):
        e=abs(new-old+0.01*net); worst=max(worst,e,abs(resid))
        if (e>2e-6 or abs(resid)>1e-7 or abs(net)>1e-8 or
            abs(left+right)>1e-8 or min(abs(left),abs(right))<=1e-4):
            raise SystemExit(f"FAIL: periodic PISO scalar mass budget step={step} e={e}")
        old=new
    print(f"M1B-1g-12 cells={n} original_native_U_Linf={du:.12g} "
          f"p_gauge_Linf={dp:.12g} scalar_Linf={linf:.12g} "
          f"scalar_RMS={rms:.12g} max_scalar_budget={worst:.12g}")
    return rms

if len(sys.argv)!=11:
    raise SystemExit("usage: check_warped_periodic_upstream_piso.py CT_COARSE REF_COARSE MESH_COARSE CTLOG_COARSE REFLOG_COARSE CT_FINE REF_FINE MESH_FINE CTLOG_FINE REFLOG_FINE")
c=audit(*[Path(x) for x in sys.argv[1:6]],960,(0.055,0.040))
f=audit(*[Path(x) for x in sys.argv[6:11]],3240,(0.040,0.030))
ratio=f/c if c>0 else float("inf")
print(f"M1B-1g-12 fine_coarse_scalar_RMS_ratio={ratio:.12g}")
if not (0<f and 0<c and ratio<=0.85):
    raise SystemExit("FAIL: periodic upstream PISO scalar spatial refinement")
print("PASS: identical-input official v2606 unforced periodic icoFoam U/p parity + native scalar")
