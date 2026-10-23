#!/usr/bin/env python3
"""Independently audit real OpenFOAM nonorthogonal mesh & ABC solver outputs.

Every exact velocity/pressure reference is computed at OpenFOAM's own
writeCellCentres centroids. Never use uniform cell-index coordinates for
a geometrically warped case.
"""
from pathlib import Path
import math
import re
import sys

VEC=re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
SCA=re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
VEC3=re.compile(r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)")
BUDGET=re.compile(r"CTEAS_M1B1 max_div_phi=([-+0-9.eE]+) boundary_net_flux=([-+0-9.eE]+) kinetic_energy=([-+0-9.eE]+)")
SCALARS=re.compile(r"CTEAS_M1A integral=([-+0-9.eE]+)")

def field(case,time,name,n,kind):
    path=Path(case)/time/name
    if not path.exists():
        raise SystemExit(f"FAIL: OpenFOAM real field not found: {path}")
    m=(VEC if kind=="vector" else SCA).search(path.read_text())
    if not m or int(m.group(1))!=n**3:
        raise SystemExit(f"FAIL: wrong actual field length/type: {path}")
    a=[tuple(map(float,z)) for z in VEC3.findall(m.group(2))] if kind=="vector" else [float(x) for x in m.group(2).split()]
    if len(a)!=n**3 or any(not math.isfinite(v) for z in a for v in (z if kind=="vector" else [z])):
        raise SystemExit(f"FAIL: malformed/nonfinite output field {path}")
    return a

def measure_mesh(log):
    s=Path(log).read_text()
    def value(pattern):
        m=re.search(pattern,s,re.I)
        if not m:
            raise SystemExit(f"FAIL: missing OpenFOAM checkMesh metric {pattern}")
        return float(m.group(1))
    if "Mesh OK." not in s:
        raise SystemExit("FAIL: official checkMesh did not report Mesh OK")
    nonorth=value(r"Mesh non-orthogonality Max:\s*([-+0-9.eE]+)")
    skew=value(r"Max skewness\s*=\s*([-+0-9.eE]+)")
    minvol=value(r"Min volume\s*=\s*([-+0-9.eE]+)")
    totalvol=value(r"Total volume\s*=\s*([-+0-9.eE]+)")
    if not all(map(math.isfinite,(nonorth,skew,minvol,totalvol))):
        raise SystemExit("FAIL: nonfinite checkMesh metrics")
    if not 1.0<nonorth<65 or not 0<=skew<3.5 or minvol<=0 or abs(totalvol-1)>1e-6:
        raise SystemExit("FAIL: warp mesh quality/volume does not meet predeclared gates")
    return nonorth,skew,minvol,totalvol

def check(case,meshlog,solverlog,n):
    nonorth,skew,minvol,totalvol=measure_mesh(meshlog)
    C=field(case,"0","C",n,"vector")
    initialU=field(case,"0","U",n,"vector")
    initialP=field(case,"0","p",n,"scalar")
    U=field(case,"0.05","U",n,"vector")
    p=field(case,"0.05","p",n,"scalar")
    k=2*math.pi
    decay=math.exp(-.01*k*k*.05)
    expectedU=[]
    expectedP=[]
    initialMismatch=0.
    centreDeformation=0.
    for index,(x,y,z) in enumerate(C):
        i=index%n;j=(index//n)%n;l=index//(n*n)
        centreDeformation=max(centreDeformation,abs(x-(i+.5)/n))
        u=math.sin(k*z)+math.cos(k*y)
        v=math.sin(k*x)+math.cos(k*z)
        w=math.sin(k*y)+math.cos(k*x)
        truth=(u,v,w)
        initialMismatch=max(initialMismatch,max(abs(a-b) for a,b in zip(initialU[index],truth)),
                            abs(initialP[index]+0.5*(u*u+v*v+w*w)))
        expectedU.append(tuple(decay*v for v in truth))
        expectedP.append(-0.5*decay*decay*(u*u+v*v+w*w))
    if not 0.02<centreDeformation<0.09 or initialMismatch>5e-12:
        raise SystemExit("FAIL: OpenFOAM actual centroids not used for exact initial fields")
    velocityError=max(abs(a-b) for v,w in zip(U,expectedU) for a,b in zip(v,w))
    meanp=sum(p)/len(p);meanref=sum(expectedP)/len(expectedP)
    pressureError=max(abs((a-meanp)-(b-meanref)) for a,b in zip(p,expectedP))
    logs=Path(solverlog).read_text()
    steps=[tuple(map(float,z)) for z in BUDGET.findall(logs)]
    scalars=[float(z) for z in SCALARS.findall(logs)]
    if len(steps)!=20 or len(scalars)!=20:
        raise SystemExit("FAIL: missing actual 20-step pressure/transport budgets")
    maxdiv=max(abs(x[0]) for x in steps)
    maxnet=max(abs(x[1]) for x in steps)
    scalarDrift=max(abs(z-1) for z in scalars)
    energy=steps[-1][2]
    print(f"M1B-1f nonorth_ABC n={n} max_nonorth_deg={nonorth:.12g} "
          f"max_skewness={skew:.12g} min_volume={minvol:.12g} "
          f"total_volume={totalvol:.12g} centre_max_dx={centreDeformation:.12g} "
          f"initial_abc_error={initialMismatch:.12g} max_U_error={velocityError:.12g} "
          f"gauge_free_p_error={pressureError:.12g} max_div_phi={maxdiv:.12g} "
          f"max_boundary_net_flux={maxnet:.12g} max_scalar_drift={scalarDrift:.12g} "
          f"final_energy={energy:.12g}")
    if maxdiv>1e-7 or maxnet>1e-7 or scalarDrift>1e-8:
        raise SystemExit("FAIL: distorted 3D pressure projection/scalar budget")
    if not .95<energy<1.55:
        raise SystemExit("FAIL: warped ABC kinetic energy outside physical range")
    if n==8 and (velocityError>0.32 or pressureError>0.80):
        raise SystemExit("FAIL: coarse warped 3D ABC analytic U/p")
    if n==16 and (velocityError>0.22 or pressureError>0.55):
        raise SystemExit("FAIL: fine warped 3D ABC analytic U/p")
    return velocityError,pressureError,nonorth

if len(sys.argv)!=7:
    raise SystemExit("usage: check_warped_abc.py CASE8 MESH8 SOLVER8 CASE16 MESH16 SOLVER16")
u8,p8,no8=check(sys.argv[1],sys.argv[2],sys.argv[3],8)
u16,p16,no16=check(sys.argv[4],sys.argv[5],sys.argv[6],16)
if u16>=.88*u8 or p16>=.92*p8:
    raise SystemExit("FAIL: actual warped mesh U/p error does not reduce on refinement")
if min(no8,no16)<=1.0:
    raise SystemExit("FAIL: warped mesh non-orthogonality not exercised")
print(f"M1B-1f warped_ABC U_refinement_ratio={u16/u8:.12g} "
      f"p_refinement_ratio={p16/p8:.12g}")
print("PASS: real nonorthogonal 3D mesh/centroids, ABC U/p, flux/scalar, refinement")
