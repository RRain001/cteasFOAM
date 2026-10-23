#!/usr/bin/env python3
"""Warp only INTERIOR OpenFOAM polyMesh points, or seed ABC at actual cell centres.

The six periodic faces are exactly preserved. Operates on real OpenFOAM ASCII
polyMesh/points emitted by blockMesh, never on an assumed point-index layout.
The seed phase reads 0/C written by OpenFOAM's writeCellCentres functionObject.
"""
from pathlib import Path
import math
import re
import sys

NUM=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VEC=re.compile(r"\(\s*("+NUM+r")\s+("+NUM+r")\s+("+NUM+r")\s*\)")
FIELD=re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S)

def vectors_from_list(raw):
    m=re.search(r"(?m)^\s*(\d+)\s*\n\s*\(\s*\n",raw)
    if not m:
        raise SystemExit("FAIL: required actual OpenFOAM ASCII counted vector list missing")
    tail=raw[m.end():]
    end=re.search(r"(?m)^\s*\)\s*$",tail)
    if not end:
        raise SystemExit("FAIL: counted vector list closing delimiter absent")
    rows=VEC.findall(tail[:end.start()])
    if len(rows)!=int(m.group(1)):
        raise SystemExit("FAIL: counted vectors mismatch")
    return m,tail,end,[tuple(map(float,row)) for row in rows]

def load_centres(case,n):
    f=case/"0"/"C"
    if not f.is_file():
        raise SystemExit("FAIL: missing actual OpenFOAM 0/C from writeCellCentres")
    m=FIELD.search(f.read_text())
    if not m or int(m.group(1))!=n**3:
        raise SystemExit("FAIL: actual OpenFOAM cell centres mismatch")
    points=[tuple(map(float,x)) for x in VEC.findall(m.group(2))]
    if len(points)!=n**3 or not all(math.isfinite(c) for xyz in points for c in xyz):
        raise SystemExit("FAIL: nonfinite or missing actual cell centre")
    return points

def warp(case,n,amp):
    path=case/"constant"/"polyMesh"/"points"
    source=path.read_text()
    m,tail,end,vertices=vectors_from_list(source)
    if len(vertices)!=(n+1)**3:
        raise SystemExit(f"FAIL: expected {(n+1)**3} mesh vertices, got {len(vertices)}")
    new=[]
    moved=0
    maxshift=0.
    for x,y,z in vertices:
        dx=amp*math.sin(math.pi*x)*math.sin(math.pi*y)*math.sin(math.pi*z)
        surface=any(abs(v)<1e-12 or abs(v-1)<1e-12 for v in (x,y,z))
        if surface and abs(dx)>1e-12:
            raise SystemExit("FAIL: attempted warp of cyclic boundary point")
        if abs(dx)>1e-12:
            moved+=1
        maxshift=max(maxshift,abs(dx))
        new.append((x+dx,y,z))
    if not moved or not maxshift>0.025 or not maxshift<0.085:
        raise SystemExit("FAIL: nontrivial bounded interior mesh warp absent")
    # Explicit boundary identity and (weak) positive x-orientation checks.
    for original,updated in zip(vertices,new):
        if any(abs(v)<1e-12 or abs(v-1)<1e-12 for v in original) and original!=updated:
            raise SystemExit("FAIL: periodic boundary vertex changed")
        if not 0<=updated[0]<=1:
            raise SystemExit("FAIL: warped vertex escaped cube")
    block="\n".join(f"({x:.17g} {y:.17g} {z:.17g})" for x,y,z in new)
    # Keep FoamFile header/count intact, and only replace the vector rows.
    path.write_text(source[:m.end()]+block+"\n"+tail[end.start():])
    print(f"M1B-1f warp n={n} points={len(vertices)} moved_interior_points={moved} max_dx={maxshift:.12g} amplitude={amp:.12g}")

def set_initial(case,name,values,n,kind):
    path=case/"0"/name
    txt=path.read_text()
    start=txt.find("internalField")
    b=txt.find("boundaryField",start)
    if start<0 or b<0 or b<=start:
        raise SystemExit(f"FAIL: original OpenFOAM {name} field layout missing")
    if kind=="vector":
        middle=f"internalField nonuniform List<vector>\n{n**3}\n(\n"+\
               "\n".join(f"({x:.17g} {y:.17g} {z:.17g})" for x,y,z in values)+"\n);\n"
    else:
        middle=f"internalField nonuniform List<scalar>\n{n**3}\n(\n"+\
               "\n".join(f"{x:.17g}" for x in values)+"\n);\n"
    path.write_text(txt[:start]+middle+txt[b:])

def seed(case,n):
    centres=load_centres(case,n)
    vel=[]
    p=[]
    k=2*math.pi
    shifts=[]
    for idx,(x,y,z) in enumerate(centres):
        if not all(-1e-12<=v<=1+1e-12 for v in (x,y,z)):
            raise SystemExit("FAIL: warped OpenFOAM centroid escaped domain")
        ix=idx%n
        iy=(idx//n)%n
        iz=idx//(n*n)
        shifts.append(abs(x-(ix+.5)/n))
        u=math.sin(k*z)+math.cos(k*y)
        v=math.sin(k*x)+math.cos(k*z)
        w=math.sin(k*y)+math.cos(k*x)
        vel.append((u,v,w))
        p.append(-0.5*(u*u+v*v+w*w))
    if not 0.02<max(shifts)<0.09:
        raise SystemExit("FAIL: centres not genuinely warped or invalid index ordering")
    set_initial(case,"U",vel,n,"vector")
    set_initial(case,"p",p,n,"scalar")
    print(f"M1B-1f actual_openfoam_centroids n={n} nCells={len(centres)} max_center_dx={max(shifts):.12g} seed=ABC_at_C")

if len(sys.argv) not in (4,5):
    raise SystemExit("usage: warp_abc_mesh.py warp|seed CASE N [amplitude]")
command,folder,n_str=sys.argv[1:4]
case=Path(folder)
n=int(n_str)
if n not in (8,16):
    raise SystemExit("FAIL: only accepted 8^3 and 16^3 ABC meshes")
if command=="warp" and len(sys.argv)==5:
    amplitude=float(sys.argv[4])
    if not 0.04<=amplitude<=0.08:
        raise SystemExit("FAIL: unexpected nonorthogonality amplitude")
    warp(case,n,amplitude)
elif command=="seed" and len(sys.argv)==4:
    seed(case,n)
else:
    raise SystemExit("FAIL: invalid warp_abc_mesh.py arguments")
