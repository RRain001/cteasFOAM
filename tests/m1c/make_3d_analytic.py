#!/usr/bin/env python3
"""M1C-11 native 3-D periodic source-free HT/HTO analytic transport case."""
from pathlib import Path
import math
import re
import sys

VELOCITY={"advect":(.2,.1,.07),"diffuse":(0.,0.,0.)}
Q=2*math.pi
ACTIVITY_DIMS="0 -3 -1 0 0 0 0"
NATIVE_VECTOR=re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
NUM=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
XYZ=re.compile(r"\(\s*("+NUM+r")\s+("+NUM+r")\s+("+NUM+r")\s*\)")

def dictionary_header(where,name):
    return (f'FoamFile\n{{\n version 2.0;\n format ascii;\n class dictionary;\n '
            f'location "{where}";\n object {name};\n}}\n')

def change(path,old,new):
    t=path.read_text()
    if t.count(old)!=1:raise SystemExit(f"FAIL: M1C-11 missing or repeated input {old}: {path}")
    path.write_text(t.replace(old,new))

def writefield(case,name,dimension,initial,cyclic):
    patch="\n".join(f" {p} {{ type cyclic; }}" for p in ("xMin","xMax","yMin","yMax","zMin","zMax"))
    kind="volVectorField" if name=="U" else "volScalarField"
    (case/"0"/name).write_text(
        f'FoamFile\n{{\n version 2.0;\n format ascii;\n class {kind};\n location "0";\n object {name};\n}}\n'
        f"dimensions [{dimension}];\ninternalField uniform {initial};\nboundaryField\n{{\n{patch}\n}}\n"
    )

def prepare(case,mode):
    if mode not in VELOCITY:raise SystemExit("FAIL: invalid M1C-11 flow mode")
    bm=case/"system/blockMeshDict"
    for axis,vec,faces in [
       ("x","1 0 0",("0 4 7 3","1 2 6 5")),
       ("y","0 1 0",("0 1 5 4","3 7 6 2")),
       ("z","0 0 1",("0 3 2 1","4 5 6 7")),
    ]:
        lo,hi=f"{axis}Min",f"{axis}Max"
        defaulttype="patch" if axis=="x" else "wall"
        change(bm,f"{lo} {{ type {defaulttype}; faces (({faces[0]})); }}",
               f"{lo} {{ type cyclic; neighbourPatch {hi}; transform translational; separationVector ({vec}); faces (({faces[0]})); }}")
        change(bm,f"{hi} {{ type {defaulttype}; faces (({faces[1]})); }}",
               f"{hi} {{ type cyclic; neighbourPatch {lo}; transform translational; separationVector ({' '.join(str(-int(v)) for v in vec.split())}); faces (({faces[1]})); }}")
    (case/"constant/transportProperties").write_text(
        dictionary_header("constant","transportProperties")+"D [0 2 -1 0 0 0 0] 0;\n"
        "momentumMode prescribed;\nactivityMode htHto;\n"
        "kOx [0 0 -1 0 0 0 0] 0.4;\n"
        "lambdaT [0 0 -1 0 0 0 0] 0.2;\n"
        "DHT [0 2 -1 0 0 0 0] 0.01;\n"
        "DHTO [0 2 -1 0 0 0 0] 0.02;\n"
        "auditActivityFlux true;\n"
    )
    ctl=case/"system/controlDict"
    change(ctl,"endTime 0.4;","endTime 0.2;")
    change(ctl,"deltaT 0.02;","deltaT 0.005;")
    change(ctl,"writeInterval 20;","writeInterval 40;")
    (case/"system/fvSchemes").write_text(
        dictionary_header("system","fvSchemes")+"ddtSchemes { default Euler; }\n"
        "gradSchemes { default leastSquares; }\n"
        "divSchemes { default none; div(phi,c) Gauss upwind;"
        " div(phi,HT) Gauss upwind; div(phi,HTO) Gauss upwind; }\n"
        "laplacianSchemes { default Gauss linear corrected; }\n"
        "interpolationSchemes { default linear; }\n"
        "snGradSchemes { default corrected; }\n"
        "fluxRequired { default no; c; HT; HTO; }\n"
    )
    (case/"system/fvSolution").write_text(
        dictionary_header("system","fvSolution")+"solvers\n{\n"
        " c { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }\n"
        " HT { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }\n"
        " HTO { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }\n"
        "}\n"
    )
    u=VELOCITY[mode]
    writefield(case,"U","0 1 -1 0 0 0 0",f"({u[0]} {u[1]} {u[2]})",True)
    writefield(case,"c","0 0 0 0 0 0 0","1",True)
    writefield(case,"HT",ACTIVITY_DIMS,"1000000",True)
    writefield(case,"HTO",ACTIVITY_DIMS,"200000",True)
    # Momentum is prescribed; no physical pressure field used or loaded.
    (case/"0/p").unlink()
    print(f"M1C-11 prepared {mode}, all six cyclic patches, no Bq source, U={u}")

def seed(case,n):
    m=NATIVE_VECTOR.search((case/"0/C").read_text())
    if not m or int(m.group(1))!=n:raise SystemExit("FAIL: missing genuine native cell centres")
    cent=[tuple(map(float,t)) for t in XYZ.findall(m.group(2))]
    if len(cent)!=n or not all(0<x<1 for row in cent for x in row):
        raise SystemExit("FAIL: native centres not physical")
    f=[math.prod(math.cos(Q*v) for v in row) for row in cent]
    for name,mean,amp in [("HT",1e6,2e5),("HTO",2e5,5e4)]:
        path=case/"0"/name
        txt=path.read_text()
        i,j=txt.index("internalField"),txt.index("boundaryField")
        vals="\n".join(f"{mean+amp*v:.17g}" for v in f)
        path.write_text(txt[:i]+f"internalField nonuniform List<scalar>\n{n}\n(\n{vals}\n);\n"+txt[j:])
    print(f"M1C-11 seeded genuine native {n} cells for 3-D cos-product modal HT/HTO")

if __name__=="__main__":
    if len(sys.argv)<3:raise SystemExit("usage: make_3d_analytic.py prepare CASE advect|diffuse | seed CASE N")
    act,case=sys.argv[1],Path(sys.argv[2])
    if act=="prepare" and len(sys.argv)==4:prepare(case,sys.argv[3])
    elif act=="seed" and len(sys.argv)==4:seed(case,int(sys.argv[3]))
    else:raise SystemExit("FAIL: M1C-11 invalid invocation")
