#!/usr/bin/env python3
"""Generate two-resolution Taylor-Green test for real OpenFOAM PISO mode."""
from pathlib import Path
import math
import sys

if len(sys.argv)!=3:
    raise SystemExit("usage: make_tgv_case.py OUTPUT_DIR N")
case=Path(sys.argv[1])
n=int(sys.argv[2])
if n not in (12,24):
    raise SystemExit("TGV reference requires n=12 or n=24")
for d in ("0","constant","system"):
    (case/d).mkdir(parents=True,exist_ok=True)

def write(directory,name,cls,body):
    prefix=(f'FoamFile\n{{\n    version 2.0;\n    format ascii;\n'
            f'    class {cls};\n    location "{directory}";\n'
            f'    object {name};\n}}\n')
    (case/directory/name).write_text(prefix+body)

write("system","blockMeshDict","dictionary",f"""
scale 1;
vertices
(
    (0 0 0) (1 0 0) (1 1 0) (0 1 0)
    (0 0 1) (1 0 1) (1 1 1) (0 1 1)
);
blocks (hex (0 1 2 3 4 5 6 7) ({n} {n} 1) simpleGrading (1 1 1));
edges ();
boundary
(
    xMin {{ type cyclic; neighbourPatch xMax; transform translational;
            separationVector (1 0 0); faces ((0 4 7 3)); }}
    xMax {{ type cyclic; neighbourPatch xMin; transform translational;
            separationVector (-1 0 0); faces ((1 2 6 5)); }}
    yMin {{ type cyclic; neighbourPatch yMax; transform translational;
            separationVector (0 1 0); faces ((0 1 5 4)); }}
    yMax {{ type cyclic; neighbourPatch yMin; transform translational;
            separationVector (0 -1 0); faces ((3 7 6 2)); }}
    emptyPlanes {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
mergePatchPairs ();
""")
write("system","controlDict","dictionary","""
application cteasAtmosphereFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 0.1;
deltaT 0.005;
writeControl timeStep;
writeInterval 20;
writeFormat ascii;
writePrecision 15;
runTimeModifiable false;
""")
write("system","fvSchemes","dictionary","""
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) Gauss upwind;
    div(phi,c) Gauss upwind;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
fluxRequired { default no; p; }
""")
write("system","fvSolution","dictionary","""
solvers
{
    U { solver PBiCGStab; preconditioner DILU; tolerance 1e-13; relTol 0; }
    p { solver PCG; preconditioner DIC; tolerance 1e-13; relTol 0; }
    c { solver PBiCGStab; preconditioner DILU; tolerance 1e-13; relTol 0; }
}
""")
write("constant","transportProperties","dictionary","""
D [0 2 -1 0 0 0 0] 0;
nu [0 2 -1 0 0 0 0] 0.01;
momentumMode pisoLaminar;
""")
uvals=[]
pvals=[]
for j in range(n):
    y=(j+0.5)/n
    for i in range(n):
        x=(i+0.5)/n
        u=math.sin(2*math.pi*x)*math.cos(2*math.pi*y)
        v=-math.cos(2*math.pi*x)*math.sin(2*math.pi*y)
        p=+0.25*(math.cos(4*math.pi*x)+math.cos(4*math.pi*y))
        uvals.append(f"({u:.17g} {v:.17g} 0)")
        pvals.append(f"{p:.17g}")
patches="""boundaryField
{
    xMin { type cyclic; }
    xMax { type cyclic; }
    yMin { type cyclic; }
    yMax { type cyclic; }
    emptyPlanes { type empty; }
}
"""
write("0","U","volVectorField",
      "dimensions [0 1 -1 0 0 0 0];\n"
      f"internalField nonuniform List<vector>\n{n*n}\n(\n"
      +"\n".join(uvals)+"\n);\n"+patches)
write("0","p","volScalarField",
      "dimensions [0 2 -2 0 0 0 0];\n"
      f"internalField nonuniform List<scalar>\n{n*n}\n(\n"
      +"\n".join(pvals)+"\n);\n"+patches)
write("0","c","volScalarField",
      "dimensions [0 0 0 0 0 0 0];\ninternalField uniform 1;\n"
      +patches)
print(f"Generated Taylor-Green vortex n={n}, nu=.01, t_end=.1")
