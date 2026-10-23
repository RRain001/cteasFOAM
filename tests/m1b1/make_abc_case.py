#!/usr/bin/env python3
"""OpenFOAM genuine 3D ABC/Beltrami periodic Navier-Stokes exact-solution case.

ABC A=B=C=1 and k=2*pi on unit cube:
  U = (sin(k*z)+cos(k*y), sin(k*x)+cos(k*z), sin(k*y)+cos(k*x))
  curl(U)=k*U; p=-0.5*|U|^2 (kinematic pressure), up to additive constant.
  U(t)=U(0)*exp(-nu*k*k*t) and p(t)=p(0)*exp(-2*nu*k*k*t).
"""
from pathlib import Path
import math
import sys

if len(sys.argv)!=3:
    raise SystemExit("usage: make_abc_case.py OUTPUT_CASE N")
case=Path(sys.argv[1])
n=int(sys.argv[2])
if n not in (8,16):
    raise SystemExit("ABC accepted verification grids are 8^3 and 16^3")
for d in ("0","constant","system"):
    (case/d).mkdir(parents=True,exist_ok=True)

def write(d,name,cls,body):
    head=(f'FoamFile\n{{\n  version 2.0;\n  format ascii;\n'
          f'  class {cls};\n  location "{d}";\n  object {name};\n}}\n')
    (case/d/name).write_text(head+body)

write("system","blockMeshDict","dictionary",f"""
scale 1;
vertices
(
 (0 0 0) (1 0 0) (1 1 0) (0 1 0)
 (0 0 1) (1 0 1) (1 1 1) (0 1 1)
);
blocks (hex (0 1 2 3 4 5 6 7) ({n} {n} {n}) simpleGrading (1 1 1));
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
 zMin {{ type cyclic; neighbourPatch zMax; transform translational;
         separationVector (0 0 1); faces ((0 3 2 1)); }}
 zMax {{ type cyclic; neighbourPatch zMin; transform translational;
         separationVector (0 0 -1); faces ((4 5 6 7)); }}
);
mergePatchPairs ();
""")
write("system","controlDict","dictionary","""
application cteasAtmosphereFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 0.05;
deltaT 0.0025;
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
 U { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }
 p { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }
 c { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }
}
""")
write("constant","transportProperties","dictionary","""
D [0 2 -1 0 0 0 0] 0;
nu [0 2 -1 0 0 0 0] 0.01;
momentumMode pisoLaminar;
nPressureCorrectors 3;
""")
vel=[]
pressure=[]
k=2*math.pi
# OpenFOAM blockMesh cell ordering: i fast, j then, k slow.
for iz in range(n):
    z=(iz+.5)/n
    for iy in range(n):
        y=(iy+.5)/n
        for ix in range(n):
            x=(ix+.5)/n
            u=math.sin(k*z)+math.cos(k*y)
            v=math.sin(k*x)+math.cos(k*z)
            w=math.sin(k*y)+math.cos(k*x)
            vel.append(f"({u:.17g} {v:.17g} {w:.17g})")
            pressure.append(f"{(-.5*(u*u+v*v+w*w)):.17g}")

patches="""boundaryField
{
 xMin { type cyclic; }
 xMax { type cyclic; }
 yMin { type cyclic; }
 yMax { type cyclic; }
 zMin { type cyclic; }
 zMax { type cyclic; }
}
"""
write("0","U","volVectorField",
      "dimensions [0 1 -1 0 0 0 0];\n"
      f"internalField nonuniform List<vector>\n{n**3}\n(\n"
      +"\n".join(vel)+"\n);\n"+patches)
write("0","p","volScalarField",
      "dimensions [0 2 -2 0 0 0 0];\n"
      f"internalField nonuniform List<scalar>\n{n**3}\n(\n"
      +"\n".join(pressure)+"\n);\n"+patches)
write("0","c","volScalarField",
      "dimensions [0 0 0 0 0 0 0];\n"
      "internalField uniform 1;\n"+patches)
print(f"Generated true 3D ABC A=B=C=1, n={n}, nCells={n**3}, dt=.0025, t=.05")
