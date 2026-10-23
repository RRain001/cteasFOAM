#!/usr/bin/env python3
"""Build a 4x40x1 driven channel, no-slip walls and periodic streamwise faces."""
import sys
from pathlib import Path
if len(sys.argv)!=2:
    raise SystemExit("usage: make_poiseuille_case.py CASE_PATH")
case=Path(sys.argv[1])
for d in ("0","constant","system"):
    (case/d).mkdir(parents=True,exist_ok=True)
def write(directory,name,cls,body):
    prefix=(f'FoamFile\n{{\n  version 2.0;\n  format ascii;\n'
            f'  class {cls};\n  location "{directory}";\n'
            f'  object {name};\n}}\n')
    (case/directory/name).write_text(prefix+body)
write("system","blockMeshDict","dictionary","""
scale 1;
vertices
(
 (0 0 0) (1 0 0) (1 1 0) (0 1 0)
 (0 0 1) (1 0 1) (1 1 1) (0 1 1)
);
blocks (hex (0 1 2 3 4 5 6 7) (4 40 1) simpleGrading (1 1 1));
edges ();
boundary
(
 xMin { type cyclic; neighbourPatch xMax; transform translational;
        separationVector (1 0 0); faces ((0 4 7 3)); }
 xMax { type cyclic; neighbourPatch xMin; transform translational;
        separationVector (-1 0 0); faces ((1 2 6 5)); }
 yMin { type wall; faces ((0 1 5 4)); }
 yMax { type wall; faces ((3 7 6 2)); }
 emptyPlanes { type empty; faces ((0 3 2 1) (4 5 6 7)); }
);
mergePatchPairs ();
""")
write("system","controlDict","dictionary","""
application cteasAtmosphereFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 1;
deltaT 0.01;
writeControl timeStep;
writeInterval 100;
writeFormat ascii;
writePrecision 15;
runTimeModifiable false;
""")
write("system","fvSchemes","dictionary","""
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss upwind; div(phi,c) Gauss upwind; }
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
nu [0 2 -1 0 0 0 0] 0.05;
momentumMode pisoLaminar;
drivingAcceleration [0 1 -2 0 0 0 0] (0.1 0 0);
""")
patchU="""boundaryField
{
 xMin { type cyclic; }
 xMax { type cyclic; }
 yMin { type fixedValue; value uniform (0 0 0); }
 yMax { type fixedValue; value uniform (0 0 0); }
 emptyPlanes { type empty; }
}
"""
patchOther="""boundaryField
{
 xMin { type cyclic; }
 xMax { type cyclic; }
 yMin { type zeroGradient; }
 yMax { type zeroGradient; }
 emptyPlanes { type empty; }
}
"""
write("0","U","volVectorField",
      "dimensions [0 1 -1 0 0 0 0];\n"
      "internalField uniform (0 0 0);\n"+patchU)
write("0","p","volScalarField",
      "dimensions [0 2 -2 0 0 0 0];\n"
      "internalField uniform 0;\n"+patchOther)
write("0","c","volScalarField",
      "dimensions [0 0 0 0 0 0 0];\n"
      "internalField uniform 1;\n"+patchOther)
print("Generated 4x40x1 no-slip Poiseuille: nu=.05, a_x=.1, t=1")
