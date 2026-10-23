#!/usr/bin/env python3
"""Generate planar Couette transient shear benchmark in OpenFOAM v2606."""
from pathlib import Path
import sys

if len(sys.argv) != 2:
    raise SystemExit("usage: make_couette_case.py OUTPUT_CASE")
case=Path(sys.argv[1])
for directory in ("0","constant","system"):
    (case/directory).mkdir(parents=True,exist_ok=True)

def write(directory,name,cls,body):
    header=(
        f'FoamFile\n{{\n    version 2.0;\n    format ascii;\n'
        f'    class {cls};\n    location "{directory}";\n'
        f'    object {name};\n}}\n'
    )
    (case/directory/name).write_text(header+body)

write("system","blockMeshDict","dictionary","""
scale 1;
vertices
(
    (0 0 0) (1 0 0) (1 1 0) (0 1 0)
    (0 0 1) (1 0 1) (1 1 1) (0 1 1)
);
blocks (hex (0 1 2 3 4 5 6 7) (1 40 1) simpleGrading (1 1 1));
edges ();
boundary
(
    yMin { type wall; faces ((0 1 5 4)); }
    yMax { type wall; faces ((3 7 6 2)); }
    xMin { type cyclic; neighbourPatch xMax; transform translational;
           separationVector (1 0 0); faces ((0 4 7 3)); }
    xMax { type cyclic; neighbourPatch xMin; transform translational;
           separationVector (-1 0 0); faces ((1 2 6 5)); }
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
deltaT 0.005;
writeControl timeStep;
writeInterval 200;
writeFormat ascii;
writePrecision 15;
runTimeModifiable false;
""")
write("system","fvSchemes","dictionary","""
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,c) Gauss upwind; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
write("system","fvSolution","dictionary","""
solvers
{
    U { solver PCG; preconditioner DIC; tolerance 1e-13; relTol 0; }
    c { solver PBiCGStab; preconditioner DILU; tolerance 1e-13; relTol 0; }
}
""")
write("constant","transportProperties","dictionary","""
D [0 2 -1 0 0 0 0] 0;
nu [0 2 -1 0 0 0 0] 0.05;
momentumMode laminarShear;
""")
write("0","U","volVectorField","""
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
    yMin { type fixedValue; value uniform (0 0 0); }
    yMax { type fixedValue; value uniform (1 0 0); }
    xMin { type cyclic; }
    xMax { type cyclic; }
    emptyPlanes { type empty; }
}
""")
write("0","c","volScalarField","""
dimensions [0 0 0 0 0 0 0];
internalField uniform 1;
boundaryField
{
    yMin { type zeroGradient; }
    yMax { type zeroGradient; }
    xMin { type cyclic; }
    xMax { type cyclic; }
    emptyPlanes { type empty; }
}
""")
print("Generated 40-cell planar Couette: nu=0.05 m2/s, top Ux=1m/s, t=1s")
