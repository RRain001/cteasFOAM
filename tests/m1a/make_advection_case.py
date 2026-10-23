#!/usr/bin/env python3
"""Generate reproducible periodic OpenFOAM transport cases, n=40 and n=80."""
import math
from pathlib import Path
import sys

if len(sys.argv) != 3:
    raise SystemExit("usage: make_advection_case.py OUTPUT N")
case, n = Path(sys.argv[1]), int(sys.argv[2])
if n not in (40, 80):
    raise SystemExit("required refinement levels: 40 and 80")
for folder in ("0", "constant", "system"):
    (case / folder).mkdir(parents=True, exist_ok=True)

def write(folder, name, cls, body):
    header = (f'FoamFile\n{{\n    version 2.0;\n    format ascii;\n'
              f'    class {cls};\n    location "{folder}";\n'
              f'    object {name};\n}}\n')
    (case / folder / name).write_text(header + body)

write("system", "blockMeshDict", "dictionary", f"""
scale 1;
vertices
(
    (0 0 0) (1 0 0) (1 1 0) (0 1 0)
    (0 0 1) (1 0 1) (1 1 1) (0 1 1)
);
blocks (hex (0 1 2 3 4 5 6 7) ({n} 1 1) simpleGrading (1 1 1));
edges ();
boundary
(
    xMin {{ type cyclic; neighbourPatch xMax; transform translational;
           faces ((0 4 7 3)); }}
    xMax {{ type cyclic; neighbourPatch xMin; transform translational;
           faces ((1 2 6 5)); }}
    emptyPlanes {{ type empty;
                   faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }}
);
mergePatchPairs ();
""")
write("system", "controlDict", "dictionary", """
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
write("system", "fvSchemes", "dictionary", """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,c) Gauss upwind; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
write("system", "fvSolution", "dictionary", """
solvers
{
    c
    {
        solver PBiCGStab;
        preconditioner DILU;
        tolerance 1e-13;
        relTol 0;
    }
}
""")
write("constant", "transportProperties", "dictionary",
      "D [0 2 -1 0 0 0 0] 0;\n")
write("0", "U", "volVectorField", """
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0.25 0 0);
boundaryField
{
    xMin { type cyclic; }
    xMax { type cyclic; }
    emptyPlanes { type empty; }
}
""")
values = "\n".join(
    f"    {1 + 0.2 * math.cos(2 * math.pi * (i + 0.5) / n):.17g}"
    for i in range(n)
)
write("0", "c", "volScalarField", f"""
dimensions [0 0 0 0 0 0 0];
internalField nonuniform List<scalar>
{n}
(
{values}
);
boundaryField
{{
    xMin {{ type cyclic; }}
    xMax {{ type cyclic; }}
    emptyPlanes {{ type empty; }}
}}
""")
print(f"Generated cyclic benchmark n={n}, U=0.25, D=0, t=1")
