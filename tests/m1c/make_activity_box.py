#!/usr/bin/env python3
"""Generate a closed three-dimensional OpenFOAM HT/HTO activity box."""
import sys
from pathlib import Path

if len(sys.argv) != 5:
    raise SystemExit("usage: make_activity_box.py CASE kOx lambdaT valid|negative-k|negative-lambda|bad-dims")
case, k, lam, mode = Path(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
if mode not in {"valid", "negative-k", "negative-lambda", "bad-dims"}:
    raise SystemExit("FAIL: invalid mode")
if mode == "negative-k":
    k = -abs(k)
if mode == "negative-lambda":
    lam = -abs(lam)
for name in ("0", "constant", "system"):
    (case / name).mkdir(parents=True, exist_ok=True)

def write(folder, name, cls, body):
    header = (f'FoamFile\n{{\n version 2.0;\n format ascii;\n'
              f' class {cls};\n location "{folder}";\n object {name};\n}}\n')
    (case / folder / name).write_text(header + body)

write("system", "blockMeshDict", "dictionary", """
scale 1;
vertices
(
 (0 0 0) (1 0 0) (1 1 0) (0 1 0)
 (0 0 1) (1 0 1) (1 1 1) (0 1 1)
);
blocks (hex (0 1 2 3 4 5 6 7) (4 3 2) simpleGrading (1 1 1));
edges ();
boundary
(
 xMin { type wall; faces ((0 4 7 3)); }
 xMax { type wall; faces ((1 2 6 5)); }
 yMin { type wall; faces ((0 1 5 4)); }
 yMax { type wall; faces ((3 7 6 2)); }
 zMin { type wall; faces ((0 3 2 1)); }
 zMax { type wall; faces ((4 5 6 7)); }
);
mergePatchPairs ();
""")
write("system", "controlDict", "dictionary", """
application cteasAtmosphereFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 1;
deltaT 0.05;
writeControl timeStep;
writeInterval 20;
writeFormat ascii;
writePrecision 17;
runTimeModifiable false;
""")
write("system", "fvSchemes", "dictionary", """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes
{
 default none;
 div(phi,c) Gauss upwind;
 div(phi,HT) Gauss upwind;
 div(phi,HTO) Gauss upwind;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
write("system", "fvSolution", "dictionary", """
solvers
{
 c { solver PBiCGStab; preconditioner DILU; tolerance 1e-13; relTol 0; }
 HT { solver PBiCGStab; preconditioner DILU; tolerance 1e-13; relTol 0; }
 HTO { solver PBiCGStab; preconditioner DILU; tolerance 1e-13; relTol 0; }
}
""")
write("constant", "transportProperties", "dictionary", f"""
D [0 2 -1 0 0 0 0] 0;
momentumMode prescribed;
activityMode htHto;
kOx [0 0 -1 0 0 0 0] {k:.17g};
lambdaT [0 0 -1 0 0 0 0] {lam:.17g};
DHT [0 2 -1 0 0 0 0] 0;
DHTO [0 2 -1 0 0 0 0] 0;
auditClosedActivity true;
""")
patches = ("xMin", "xMax", "yMin", "yMax", "zMin", "zMax")
u_boundary = "".join(f" {p} {{ type fixedValue; value uniform (0 0 0); }}\n" for p in patches)
scalar_boundary = "".join(f" {p} {{ type zeroGradient; }}\n" for p in patches)
write("0", "U", "volVectorField",
      "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0 0 0);\n"
      f"boundaryField\n{{\n{u_boundary}}}\n")
write("0", "c", "volScalarField",
      "dimensions [0 0 0 0 0 0 0];\ninternalField uniform 0;\n"
      f"boundaryField\n{{\n{scalar_boundary}}}\n")
for species, value in (("HT", 1e6), ("HTO", 2e5)):
    dims = "[0 0 0 0 0 0 0]" if mode == "bad-dims" and species == "HT" else "[0 -3 -1 0 0 0 0]"
    write("0", species, "volScalarField",
          f"dimensions {dims};\ninternalField uniform {value:.17g};\n"
          f"boundaryField\n{{\n{scalar_boundary}}}\n")
print(f"M1C activity box k={k} lambda={lam} mode={mode} 3D=4x3x2")
