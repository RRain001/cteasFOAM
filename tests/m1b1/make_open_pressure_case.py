#!/usr/bin/env python3
"""Create an actual static-pressure inlet/outlet laminar channel, plus zero-dp control."""
from pathlib import Path
import math
import sys

if len(sys.argv)!=3:
    raise SystemExit("usage: make_open_pressure_case.py CASE driven|zero")
case=Path(sys.argv[1])
kind=sys.argv[2]
if kind not in ("driven","zero"):
    raise SystemExit("case type must be driven or zero")
dp=0.1 if kind=="driven" else 0.0
for d in ("0","constant","system"):
    (case/d).mkdir(parents=True,exist_ok=True)
def write(d,name,cls,body):
    h=(f'FoamFile\n{{\n version 2.0;\n format ascii;\n'
       f' class {cls};\n location "{d}";\n object {name};\n}}\n')
    (case/d/name).write_text(h+body)

write("system","blockMeshDict","dictionary","""
scale 1;
vertices
(
 (0 0 0) (1 0 0) (1 1 0) (0 1 0)
 (0 0 1) (1 0 1) (1 1 1) (0 1 1)
);
blocks (hex (0 1 2 3 4 5 6 7) (8 40 1) simpleGrading (1 1 1));
edges ();
boundary
(
 xMin { type patch; faces ((0 4 7 3)); }
 xMax { type patch; faces ((1 2 6 5)); }
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
nu [0 2 -1 0 0 0 0] 0.05;
momentumMode pisoLaminar;
nPressureCorrectors 3;
auditOpenPressureFlux true;
""")
uPatch="""boundaryField
{
 xMin { type pressureInletOutletVelocity; value uniform (0 0 0); }
 xMax { type pressureInletOutletVelocity; value uniform (0 0 0); }
 yMin { type fixedValue; value uniform (0 0 0); }
 yMax { type fixedValue; value uniform (0 0 0); }
 emptyPlanes { type empty; }
}
"""
pPatch=f"""boundaryField
{{
 xMin {{ type fixedValue; value uniform {dp:.17g}; }}
 xMax {{ type fixedValue; value uniform 0; }}
 yMin {{ type fixedFluxPressure; value uniform 0; }}
 yMax {{ type fixedFluxPressure; value uniform 0; }}
 emptyPlanes {{ type empty; }}
}}
"""
cPatch="""boundaryField
{
 xMin { type fixedValue; value uniform 1; }
 xMax { type zeroGradient; }
 yMin { type zeroGradient; }
 yMax { type zeroGradient; }
 emptyPlanes { type empty; }
}
"""
p0=[]
for j in range(40):
    for i in range(8):
        x=(i+.5)/8
        p0.append(f"{dp*(1-x):.17g}")
write("0","U","volVectorField",
      "dimensions [0 1 -1 0 0 0 0];\n"
      "internalField uniform (0 0 0);\n"+uPatch)
write("0","p","volScalarField",
      "dimensions [0 2 -2 0 0 0 0];\n"
      f"internalField nonuniform List<scalar>\n{len(p0)}\n(\n"
      +"\n".join(p0)+"\n);\n"+pPatch)
write("0","c","volScalarField",
      "dimensions [0 0 0 0 0 0 0];\n"
      "internalField uniform 1;\n"+cPatch)
print(f"Generated actual two-ended fixed pressure case={kind} Δp={dp} m²/s² 8x40x1")
