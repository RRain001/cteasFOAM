#!/usr/bin/env python3
"""Independent, real 3-D rectangular-duct pressure benchmark case builder.

x=0/1: non-periodic physical fixed static-kinematic-pressure patches;
y=0/1, z=0/1: no-slip walls. The interior-only x warp preserves ALL
six physical boundary surfaces, including inlet/outlet face geometry.
"""
from pathlib import Path
import math
import re
import sys

NX, NY, NZ = 8, 20, 6
AMP = 0.04
NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VEC = re.compile(r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)")
CENTRES = re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", re.S)


def write(case, directory, name, cls, body):
    header = (
        f'FoamFile\n{{\n version 2.0;\n format ascii;\n'
        f' class {cls};\n location "{directory}";\n object {name};\n}}\n'
    )
    (case / directory / name).write_text(header + body)


def make(case, mode):
    if mode not in ("flat", "warped", "zero", "reverse"):
        raise SystemExit("FAIL: expected flat, warped, zero or reverse")
    for dirname in ("0", "constant", "system"):
        (case / dirname).mkdir(parents=True, exist_ok=True)
    dp = 0.0 if mode == "zero" else (-0.1 if mode == "reverse" else 0.1)
    write(case, "system", "blockMeshDict", "dictionary", f"""
scale 1;
vertices
(
 (0 0 0) (1 0 0) (1 1 0) (0 1 0)
 (0 0 1) (1 0 1) (1 1 1) (0 1 1)
);
blocks (hex (0 1 2 3 4 5 6 7) ({NX} {NY} {NZ}) simpleGrading (1 1 1));
edges ();
boundary
(
 xMin {{ type patch; faces ((0 4 7 3)); }}
 xMax {{ type patch; faces ((1 2 6 5)); }}
 yMin {{ type wall; faces ((0 1 5 4)); }}
 yMax {{ type wall; faces ((3 7 6 2)); }}
 zMin {{ type wall; faces ((0 3 2 1)); }}
 zMax {{ type wall; faces ((4 5 6 7)); }}
);
mergePatchPairs ();
""")
    write(case, "system", "controlDict", "dictionary", """
application cteasAtmosphereFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 0.4;
deltaT 0.02;
writeControl timeStep;
writeInterval 20;
writeFormat ascii;
writePrecision 15;
runTimeModifiable false;
""")
    write(case, "system", "fvSchemes", "dictionary", """
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
    write(case, "system", "fvSolution", "dictionary", """
solvers
{
 U { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }
 p { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }
 c { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }
}
""")
    write(case, "constant", "transportProperties", "dictionary", """
D [0 2 -1 0 0 0 0] 0;
nu [0 2 -1 0 0 0 0] 0.05;
momentumMode pisoLaminar;
nPressureCorrectors 3;
nNonOrthogonalCorrectors 2;
auditOpenPressureFlux true;
""")
    walls = (
        " yMin { type fixedValue; value uniform (0 0 0); }\n"
        " yMax { type fixedValue; value uniform (0 0 0); }\n"
        " zMin { type fixedValue; value uniform (0 0 0); }\n"
        " zMax { type fixedValue; value uniform (0 0 0); }\n"
    )
    write(case, "0", "U", "volVectorField",
          "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0 0 0);\n"
          "boundaryField\n{\n"
          " xMin { type pressureInletOutletVelocity; value uniform (0 0 0); }\n"
          " xMax { type pressureInletOutletVelocity; value uniform (0 0 0); }\n"
          + walls + "}\n")
    write(case, "0", "p", "volScalarField",
          "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\n"
          "boundaryField\n{\n"
          f" xMin {{ type fixedValue; value uniform {dp:.17g}; }}\n"
          " xMax { type fixedValue; value uniform 0; }\n"
          + "".join(f" {patch} {{ type fixedFluxPressure; value uniform 0; }}\n"
                    for patch in ("yMin", "yMax", "zMin", "zMax"))
          + "}\n")
    write(case, "0", "c", "volScalarField",
          "dimensions [0 0 0 0 0 0 0];\ninternalField uniform 1;\n"
          "boundaryField\n{\n"
          " xMin { type fixedValue; value uniform 1; }\n"
          " xMax { type zeroGradient; }\n"
          + "".join(f" {patch} {{ type zeroGradient; }}\n"
                    for patch in ("yMin", "yMax", "zMin", "zMax"))
          + "}\n")
    print(f"M1B-1g generated mode={mode} cells={NX*NY*NZ} dp={dp} real3D")


def warp(case):
    path = case / "constant" / "polyMesh" / "points"
    source = path.read_text()
    m = re.search(r"(?m)^\s*(\d+)\s*\n\s*\(\s*\n", source)
    if m is None:
        raise SystemExit("FAIL: missing actual OpenFOAM counted points")
    tail = source[m.end():]
    ending = re.search(r"(?m)^\s*\)\s*$", tail)
    if ending is None:
        raise SystemExit("FAIL: malformed OpenFOAM points terminator")
    points = [tuple(map(float, v)) for v in VEC.findall(tail[:ending.start()])]
    if len(points) != int(m.group(1)) or len(points) != (NX+1)*(NY+1)*(NZ+1):
        raise SystemExit("FAIL: wrong real 3-D mesh point count")
    changed, updated = 0, []
    for x, y, z in points:
        boundary = any(abs(v) < 1e-12 or abs(v-1) < 1e-12 for v in (x, y, z))
        dx = 0.0 if boundary else AMP*math.sin(math.pi*x)*math.sin(math.pi*y)*math.sin(math.pi*z)
        changed += abs(dx) > 1e-12
        new = (x+dx, y, z)
        if not 0 <= new[0] <= 1:
            raise SystemExit("FAIL: warped vertex leaves domain")
        updated.append(new)
    if changed < 100:
        raise SystemExit("FAIL: genuine interior 3-D mesh deformation missing")
    block = "\n".join(f"({x:.17g} {y:.17g} {z:.17g})" for x, y, z in updated)
    path.write_text(source[:m.end()] + block + "\n" + tail[ending.start():])
    print(f"M1B-1g real polyMesh points={len(points)} moved={changed} amplitude={AMP}")


def seed(case, mode):
    field = case / "0" / "C"
    if not field.is_file():
        raise SystemExit("FAIL: must run OpenFOAM writeCellCentres on actual mesh")
    m = CENTRES.search(field.read_text())
    if m is None or int(m.group(1)) != NX*NY*NZ:
        raise SystemExit("FAIL: invalid native cell centres")
    coords = [tuple(map(float, x)) for x in VEC.findall(m.group(2))]
    if len(coords) != NX*NY*NZ or not all(math.isfinite(v) for xyz in coords for v in xyz):
        raise SystemExit("FAIL: incomplete/nonfinite native centres")
    if mode == "warped" and max(abs(x - ((i % NX)+.5)/NX) for i, (x,y,z) in enumerate(coords)) < .01:
        raise SystemExit("FAIL: warped mesh centres not physically displaced")
    dp = 0.0 if mode == "zero" else (-.1 if mode == "reverse" else .1)
    vals = [dp*(1-x) for x, y, z in coords]
    path = case / "0" / "p"
    txt = path.read_text()
    first, last = txt.find("internalField"), txt.find("boundaryField")
    if first < 0 or last < first:
        raise SystemExit("FAIL: malformed initial p field")
    seedtxt = (f"internalField nonuniform List<scalar>\n{len(vals)}\n(\n"
               + "\n".join(f"{v:.17g}" for v in vals) + "\n);\n")
    path.write_text(txt[:first] + seedtxt + txt[last:])
    print(f"M1B-1g seeded p(x) from actual OpenFOAM C cells={len(vals)} mode={mode}")


if len(sys.argv) not in (4, 5):
    raise SystemExit("usage: make_warped_open_pressure_case.py make|warp|seed CASE flat|warped|zero|reverse [NXxNYxNZ]")
action, directory, mode = sys.argv[1:4]
if len(sys.argv) == 5:
    match = re.fullmatch(r"([1-9][0-9]*)x([1-9][0-9]*)x([1-9][0-9]*)", sys.argv[4])
    if not match:
        raise SystemExit("FAIL: invalid three-dimensional resolution specification")
    NX, NY, NZ = map(int, match.groups())
    if min(NX, NY, NZ) < 6 or max(NX, NY, NZ) > 80:
        raise SystemExit("FAIL: resolution is outside bounded canonical mesh range")
target = Path(directory)
if action == "make":
    make(target, mode)
elif action == "warp" and mode in ("warped", "zero", "reverse"):
    warp(target)
elif action == "seed":
    seed(target, mode)
else:
    raise SystemExit("FAIL: invalid M1B-1g generator action/mode")
