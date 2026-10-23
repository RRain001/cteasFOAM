#!/usr/bin/env python3
"""M1B-1g-10: periodic, prescribed nonzero-U analytic scalar experiment."""
from pathlib import Path
import math
import re
import sys

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VECTOR = re.compile(r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)")
CENTRES = re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", re.S)


def replace_once(path, before, after):
    content = path.read_text()
    if content.count(before) != 1:
        raise SystemExit(f"FAIL: expected one occurrence of {before!r} in {path}")
    path.write_text(content.replace(before, after))


def replace_body(path, marker, body):
    content = path.read_text()
    i = content.find(marker)
    if i < 0:
        raise SystemExit(f"FAIL: missing {marker} in {path}")
    path.write_text(content[:i] + body)


def prepare(case):
    mesh = case / "system/blockMeshDict"
    replace_once(mesh, "xMin { type patch; faces ((0 4 7 3)); }",
                 "xMin { type cyclic; neighbourPatch xMax; transform translational; separationVector (1 0 0); faces ((0 4 7 3)); }")
    replace_once(mesh, "xMax { type patch; faces ((1 2 6 5)); }",
                 "xMax { type cyclic; neighbourPatch xMin; transform translational; separationVector (-1 0 0); faces ((1 2 6 5)); }")
    replace_body(case / "constant/transportProperties", "D [0 2 -1 0 0 0 0]",
                 "D [0 2 -1 0 0 0 0] 0.01;\nmomentumMode prescribed;\nauditPassiveScalarFlux true;\n")
    replace_once(case / "system/controlDict", "deltaT 0.02;", "deltaT 0.01;")
    replace_once(case / "system/controlDict", "writeInterval 20;", "writeInterval 40;")
    replace_once(case / "system/fvSchemes", "gradSchemes { default Gauss linear; }",
                 "gradSchemes { default leastSquares; }")
    replace_once(case / "system/fvSchemes", "fluxRequired { default no; p; }",
                 "fluxRequired { default no; p; c; }")
    patches_u = """boundaryField
{
 xMin { type cyclic; }
 xMax { type cyclic; }
 yMin { type fixedValue; value uniform (0.25 0 0); }
 yMax { type fixedValue; value uniform (0.25 0 0); }
 zMin { type fixedValue; value uniform (0.25 0 0); }
 zMax { type fixedValue; value uniform (0.25 0 0); }
}
"""
    patches_c = """boundaryField
{
 xMin { type cyclic; }
 xMax { type cyclic; }
 yMin { type zeroGradient; }
 yMax { type zeroGradient; }
 zMin { type zeroGradient; }
 zMax { type zeroGradient; }
}
"""
    replace_body(case / "0/U", "internalField", "internalField uniform (0.25 0 0);\n" + patches_u)
    replace_body(case / "0/c", "internalField", "internalField uniform 0.5;\n" + patches_c)
    (case / "0/p").unlink()
    print("M1B-1g-10 configured U=(0.25,0,0) D=0.01 dt=0.01 periodic x")


def seed(case, n):
    if n not in (960, 3240):
        raise SystemExit("FAIL: unexpected cell count")
    match = CENTRES.search((case / "0/C").read_text())
    if not match or int(match.group(1)) != n:
        raise SystemExit("FAIL: genuine OpenFOAM cell centres missing")
    centres = [tuple(map(float, xyz)) for xyz in VECTOR.findall(match.group(2))]
    if len(centres) != n or not all(0 < v < 1 for xyz in centres for v in xyz):
        raise SystemExit("FAIL: invalid cell centres")
    vals = [0.5 + 0.25*math.cos(2*math.pi*x) for x,y,z in centres]
    data = "internalField nonuniform List<scalar>\n" + str(n) + "\n(\n"
    data += "\n".join(f"{v:.17g}" for v in vals) + "\n);\n"
    c = case / "0/c"
    source = c.read_text()
    i, j = source.find("internalField"), source.find("boundaryField")
    if i < 0 or j <= i:
        raise SystemExit("FAIL: malformed scalar field")
    c.write_text(source[:i] + data + source[j:])
    print(f"M1B-1g-10 initialized {n} native centres")


if len(sys.argv) < 3 or sys.argv[1] not in ("prepare", "seed"):
    raise SystemExit("usage: make_warped_periodic_advection_diffusion.py prepare CASE | seed CASE N")
case = Path(sys.argv[2])
if sys.argv[1] == "prepare" and len(sys.argv) == 3:
    prepare(case)
elif sys.argv[1] == "seed" and len(sys.argv) == 4:
    seed(case, int(sys.argv[3]))
else:
    raise SystemExit("FAIL: wrong generator arguments")
