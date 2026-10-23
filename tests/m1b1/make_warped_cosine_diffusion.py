#!/usr/bin/env python3
"""M1B-1g-9: seed a real warped mesh with the transient Neumann eigenmode."""
from pathlib import Path
import math
import re
import sys

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VECTOR = re.compile(r"\(\s*(" + NUMBER + r")\s+(" + NUMBER + r")\s+(" + NUMBER + r")\s*\)")
CENTRES = re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", re.S)

if len(sys.argv) != 3:
    raise SystemExit("usage: make_warped_cosine_diffusion.py CASE CELL_COUNT")
case, n = Path(sys.argv[1]), int(sys.argv[2])
if n not in (960, 3240):
    raise SystemExit("FAIL: unexpected canonical resolution")
m = CENTRES.search((case / "0/C").read_text())
if not m or int(m.group(1)) != n:
    raise SystemExit("FAIL: genuine OpenFOAM cell centres missing")
centres = [tuple(map(float, v)) for v in VECTOR.findall(m.group(2))]
if len(centres) != n or not all(0 < x < 1 and 0 < y < 1 and 0 < z < 1 for x,y,z in centres):
    raise SystemExit("FAIL: malformed physical cell centres")
values = [0.5 + 0.25 * math.cos(math.pi * x) for x,y,z in centres]
field = case / "0/c"
source = field.read_text()
a, b = source.find("internalField"), source.find("boundaryField")
if a < 0 or b <= a:
    raise SystemExit("FAIL: malformed native scalar field")
replacement = ("internalField nonuniform List<scalar>\n" + str(n) + "\n(\n"
               + "\n".join(f"{v:.17g}" for v in values) + "\n);\n")
source = source[:a] + replacement + source[b:]
for patch in ("xMin", "xMax", "yMin", "yMax", "zMin", "zMax"):
    source, hits = re.subn(rf"\b{patch}\s*\{{[^{{}}]*\}}",
                            f"{patch} {{ type zeroGradient; }}", source, count=1)
    if hits != 1:
        raise SystemExit(f"FAIL: missing native patch {patch}")
field.write_text(source)
properties = case / "constant/transportProperties"
data = properties.read_text()
old = "D [0 2 -1 0 0 0 0] 0;"
if data.count(old) != 1:
    raise SystemExit("FAIL: native D coefficient not unique")
properties.write_text(data.replace(old, "D [0 2 -1 0 0 0 0] 0.01;") + "auditPassiveScalarFlux true;\n")
schemes = case / "system/fvSchemes"
data = schemes.read_text()
for old, new in (("fluxRequired { default no; p; }", "fluxRequired { default no; p; c; }"),
                 ("gradSchemes { default Gauss linear; }", "gradSchemes { default leastSquares; }")):
    if data.count(old) != 1:
        raise SystemExit("FAIL: native discretization unexpectedly changed")
    data = data.replace(old, new)
schemes.write_text(data)
print(f"M1B-1g-9 native cosine Neumann mode cells={n} D=0.01")
