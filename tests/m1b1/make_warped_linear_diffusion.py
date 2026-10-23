#!/usr/bin/env python3
"""M1B-1g-8: independent affine steady diffusion reference on warped native mesh."""
from pathlib import Path
import re
import sys
from check_warped_open_pressure import read_field

if len(sys.argv) != 2:
    raise SystemExit("usage: make_warped_linear_diffusion.py CASE")
case = Path(sys.argv[1])
centres = read_field(case, "0", "C", True)
if len(centres) != 960 or not all(0 < xyz[0] < 1 for xyz in centres):
    raise SystemExit("FAIL: missing real 960-cell native centroids")
c0 = case / "0/c"
src = c0.read_text()
start, end = src.index("internalField"), src.index("boundaryField")
if end <= start:
    raise SystemExit("FAIL: invalid native scalar field")
values = [1 - xyz[0] for xyz in centres]
field = "internalField nonuniform List<scalar>\n960\n(\n" + "\n".join(f"{v:.17g}" for v in values) + "\n);\n"
src = src[:start] + field + src[end:]
for name, val in (("xMin", 1), ("xMax", 0)):
    pattern = rf"\b{name}\s*\{{[^{{}}]*\}}"
    src, n = re.subn(pattern, f"{name} {{ type fixedValue; value uniform {val}; }}", src, count=1)
    if n != 1:
        raise SystemExit(f"FAIL: missing physical {name} patch")
for name in ("yMin", "yMax", "zMin", "zMax"):
    if not re.search(rf"\b{name}\s*\{{\s*type\s+zeroGradient;", src):
        raise SystemExit(f"FAIL: wall {name} not insulating")
c0.write_text(src)
transport = case / "constant/transportProperties"
data = transport.read_text()
needle = "D [0 2 -1 0 0 0 0] 0;"
if data.count(needle) != 1:
    raise SystemExit("FAIL: expected unique zero-D entry")
transport.write_text(data.replace(needle, "D [0 2 -1 0 0 0 0] 0.01;") + "auditPassiveScalarFlux true;\n")
schemes = case / "system/fvSchemes"
text = schemes.read_text()
old = "fluxRequired { default no; p; }"
if text.count(old) != 1:
    raise SystemExit("FAIL: missing original fluxRequired")
schemes.write_text(text.replace(old, "fluxRequired { default no; p; c; }"))
print(f"M1B-1g-8 configured native-centre affine c=1-x on {len(values)} warped cells; D=0.01")
