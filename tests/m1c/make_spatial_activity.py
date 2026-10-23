#!/usr/bin/env python3
"""M1C-1: augment native warped periodic mesh case with HT/HTO fields."""
import math
import re
import sys
from pathlib import Path

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VECTOR = re.compile(r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)")
CENTRES = re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", re.S)

def replace_once(path, before, after):
    text = path.read_text()
    if text.count(before) != 1:
        raise SystemExit(f"FAIL: expected exactly one {before!r} in {path}")
    path.write_text(text.replace(before, after))

def prepare(case):
    props = case / "constant/transportProperties"
    replace_once(props, "auditPassiveScalarFlux true;",
                 """auditPassiveScalarFlux true;
activityMode htHto;
kOx [0 0 -1 0 0 0 0] 0.4;
lambdaT [0 0 -1 0 0 0 0] 0.2;
DHT [0 2 -1 0 0 0 0] 0.01;
DHTO [0 2 -1 0 0 0 0] 0.02;
auditActivityFlux true;""")
    scheme = case / "system/fvSchemes"
    replace_once(scheme, "div(phi,c) Gauss upwind;",
                 "div(phi,c) Gauss upwind;\n div(phi,HT) Gauss upwind;\n div(phi,HTO) Gauss upwind;")
    replace_once(scheme, "fluxRequired { default no; p; c; }",
                 "fluxRequired { default no; p; c; HT; HTO; }")
    solution = case / "system/fvSolution"
    text = solution.read_text()
    if text.count(" c { solver PBiCGStab;") != 1:
        raise SystemExit("FAIL: original c solver entry missing")
    c_line = next(line for line in text.splitlines(True) if line.startswith(" c { solver PBiCGStab;"))
    replace_once(solution, c_line, c_line +
                 c_line.replace(" c {", " HT {") +
                 c_line.replace(" c {", " HTO {"))
    template = (case / "0/c").read_text()
    if "dimensions [0 0 0 0 0 0 0];" not in template:
        raise SystemExit("FAIL: original scalar dimensions")
    for species, mean in (("HT", 1e6), ("HTO", 2e5)):
        field = template.replace("object c;", f"object {species};")
        field = field.replace("dimensions [0 0 0 0 0 0 0];",
                              "dimensions [0 -3 -1 0 0 0 0];")
        field = field.replace("internalField uniform 0.5;",
                              f"internalField uniform {mean:.17g};")
        if "internalField uniform" not in field:
            raise SystemExit("FAIL: original scalar internalField missing")
        (case / "0" / species).write_text(field)
    print("M1C-1 prepared k=0.4 lambda=0.2 DHT=0.01 DHTO=0.02")

def seed(case, n):
    m = CENTRES.search((case / "0/C").read_text())
    if not m or int(m.group(1)) != n:
        raise SystemExit("FAIL: native centres missing or wrong length")
    centres = [tuple(map(float, v)) for v in VECTOR.findall(m.group(2))]
    if len(centres) != n or not all(0 < v < 1 for row in centres for v in row):
        raise SystemExit("FAIL: invalid native centres")
    for species, mean, amp in (("HT", 1e6, 2e5), ("HTO", 2e5, 5e4)):
        vals = [mean + amp*math.cos(2*math.pi*x) for x,y,z in centres]
        path = case / "0" / species
        data = path.read_text()
        i, j = data.find("internalField"), data.find("boundaryField")
        if i < 0 or j <= i:
            raise SystemExit("FAIL: missing species field sections")
        values = "internalField nonuniform List<scalar>\n" + str(n) + "\n(\n"
        values += "\n".join(f"{v:.17g}" for v in vals) + "\n);\n"
        path.write_text(data[:i] + values + data[j:])
    print(f"M1C-1 initialized HT and HTO on {n} real native centres")

if len(sys.argv) < 3 or sys.argv[1] not in ("prepare", "seed"):
    raise SystemExit("usage: make_spatial_activity.py prepare CASE | seed CASE N")
case = Path(sys.argv[2])
if sys.argv[1] == "prepare" and len(sys.argv) == 3:
    prepare(case)
elif sys.argv[1] == "seed" and len(sys.argv) == 4:
    seed(case, int(sys.argv[3]))
else:
    raise SystemExit("FAIL: wrong arguments")
