#!/usr/bin/env python3
"""M1C-7 pressure-open PISO warped 3D with internal Gaussian timed HT/HTO."""
from pathlib import Path
import math
import re
import sys
if len(sys.argv) != 3 or sys.argv[2] not in ("release", "zero_rate"):
    raise SystemExit("usage: make_open_gaussian_release.py CASE release|zero_rate")
case, mode = Path(sys.argv[1]), sys.argv[2]
if not (case/"0/HT").exists():
    raise SystemExit("FAIL: first build real 960-cell pressure case with make_open_activity")
props=case/"constant/transportProperties"
rate_ht, rate_hto = (60000, 12000) if mode=="release" else (0,0)
with props.open("a") as f:
    f.write("sourceShape gaussian;\n")
    f.write("sourceCenter [0 1 0 0 0 0 0] (0.38 0.62 0.55);\n")
    f.write("sourceSigma [0 1 0 0 0 0 0] 0.13;\n")
    f.write(f"sourceHT_Bq_s [0 0 -2 0 0 0 0] {rate_ht};\n")
    f.write(f"sourceHTO_Bq_s [0 0 -2 0 0 0 0] {rate_hto};\n")
    f.write("sourceStart [0 0 1 0 0 0 0] 0.075;\n")
    f.write("sourceEnd [0 0 1 0 0 0 0] 0.225;\n")
for species in ("HT","HTO"):
    path=case/"0"/species
    s=path.read_text()
    for patch in ("xMin","xMax"):
        pattern=rf"({patch}\s*\{{\s*type inletOutlet;\s*inletValue uniform )([-+0-9.eE]+)(;)"
        s,n=re.subn(pattern,r"\g<1>0\g<3>",s)
        if n!=1:
            raise SystemExit(f"FAIL: missing native {species} {patch} inletOutlet")
    path.write_text(s)
# Independently reproducible, streamwise-divergence-free startup profile,
# on genuine warped OpenFOAM cell centres (and zero side-wall BCs).
# This avoids an incorrect demand for fully developed flow at t=0 from rest.
cell_file=(case/"0/C").read_text()
cm=re.search(r"internalField\s+nonuniform\s+List<vector>\s+960\s*\((.*?)\)\s*;",cell_file,re.S)
if not cm: raise SystemExit("FAIL: missing 960 native OpenFOAM cell centres")
num=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
centres=[tuple(map(float,row)) for row in re.findall(
    r"\(\s*("+num+r")\s+("+num+r")\s+("+num+r")\s*\)",cm.group(1))]
if len(centres)!=960:raise SystemExit("FAIL: cannot seed 960 local U cells")
u0=[(0.05*math.sin(math.pi*y)*math.sin(math.pi*z),0,0) for x,y,z in centres]
u_path=case/"0/U"
u_file=u_path.read_text()
old="internalField uniform (0 0 0);"
if u_file.count(old)!=1:raise SystemExit("FAIL: unknown PISO initial U template")
initial="internalField nonuniform List<vector>\n960\n(\n"+"".join(
    f"({x:.17g} {y:.17g} {z:.17g})\n" for x,y,z in u0)+");"
u_path.write_text(u_file.replace(old,initial))
print(f"M1C-7 prepared mode={mode} native pressure-open 3D local source HT={rate_ht} HTO={rate_hto} Bq/s")
