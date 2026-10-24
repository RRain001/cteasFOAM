#!/usr/bin/env python3
"""Build M1C-8 fixed physical source/IC with variable true native cell count."""
from pathlib import Path
import math
import re
import sys

CONFIG = {
    "space_lo": ("6x16x6", 576, ".02", 20),
    "baseline": ("8x20x6", 960, ".02", 20),
    "space_hi": ("12x30x8", 2880, ".02", 20),
    "time_lo": ("8x20x6", 960, ".04", 10),
    "time_hi": ("8x20x6", 960, ".01", 40),
}
if len(sys.argv)!=3 or sys.argv[2] not in CONFIG:
    raise SystemExit("usage: make_convergence_case.py CASE space_lo|baseline|space_hi|time_lo|time_hi")
case,mode=Path(sys.argv[1]),sys.argv[2]
grid,n,dt,steps=CONFIG[mode]
if not (case/"0/HT").is_file() or not (case/"0/C").is_file():
    raise SystemExit("FAIL: native mesh geometry/species must exist before M1C-8 setup")
props=case/"constant/transportProperties"
with props.open("a") as f:
    f.write("""
sourceShape gaussian;
sourceCenter [0 1 0 0 0 0 0] (0.38 0.62 0.55);
sourceSigma [0 1 0 0 0 0 0] 0.13;
sourceHT_Bq_s [0 0 -2 0 0 0 0] 60000;
sourceHTO_Bq_s [0 0 -2 0 0 0 0] 12000;
sourceStart [0 0 1 0 0 0 0] 0.075;
sourceEnd [0 0 1 0 0 0 0] 0.225;
""")
for species in ("HT","HTO"):
    p=case/"0"/species
    old=p.read_text()
    for patch in ("xMin","xMax"):
        pattern=rf"({patch}\s*\{{\s*type inletOutlet;\s*inletValue uniform )([-+0-9.eE]+)(;)"
        old,nsub=re.subn(pattern,r"\g<1>0\g<3>",old)
        if nsub!=1: raise SystemExit(f"FAIL: missing zero-Bq inlet {species}/{patch}")
    p.write_text(old)
ctrl=case/"system/controlDict"
s=ctrl.read_text()
if s.count("deltaT 0.02;")!=1 or s.count("writeInterval 20;")!=1:
    raise SystemExit("FAIL: unsupported frozen time schedule template")
ctrl.write_text(s.replace("deltaT 0.02;",f"deltaT {dt};").replace("writeInterval 20;",f"writeInterval {steps};"))
s=(case/"0/C").read_text()
m=re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",s,re.S)
if not m or int(m.group(1))!=n:raise SystemExit("FAIL: native centres do not match grid")
num=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
centres=[tuple(map(float,t)) for t in re.findall(r"\(\s*("+num+r")\s+("+num+r")\s+("+num+r")\s*\)",m.group(2))]
if len(centres)!=n:raise SystemExit("FAIL: incomplete native centres")
u0=[(0.05*math.sin(math.pi*y)*math.sin(math.pi*z),0.,0.) for x,y,z in centres]
p=case/"0/U"
s=p.read_text();old="internalField uniform (0 0 0);"
if s.count(old)!=1:raise SystemExit("FAIL: invalid native initial velocity template")
native="internalField nonuniform List<vector>\n"+str(n)+"\n(\n"+''.join(
    f"({a:.17g} {b:.17g} {c:.17g})\n" for a,b,c in u0)+");"
p.write_text(s.replace(old,native))
print(f"M1C-8 generated {mode} native_cells={n} grid={grid} dt={dt} steps={steps} source_total_HT=9000Bq HTO=1800Bq")
