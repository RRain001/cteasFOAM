#!/usr/bin/env python3
"""M1C-3 exact-native warped open pressure duct HT/HTO activity inlet cases."""
from pathlib import Path
import re
import sys

if len(sys.argv) != 3 or sys.argv[2] not in ("warped","reverse","zero"):
    raise SystemExit("usage: make_open_activity.py CASE warped|reverse|zero")
case, mode = Path(sys.argv[1]), sys.argv[2]
active = "xMax" if mode == "reverse" else "xMin"

def change(path, before, after):
    data=path.read_text()
    if data.count(before)!=1: raise SystemExit(f"FAIL: unexpected template {path}: {before!r}")
    path.write_text(data.replace(before,after))

props=case/"constant/transportProperties"
with props.open("a") as f:
    f.write("""activityMode htHto;
kOx [0 0 -1 0 0 0 0] 0.4;
lambdaT [0 0 -1 0 0 0 0] 0.2;
DHT [0 2 -1 0 0 0 0] 0.01;
DHTO [0 2 -1 0 0 0 0] 0.02;
auditActivityFlux true;
""")
scheme=case/"system/fvSchemes"
change(scheme,"div(phi,c) Gauss upwind;",
       "div(phi,c) Gauss upwind;\n div(phi,HT) Gauss upwind;\n div(phi,HTO) Gauss upwind;")
change(scheme,"fluxRequired { default no; p; }",
       "fluxRequired { default no; p; HT; HTO; }")
solution=case/"system/fvSolution"
data=solution.read_text()
line=next((line for line in data.splitlines(True) if line.startswith(" c { solver PBiCGStab;")),None)
if line is None: raise SystemExit("FAIL: required native c equation solver")
change(solution,line,line+line.replace(" c {"," HT {")+line.replace(" c {"," HTO {"))
for name,inlet in (("HT",1000000),("HTO",200000)):
    field=(case/"0/c").read_text()
    field=field.replace("object c;","object "+name+";")
    field=field.replace("dimensions [0 0 0 0 0 0 0];","dimensions [0 -3 -1 0 0 0 0];")
    field=field.replace("internalField uniform 1;","internalField uniform 0;")
    for patch,val in (("xMin",inlet if active=="xMin" else 0),
                      ("xMax",inlet if active=="xMax" else 0)):
        field,n=re.subn(rf"\b{patch}\s*\{{[^{{}}]*\}}",
                      f"{patch} {{ type inletOutlet; inletValue uniform {val}; value uniform 0; }}",
                      field,count=1)
        if n!=1: raise SystemExit(f"FAIL: missing {patch} species patch")
    if (field.count("inletOutlet")!=2 or "internalField uniform 0;" not in field
        or "dimensions [0 -3 -1 0 0 0 0];" not in field):
        raise SystemExit("FAIL: incomplete species initialization")
    (case/"0"/name).write_text(field)
print(f"M1C-3 prepared mode={mode} active={active} HT=1000000 HTO=200000 Bq/m3 DHT=0.01 DHTO=0.02")
