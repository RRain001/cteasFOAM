#!/usr/bin/env python3
"""Unforced, identical-input periodic PISO benchmark for upstream icoFoam."""
from pathlib import Path
import sys

def replace_once(path, old, new):
    s=path.read_text()
    if s.count(old)!=1:
        raise SystemExit(f"FAIL: expected one {old!r} in {path}")
    path.write_text(s.replace(old,new))

if len(sys.argv)!=2:
    raise SystemExit("usage: make_warped_periodic_upstream_piso.py CASE")
case=Path(sys.argv[1])
props=case/"constant/transportProperties"
replace_once(props,"drivingAcceleration [0 1 -2 0 0 0 0] (0.5 0 0);",
             "drivingAcceleration [0 1 -2 0 0 0 0] (0 0 0);")
replace_once(case/"0/U","internalField uniform (0 0 0);",
             "internalField uniform (0.25 0 0);")
sol=case/"system/fvSolution"
s=sol.read_text()
if "PISO" in s or "pFinal" in s: raise SystemExit("FAIL: already configured reference PISO")
i=s.index("solvers"); b=s.index("{",i)
depth=0
for j in range(b,len(s)):
    if s[j]=="{": depth+=1
    elif s[j]=="}":
        depth-=1
        if depth==0: break
else: raise SystemExit("FAIL: malformed fvSolution")
s=s[:j] + "pFinal { $p; relTol 0; }\n" + s[j:]
s += """
PISO
{
    momentumPredictor yes;
    nCorrectors 3;
    nNonOrthogonalCorrectors 2;
    pRefCell 0;
    pRefValue 0;
}
"""
sol.write_text(s)
print("M1B-1g-12 unforced PISO initial U=(0.25,0,0), 3x3 pressure reference")
