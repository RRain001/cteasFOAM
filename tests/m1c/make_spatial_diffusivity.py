#!/usr/bin/env python3
"""Real user-specified OpenFOAM spatial HT/HTO diffusion scalar field inputs."""
import math, re, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from make_3d_analytic import NATIVE_VECTOR, XYZ

def write(case,n,mode,open_case=False):
    text=(case/"0/C").read_text()
    m=NATIVE_VECTOR.search(text)
    if not m or int(m.group(1))!=n:raise SystemExit("FAIL: M1C22 genuine native mesh cell centres absent")
    centres=[tuple(map(float,x)) for x in XYZ.findall(m.group(2))]
    if len(centres)!=n:raise SystemExit("FAIL: M1C22 malformed genuine native mesh centres")
    if mode not in ("constant","variable","negative"):raise SystemExit("invalid M1C22 test mode")
    for name,reference in (("DHTSpatial",0.01),("DHTOSpatial",0.02)):
        vals=[]
        for x,y,z in centres:
            if mode=="constant":
                val=reference
            else:
                val=reference*(1+0.4*math.sin(2*math.pi*x)*math.cos(2*math.pi*z))
            vals.append(val)
        if mode=="negative":vals[0]=-reference
        if open_case:
            patch="\n".join(f" {p} {{ type zeroGradient; }}" for p in
                            ("xMin","xMax","yMin","yMax","zMin","zMax"))
        else:
            patch="\n".join(f" {p} {{ type cyclic; }}" for p in
                            ("xMin","xMax","yMin","yMax","zMin","zMax"))
        values="\n".join(f" {v:.17g}" for v in vals)
        contents=(f'FoamFile {{ version 2.0; format ascii; class volScalarField; location "0"; object {name}; }}\n'
             f'dimensions [0 2 -1 0 0 0 0];\n'
             f'internalField nonuniform List<scalar>\n{n}\n(\n{values}\n);\n'
             f'boundaryField {{\n{patch}\n}}\n')
        (case/"0"/name).write_text(contents)
        print(f"M1C22_INPUT species={name} mode={mode} min={min(vals):.17g} max={max(vals):.17g}")

if __name__=="__main__":
    if len(sys.argv) not in (4,5):raise SystemExit("usage CASE NCELLS constant|variable|negative [open]")
    write(Path(sys.argv[1]),int(sys.argv[2]),sys.argv[3],len(sys.argv)==5 and sys.argv[4]=="open")
