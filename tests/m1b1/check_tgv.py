#!/usr/bin/env python3
"""Independent analytic 2D Taylor-Green checker, real OpenFOAM output only."""
import math
from pathlib import Path
import re
import sys

SCALAR_RE=re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
VECTOR_RE=re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
CONT_RE=re.compile(
    r"CTEAS_M1B1 max_div_phi=([-+0-9.eE]+) "
    r"boundary_net_flux=([-+0-9.eE]+) kinetic_energy=([-+0-9.eE]+)"
)
def load_field(case,name,n,kind):
    path=Path(case)/"0.1"/name
    if not path.is_file():
        raise SystemExit(f"FAIL: actual OpenFOAM output absent: {path}")
    txt=path.read_text()
    m=(VECTOR_RE if kind=="vector" else SCALAR_RE).search(txt)
    if not m or int(m.group(1))!=n*n:
        raise SystemExit(f"FAIL: malformed full {name} field n={n}")
    if kind=="vector":
        v=[tuple(map(float,w.split())) for w in re.findall(r"\(([^()]*)\)",m.group(2))]
        if len(v)!=n*n or any(len(z)!=3 for z in v):
            raise SystemExit("FAIL: malformed OpenFOAM velocity data")
        flat=[a for sub in v for a in sub]
    else:
        v=list(map(float,m.group(2).split()))
        flat=v
    if len(v)!=n*n or not all(math.isfinite(a) for a in flat):
        raise SystemExit("FAIL: invalid or nonfinite output")
    return v

def check(case,log,n):
    U=load_field(case,"U",n,"vector")
    p=load_field(case,"p",n,"scalar")
    decay=math.exp(-2*(2*math.pi)**2*.01*.1)
    p_decay=decay**2
    exactU=[]
    exactP=[]
    for j in range(n):
        y=(j+.5)/n
        for i in range(n):
            x=(i+.5)/n
            exactU.append((math.sin(2*math.pi*x)*math.cos(2*math.pi*y)*decay,
                 -math.cos(2*math.pi*x)*math.sin(2*math.pi*y)*decay,0))
            exactP.append(+.25*(math.cos(4*math.pi*x)+math.cos(4*math.pi*y))*p_decay)
    uerr=max(abs(a-b) for v,w in zip(U,exactU) for a,b in zip(v,w))
    pmean=sum(p)/(n*n)
    pref=sum(exactP)/(n*n)
    perr=max(abs((a-pmean)-(b-pref)) for a,b in zip(p,exactP))
    energy=sum(.5*(v[0]**2+v[1]**2+v[2]**2) for v in U)/(n*n)
    hist=[tuple(map(float,z)) for z in CONT_RE.findall(Path(log).read_text())]
    if len(hist)!=20:
        raise SystemExit(f"FAIL: expected 20 actual PISO steps, got {len(hist)}")
    max_div=max(abs(z[0]) for z in hist)
    max_net=max(abs(z[1]) for z in hist)
    if not all(math.isfinite(z) for row in hist for z in row):
        raise SystemExit("FAIL: nonfinite continuity or energy logs")
    print(f"M1B-1a TGV n={n} velocity_error={uerr:.12g} "
          f"pressure_error_gauge_free={perr:.12g} "
          f"max_div_phi={max_div:.12g} max_boundary_flux={max_net:.12g} "
          f"final_energy={energy:.12g}")
    if max_div>1e-7 or max_net>1e-7:
        raise SystemExit("FAIL: corrected face flux violates continuity")
    if not 0.13<energy<0.27:
        raise SystemExit("FAIL: nonphysical Taylor-Green kinetic energy")
    if uerr>0.16 or perr>0.32:
        raise SystemExit("FAIL: Taylor-Green velocity/pressure solution deviates excessively")
    return uerr,perr

if len(sys.argv)!=5:
    raise SystemExit("usage: check_tgv.py CASE12 LOG12 CASE24 LOG24")
a,pa=check(sys.argv[1],sys.argv[2],12)
b,pb=check(sys.argv[3],sys.argv[4],24)
if b>=.9*a:
    raise SystemExit("FAIL: velocity error does not improve on refined grid")
if pb>=1.1*pa:
    raise SystemExit("FAIL: pressure error worsens under refinement")
print("PASS: real pressure-corrected TGV vortex, conservation and grid refinement")
