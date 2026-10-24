#!/usr/bin/env python3
"""Independent actual native variable-D source, Bq, cells and uniform-limit check."""
import math,re,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_robustness import records
from check_gaussian_green import field,DIM,NUM,INV
from check_fixed_receptors import grid
from check_convection_scheme import EXTREMA

def fail(msg):raise SystemExit("FAIL: M1C22 "+msg)

def check(case,n,mode):
    s=(case/"solve.log").read_text()
    inv,budget=records(case,80,.005)
    extrema=[tuple(map(float,x)) for x in EXTREMA.findall(s)]
    if len(extrema)!=80 or any(not all(map(math.isfinite,row)) or row[0]<-1e-8 or row[2]<-1e-8 for row in extrema):
        fail("invalid native real-species per-step min/max")
    centres,volumes,lookup=grid(case,n)
    for i,name in enumerate(("HT","HTO")):
        a=field(case/"0.4"/name,n**3,'scalar',DIM)
        if not all(map(math.isfinite,a)) or min(a)<-1e-8:fail("nonfinite/negative genuine primitive cells "+name)
        if abs(sum(x*v for x,v in zip(a,volumes))-inv[-1][i])>1e-5:
            fail("native logged Bq vs real integrated primitive cells "+name)
    pat=r"CTEAS_M1C22_SPATIAL_DIFFUSIVITY species=(DHTSpatial|DHTOSpatial) min=("+NUM+r") max=("+NUM+r")"
    values=re.findall(pat,s)
    if mode=="base":
        if values:fail("opt-out generated variable-coefficient logs")
    else:
        if len(values)!=2:fail("actual physical spatial diffusivity input logs missing")
        for name,low,high in values:
            a,b=float(low),float(high)
            if not all(map(math.isfinite,(a,b))) or a<=0 or b<a:
                fail("invalid native spatial diffusivity global extrema")
            if mode=="variable" and (b/a<1.2):fail("claimed spatial diffusivity insufficiently varying")
            if mode=="uniform" and (abs(b/a-1)>1e-12):fail("claimed uniform field actually nonuniform")
    print(f"M1C22_NATIVE mode={mode} N={n} max_signed_Bq={budget:.16g}")
    return inv

def compare(root,uniform,variable,n):
    b=check(root,n,"base");u=check(uniform,n,"uniform");check(variable,n,"variable")
    for step,(v,w) in enumerate(zip(b,u),1):
        if max(abs(a-c) for a,c in zip(v,w))>1e-5:
            fail(f"constant input spatial field native Bq not equal scalar baseline step={step}")
    for name in ("HT","HTO"):
        a=field(root/"0.4"/name,n**3,"scalar",DIM)
        b=field(uniform/"0.4"/name,n**3,"scalar",DIM)
        c=field(variable/"0.4"/name,n**3,"scalar",DIM)
        equal=max(abs(x-y) for x,y in zip(a,b))
        changed=max(abs(x-y) for x,y in zip(a,c))
        if equal>3e-4:fail(f"uniform field not equal original scalar branch {name}")
        if changed<0.01:fail(f"nonconstant diffusivity had no true numerical effect {name}")
        print(f"M1C22_EQUIVALENCE N={n} {name} constant_field_delta_Bq_m3={equal:.12g} variable_D_effect_Bq_m3={changed:.12g}")
    print(f"PASS: M1C22 actual N={n} 3D native variable-D conservative HT/HTO transport and uniform-field scalar equivalence")
if __name__=="__main__":
    a=sys.argv[1:]
    if len(a)==3:check(Path(a[0]),int(a[1]),a[2])
    elif len(a)==5 and a[0]=="--compare":compare(Path(a[1]),Path(a[2]),Path(a[3]),int(a[4]))
    else:raise SystemExit("usage CASE N MODE | --compare BASE UNIFORM VARIABLE N")
