#!/usr/bin/env python3
"""M1C-14 independent structured-cell-centre 3D fixed-position reconstruction."""
import json
import math
import re
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_gaussian_green import (
    POS, green, GL24, GL32, tail_mass, field, DIM, INV, HTF, ZTF, EM,
    NUM, START, END, RATES, WIND, SOURCE_CENTER, SIGMA, nearest
)

def fail(s):
    raise SystemExit("FAIL: M1C-14 "+s)

def grid(case,n):
    centres=field(case/"0/C",n**3,"vector","0 1 0 0 0 0 0")
    volumes=field(case/"0/V",n**3,"scalar","0 3 0 0 0 0 0")
    if min(volumes)<=0 or abs(sum(volumes)-1)>1e-9:
        fail("native warped cell-volume closure")
    columns={}
    for a,(x,y,z) in enumerate(centres):
        j=round(n*y-.5); k=round(n*z-.5)
        if not (0<=j<n and 0<=k<n) or max(abs(y-(j+.5)/n),abs(z-(k+.5)/n))>0.01:
            fail("native physical cell-centre y/z lattice geometry corrupted")
        if not 0<x<1:
            fail("native cell-centre x out of cube")
        columns.setdefault((j,k),[]).append((x,a))
    if len(columns)!=n*n:
        fail("native structured y/z columns incomplete")
    lookup={}
    for (j,k),col in columns.items():
        if len(col)!=n:
            fail("native structured x-column size corrupted")
        col.sort()
        if any(col[i+1][0]-col[i][0]<1e-9 for i in range(n-1)):
            fail("native physical x-column order/uniqueness")
        for i,(x,index) in enumerate(col):
            nominal=(i+.5)/n
            if abs(x-nominal)>.041:
                fail("native warped cell-centre x displacement unexpectedly large")
            lookup[(i,j,k)]=index
    if len(lookup)!=n**3:
        fail("native three-dimensional cell index map incomplete")
    return centres,volumes,lookup

def invert_warp_x(x,y,z):
    amp=.04*math.sin(math.pi*y)*math.sin(math.pi*z)
    lo,hi=0.,1.
    for _ in range(60):
        mid=(lo+hi)*.5
        if mid+amp*math.sin(math.pi*mid)<x:
            lo=mid
        else:
            hi=mid
    return (lo+hi)*.5

def weights(target,centres,lookup,n):
    # Structured logical indexing is for finding the native stencil only.
    # Physical weights come from the actual warped 8-centre trilinear map.
    x,y,z=target
    qx=invert_warp_x(x,y,z)
    i=max(0,min(n-2,math.floor(n*qx-.5)))
    j=max(0,min(n-2,math.floor(n*y-.5)))
    k=max(0,min(n-2,math.floor(n*z-.5)))
    nodes={(a,b,c):lookup[(i+a,j+b,k+c)]
           for a in (0,1) for b in (0,1) for c in (0,1)}
    if len(set(nodes.values()))!=8:
        fail("native 8-centre trilinear cell topology invalid")
    coords={ijk:centres[index] for ijk,index in nodes.items()}

    def mapping(t):
        u,v,w=t
        p=[0.,0.,0.]
        jac=[[0.,0.,0.] for _ in range(3)]
        for (a,b,c),xyz in coords.items():
            ux=u if a else 1-u
            vy=v if b else 1-v
            wz=w if c else 1-w
            product=ux*vy*wz
            derivatives=((1 if a else -1)*vy*wz,
                         ux*(1 if b else -1)*wz,
                         ux*vy*(1 if c else -1))
            for axis in range(3):
                p[axis]+=product*xyz[axis]
                for derivative in range(3):
                    jac[axis][derivative]+=xyz[axis]*derivatives[derivative]
        return p,jac

    def cross(a,b):
        return (a[1]*b[2]-a[2]*b[1],
                a[2]*b[0]-a[0]*b[2],
                a[0]*b[1]-a[1]*b[0])
    def dot(a,b):
        return sum(x*y for x,y in zip(a,b))
    t=[n*qx-i-.5,n*y-j-.5,n*z-k-.5]
    for _ in range(18):
        physical,J=mapping(t)
        delta=[q-p for q,p in zip(target,physical)]
        if max(map(abs,delta))<=2e-13:
            break
        col=[tuple(J[axis][idx] for axis in range(3)) for idx in range(3)]
        det=dot(col[0],cross(col[1],col[2]))
        if not math.isfinite(det) or det<=1e-12:
            fail("degenerate or inverted physical native trilinear interpolation map")
        correction=(dot(delta,cross(col[1],col[2]))/det,
                    dot(col[0],cross(delta,col[2]))/det,
                    dot(col[0],cross(col[1],delta))/det)
        t=[q+d for q,d in zip(t,correction)]
    physical,J=mapping(t)
    if max(abs(p-q) for p,q in zip(physical,target))>1e-9:
        fail("nonorthogonal native trilinear physical inverse residual")
    if not all(-1e-12<=q<=1+1e-12 for q in t):
        fail("physical receptor outside native 8-centre trilinear stencil")
    u,v,w=t
    anchors=[]
    for (a,b,c),idx in nodes.items():
        weight=(u if a else 1-u)*(v if b else 1-v)*(w if c else 1-w)
        anchors.append((idx,weight))
    if len({a for a,_ in anchors})!=8 or min(weight for _,weight in anchors)<-1e-12:
        fail("native eight-centre weights invalid")
    if abs(sum(weight for _,weight in anchors)-1)>1e-12:
        fail("native trilinear weights not conservative")
    def interpolate(data):
        return sum(weight*data[index] for index,weight in anchors)
    affine=((1.,[1.]*len(centres)),
            (x,[c[0] for c in centres]),
            (y,[c[1] for c in centres]),
            (z,[c[2] for c in centres]),
            (5+2*x-3*y+4*z,[5+2*a-3*b+4*c for a,b,c in centres]))
    if any(abs(interpolate(values)-answer)>1e-9 for answer,values in affine):
        fail("native nonorthogonal 3-D physical affine interpolation failed")
    return anchors,interpolate

def audit(case,n):
    if n not in (8,12,16,24):
        fail("unsupported frozen mesh")
    centres,volumes,lookup=grid(case,n)
    props=(case/"constant/transportProperties").read_text()
    for token in (
        "sourceShape gaussian;",
        "sourceCenter [0 1 0 0 0 0 0] (0.42 0.55 0.5);",
        "sourceSigma [0 1 0 0 0 0 0] 0.08;",
        "sourceHT_Bq_s [0 0 -2 0 0 0 0] 60000;",
        "sourceHTO_Bq_s [0 0 -2 0 0 0 0] 12000;",
        "sourceStart [0 0 1 0 0 0 0] 0.073;",
        "sourceEnd [0 0 1 0 0 0 0] 0.137;",
        "kOx [0 0 -1 0 0 0 0] 0.4;",
        "lambdaT [0 0 -1 0 0 0 0] 0.2;",
        "DHT [0 2 -1 0 0 0 0] 0.01;",
        "DHTO [0 2 -1 0 0 0 0] 0.02;",
        "momentumMode prescribed;",
        "activityMode htHto;",
    ):
        if token not in props:
            fail("native physical source/wind/reactive parameter modified "+token)
    ctl=(case/"system/controlDict").read_text()
    if not all(v in ctl for v in ("endTime 0.4;","deltaT 0.005;","writeInterval 80;")):
        fail("native 80-step time contract changed")
    for part in ("0/HT","0/HTO","0/U","0.4/HT","0.4/HTO","0.4/U"):
        t=(case/part).read_text()
        for patch in ("xMin","xMax","yMin","yMax","zMin","zMax"):
            if not re.search(rf"\b{patch}\s*\{{\s*type\s+cyclic;",t):
                fail("physical periodic species patch corrupted")
    data={species:field(case/"0.4"/species,n**3,"scalar",DIM) for species in ("HT","HTO")}
    ini={species:field(case/"0"/species,n**3,"scalar",DIM) for species in ("HT","HTO")}
    if any(min(v)<-1e-12 or not all(math.isfinite(x) for x in v) for v in data.values()):
        fail("native real Gaussian activity contains nonfinite/negative field values")
    if max(max(map(abs,ini[s])) for s in ini)>1e-8:
        fail("native initial state not zero")
    for p in ("0/U","0.4/U"):
        vv=field(case/p,n**3,"vector","0 1 -1 0 0 0 0")
        if max(abs(a-b) for v in vv for a,b in zip(v,WIND))>1e-12:
            fail("native physical prescribed vector U changed")
    s=(case/"solve.log").read_text()
    times=[float(q) for q in re.findall(r"(?m)^Time = ("+NUM+r")\s*$",s)]
    logs=[tuple(map(float,m)) for m in INV.findall(s)]
    ff=[tuple(map(float,m)) for m in HTF.findall(s)]
    zz=[tuple(map(float,m)) for m in ZTF.findall(s)]
    em=[tuple(map(float,m)) for m in EM.findall(s)]
    if any(len(v)!=80 for v in (times,logs,ff,zz,em)):
        fail("native every-step Gaussian activity evidence absent")
    if any(abs(q-.005*i)>2e-12 for i,q in enumerate(times,1)):
        fail("native real simulation clock incorrect")
    if any(len(re.findall(r"Solving for "+sp+r",",s))!=80 for sp in ("HT","HTO","c")):
        fail("native C++ HT/HTO/passive 80-step solve incomplete")
    prev=[0.,0.]
    totals=[0.,0.]
    max_budget=0.
    for i,((h,z),(eh,ez,sec),(lh,rh,fh,bh),(lz,rz,fz,bz,bt)) in enumerate(zip(logs,em,ff,zz),1):
        overlap=max(0.,min(i*.005,END)-max((i-1)*.005,START))
        source_h,source_z=60000*overlap,12000*overlap
        totals[0]+=source_h; totals[1]+=source_z
        if max(abs(eh-source_h),abs(ez-source_z),abs(sec-overlap))>1e-8:
            fail("genuine Gaussian finite-window emission mismatch")
        if max(abs(fh),abs(fz),abs(lh+rh),abs(lz+rz))>1e-5:
            fail("native cyclic net boundary species flux invalid")
        rh_res=h-prev[0]+.005*fh+.005*.6*h-source_h
        rz_res=z-prev[1]+.005*fz+.005*.2*z-.005*.4*h-source_z
        tot_res=h+z-sum(prev)+.005*(fh+fz)+.005*.2*(h+z)-source_h-source_z
        max_budget=max(max_budget,*map(abs,(rh_res,rz_res,tot_res,bh,bz,bt)))
        if max_budget>1e-5:
            fail("native signed per-step physical HT/HTO Bq budget")
        prev=[h,z]
    if max(abs(totals[0]-3840),abs(totals[1]-768))>1e-8:
        fail("genuine finite source Bq total mismatch")
    for j,sp in enumerate(("HT","HTO")):
        if abs(sum(vol*x for vol,x in zip(volumes,data[sp]))-logs[-1][j])>1e-5:
            fail("native final source field integrated Bq differs from solver log "+sp)
    records={}
    for name,probe in POS.items():
        anchors,interpolator=weights(probe,centres,lookup,n)
        reference24=green(probe,GL24)
        reference32=green(probe,GL32)
        if max(abs(a-b) for a,b in zip(reference24,reference32))>0.05:
            fail("fixed physical reference quadrature not converged")
        per_species={}
        for si,sp in enumerate(("HT","HTO")):
            native=interpolator(data[sp])
            anchor_green=sum(w*green(centres[index],GL24)[si] for index,w in anchors)
            exact=reference32[si]
            physical_difference=native-anchor_green
            observation_difference=anchor_green-exact
            total=native-exact
            if abs(total-physical_difference-observation_difference)>1e-7:
                fail("fixed physical 3D interpolation error decomposition not closed")
            if not math.isfinite(native) or native<=0 or exact<=10:
                fail("fixed physical probe nonfinite, nonpositive or below reliable reference")
            relative=abs(total)/exact
            if relative>1.5:
                fail("fixed physical probe pathological error exceeds frozen bound")
            per_species[sp]={"interpolated_native_Bq_m3":native,
                            "fixed_analytic_Bq_m3":exact,
                            "signed_error_Bq_m3":total,
                            "relative_error":relative,
                            "interpolated_FVM_minus_analytic_anchors_Bq_m3":physical_difference,
                            "interpolation_of_analytic_minus_exact_fixed_Bq_m3":observation_difference}
            print(f"M1C14_FIXED N={n} {name} {sp} native={native:.9g} green={exact:.9g} rel={relative:.8g} FVM={physical_difference:.9g} observation={observation_difference:.9g}")
        nearest_i,drift=nearest(centres,probe)
        records[name]={"point_m":probe,"native_nearest_index":nearest_i,"native_nearest_distance_m":drift,
                       "anchor_indices":[a for a,w in anchors],
                       "anchor_weights":[w for a,w in anchors],
                       "species":per_species}
    print(f"M1C14_CASE N={n} cells={n**3} Bqmax={max_budget:.9g} sourceHT={totals[0]:.8g} sourceHTO={totals[1]:.8g}")
    return {"N":n,"cells":n**3,"receptors":records,"max_budget_Bq":max_budget,"emitted_Bq":totals}

def matrix(root):
    cases={}
    for n in (8,12,16,24):
        cases[n]=audit(root/f"{n}_normal",n)
    for name in POS:
        for sp in ("HT","HTO"):
            errs=[cases[n]["receptors"][name]["species"][sp]["relative_error"] for n in (8,12,16,24)]
            print(f"M1C14_TREND {name} {sp} rel_errors_8_12_16_24={errs} finest_improves_vs_8={errs[-1]<errs[0]} strictly_monotone={all(a>b for a,b in zip(errs,errs[1:]))}")
    report={"scope":"3-D independent fixed-coordinate cell-centre reconstruction; not OpenFOAM interpolationCellPoint; native warped periodic synthetic Gaussian",
            "parent_sha":"39b6393113f4f1f63af449d927a2fc2d0f99f6a6",
            "cases":cases,"tail_outside_cube":tail_mass(),
            "conclusion":"SCOPED_FIXED_POSITION_RECEPTOR_RECONSTRUCTION_DIAGNOSTIC"}
    (root/"m1c14-report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("M1C14_REPORT_JSON="+json.dumps(report,separators=(",",":")))
    print("PASS: M1C-14 four native source-reaction cases, real fixed-location interpolation/affine validation, independent Green error decomposition and strict 3D Bq budgets")

if __name__=="__main__":
    if len(sys.argv)==4 and sys.argv[1]=="--single":
        audit(Path(sys.argv[2]),int(sys.argv[3]))
        print("PASS: M1C-14 native fixed physical point audit")
    elif len(sys.argv)==3 and sys.argv[1]=="--matrix":
        matrix(Path(sys.argv[2]))
    else:raise SystemExit("usage: check_fixed_receptors.py --single CASE N | --matrix ROOT")
