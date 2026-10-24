#!/usr/bin/env python3
"""Independent native-cell finite-window Gaussian source Bq budget checker."""
from pathlib import Path
import math
import re
import sys

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
INV = re.compile(r"CTEAS_M1C activity_HT_Bq=(" + NUM + r") activity_HTO_Bq=(" + NUM + r")")
EM = re.compile(r"CTEAS_M1C emitted_HT_Bq=(" + NUM + r") emitted_HTO_Bq=(" + NUM + r") source_overlap_s=(" + NUM + r")")
BUD = re.compile(r"CTEAS_M1C closed_HT_residual_Bq=(" + NUM + r") closed_HTO_residual_Bq=(" + NUM + r") closed_total_residual_Bq=(" + NUM + r")")

def fail(msg):
    raise SystemExit("FAIL: M1C-6 " + msg)

def field(path, count, kind, dimensions):
    if not path.is_file():
        fail(f"missing actual native OpenFOAM field {path}")
    s = path.read_text()
    if not re.search(r"dimensions\s+\[" + re.escape(dimensions) + r"\]", s):
        fail(f"wrong native units {path}")
    if kind == "scalar":
        mat = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", s, re.S)
        if mat:
            if int(mat.group(1)) != count: fail("native scalar length mismatches mesh")
            vals = [float(x) for x in mat.group(2).split()]
        else:
            mat = re.search(r"internalField\s+uniform\s+("+NUM+r")\s*;", s)
            vals = [float(mat.group(1))]*count if mat else []
    else:
        mat = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", s, re.S)
        if mat:
            if int(mat.group(1)) != count: fail("native vector length mismatches mesh")
            vals = [tuple(map(float,v)) for v in re.findall(
                r"\(\s*("+NUM+r")\s+("+NUM+r")\s+("+NUM+r")\s*\)",mat.group(2))]
        else:
            mat = re.search(r"internalField\s+uniform\s+\(("+NUM+r")\s+("+NUM+r")\s+("+NUM+r")\)\s*;",s)
            vals = [tuple(map(float,mat.groups()))]*count if mat else []
    if len(vals)!=count or not all(
        math.isfinite(v) for row in vals for v in (row if kind=="vector" else (row,))
    ):
        fail(f"missing/invalid actual native field {path}")
    return vals

def history(log):
    s=Path(log).read_text()
    times=[float(x) for x in re.findall(r"(?m)^Time = ("+NUM+r")\s*$",s)]
    ih=[tuple(map(float,v)) for v in INV.findall(s)]
    em=[tuple(map(float,v)) for v in EM.findall(s)]
    bu=[tuple(map(float,v)) for v in BUD.findall(s)]
    if any(len(v)!=20 for v in (times,ih,em,bu)):
        fail("20 genuine time/inventory/emission/budget steps missing")
    if len(re.findall("Solving for HT,",s))!=20 or len(re.findall("Solving for HTO,",s))!=20:
        fail("native 20 HT/HTO matrix solves missing")
    if not all(math.isfinite(x) for row in ih+em+bu for x in row):
        fail("nonfinite solver activity logs")
    if "CTEAS_M1C sourceShape=gaussian" not in s:
        fail("native Gaussian spatial branch was not run")
    return times,ih,em,bu

def check(case,log,mode):
    if mode not in ("coarse","fine"): fail("unknown case")
    count = 24 if mode=="coarse" else 192
    case=Path(case)
    s=(case/"mesh-check.log").read_text()
    if ("Mesh OK." not in s or not re.search(
        rf"(?m)^\s*cells:\s*{count}\s*$",s)):
        fail("wrong genuine cell-count or native checkMesh")
    if not (case/"constant/polyMesh/owner").is_file():
        fail("missing native mesh owner")
    props=(case/"constant/transportProperties").read_text()
    for token in ("sourceShape gaussian;", "sourceCenter [0 1 0 0 0 0 0] (0.38 0.62 0.55);",
                  "sourceSigma [0 1 0 0 0 0 0] 0.22;",
                  "sourceHT_Bq_s [0 0 -2 0 0 0 0] 60000;",
                  "sourceHTO_Bq_s [0 0 -2 0 0 0 0] 12000;",
                  "sourceStart [0 0 1 0 0 0 0] 0.074999999999999997;",
                  "sourceEnd [0 0 1 0 0 0 0] 0.22500000000000001;",
                  "auditClosedActivity true;"):
        if token not in props: fail(f"source input contract mismatch {token}")
    centers=field(case/"0/C",count,"vector","0 1 0 0 0 0 0")
    vols=field(case/"0/V",count,"scalar","0 3 0 0 0 0 0")
    if min(vols)<=0 or abs(sum(vols)-1)>1e-10:
        fail("invalid positive actual native cell volumes")
    center=(.38,.62,.55); sigma=.22; dt=.05
    raw=[math.exp(-sum((x-y)**2 for x,y in zip(pos,center))/(2*sigma**2)) for pos in centers]
    W=sum(v*w for v,w in zip(vols,raw))
    if not math.isfinite(W) or W<=0 or min(raw)<=0:
        fail("empty/invalid Gaussian support")
    if max(raw)/min(raw)<1.5:
        fail("nonlocal uniform surrogate instead of Gaussian")
    centroid=[sum(v*w*p[j] for v,w,p in zip(vols,raw,centers))/W for j in range(3)]
    drift=math.sqrt(sum((a-b)**2 for a,b in zip(centroid,center)))
    if drift>.15 or not all(0<x<1 for x in centroid):
        fail(f"off-target native Gaussian moment {centroid} drift {drift}")
    zero=[0.]*count
    for name in ("HT","HTO"):
        vals=field(case/"0"/name,count,"scalar","0 -3 -1 0 0 0 0")
        if max(map(abs,vals))>1e-12: fail("source benchmark not initially empty")
    h=list(zero);z=list(zero)
    trows,ih,em,bu=history(log)
    maxfield=maxbudget=maxtotal=maxinv=maxem=0.
    released=[0.,0.]
    for step,(t,(reported_h,reported_z),(emh,emz,overlap),(bh,bz,bt)) in enumerate(zip(trows,ih,em,bu),1):
        if abs(t-step*dt)>2e-12: fail("native time grid changed")
        tau=max(0.,min(t,.225)-max(t-dt,.075))
        jh,jz=60000*tau,12000*tau
        maxem=max(maxem,abs(overlap-tau),abs(emh-jh),abs(emz-jz))
        released[0]+=jh; released[1]+=jz
        oldh=sum(x*v for x,v in zip(h,vols));oldz=sum(x*v for x,v in zip(z,vols))
        h=[(x+jh*w/W)/(1+dt*.6) for x,w in zip(h,raw)]
        z=[(x+jz*w/W+dt*.4*y)/(1+dt*.2) for x,w,y in zip(z,raw,h)]
        newh=sum(x*v for x,v in zip(h,vols));newz=sum(x*v for x,v in zip(z,vols))
        maxinv=max(maxinv,abs(newh-reported_h),abs(newz-reported_z))
        ch=reported_h-(ih[step-2][0] if step>1 else 0)+dt*.6*reported_h-jh
        cz=reported_z-(ih[step-2][1] if step>1 else 0)+dt*.2*reported_z-dt*.4*reported_h-jz
        ct=(reported_h+reported_z)-(sum(ih[step-2]) if step>1 else 0)+dt*.2*(reported_h+reported_z)-jh-jz
        maxbudget=max(maxbudget,*map(abs,(ch,cz,bh,bz,ch-bh,cz-bz)))
        maxtotal=max(maxtotal,*map(abs,(ct,bt,ct-bt)))
        if min(reported_h,reported_z)<-1e-12:fail("negative reported inventories")
    final_h=field(case/"1/HT",count,"scalar","0 -3 -1 0 0 0 0")
    final_z=field(case/"1/HTO",count,"scalar","0 -3 -1 0 0 0 0")
    if min(final_h+final_z)<-1e-12:fail("negative native final field")
    maxfield=max(max(abs(x-y) for x,y in zip(h,final_h)),
                 max(abs(x-y) for x,y in zip(z,final_z)))
    maxinv=max(maxinv,abs(sum(x*v for x,v in zip(final_h,vols))-ih[-1][0]),
               abs(sum(x*v for x,v in zip(final_z,vols))-ih[-1][1]))
    if maxfield>1e-6 or maxinv>1e-6 or maxbudget>1e-6 or maxtotal>2e-6 or maxem>1e-9:
        fail(f"Gauss BE source mismatch: field {maxfield}, inventory {maxinv}, budget {maxbudget}, total {maxtotal}, injection {maxem}")
    if abs(released[0]-9000)>1e-9 or abs(released[1]-1800)>1e-9:
        fail("total integrated finite-window release is not exact")
    # Independent integration over all native cell volumes.
    if any(abs(sum(j*w*v/W for w,v in zip(raw,vols))-j)>1e-9 for j in released):
        fail("spatial normalization changed total Bq")
    print(f"PASS: M1C-6 mode={mode} cells={count} center_drift_m={drift:.12g} gaussian_ratio={max(raw)/min(raw):.12g} "
          f"max_field_Bq_m3={maxfield:.12g} max_inventory_Bq={maxinv:.12g} "
          f"max_budget_Bq={maxbudget:.12g} max_total_Bq={maxtotal:.12g} emittedHT_Bq={released[0]:.12g} emittedHTO_Bq={released[1]:.12g}")
    return ih,final_h,final_z

if __name__=="__main__":
    if len(sys.argv)!=4:raise SystemExit("usage: check_gaussian_source.py CASE LOG coarse|fine")
    check(*sys.argv[1:])
