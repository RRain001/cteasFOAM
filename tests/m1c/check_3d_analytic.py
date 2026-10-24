#!/usr/bin/env python3
"""Independent continuum 3-D analytic eigenmode / native signed Bq verification."""
import json
import math
import re
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_gaussian_source import field

NUM=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VEC=re.compile(r"\(\s*("+NUM+r")\s+("+NUM+r")\s+("+NUM+r")\s*\)")
INV=re.compile(r"CTEAS_M1C activity_HT_Bq=("+NUM+r") activity_HTO_Bq=("+NUM+r")")
HTF=re.compile(r"CTEAS_M1C flux_HT_xMin_Bq_s=("+NUM+r") flux_HT_xMax_Bq_s=("+NUM+r") flux_HT_net_Bq_s=("+NUM+r") flux_HT_budget_residual_Bq=("+NUM+r")")
ZTF=re.compile(r"CTEAS_M1C flux_HTO_xMin_Bq_s=("+NUM+r") flux_HTO_xMax_Bq_s=("+NUM+r") flux_HTO_net_Bq_s=("+NUM+r") flux_HTO_budget_residual_Bq=("+NUM+r") flux_total_budget_residual_Bq=("+NUM+r")")
U={"advect":(.2,.1,.07),"diffuse":(0.,0.,0.)}
DIM="0 -3 -1 0 0 0 0"
PI2=2*math.pi
DT=.005
T=.2
THRESHOLDS={8:(80000,30000),12:(60000,25000),16:(50000,20000)}

def fail(msg):raise SystemExit("FAIL: M1C-11 "+msg)

def physical_geometry(case,n,nt):
    meshtext=(case/"mesh-check.log").read_text()
    if "Mesh OK." not in meshtext or not re.search(rf"(?m)^\s*cells:\s*{n}\s*$",meshtext):
        fail("real native checkMesh cell count/topology")
    if not (case/"constant/polyMesh/owner").is_file():
        fail("native OpenFOAM owner file absent")
    boundary=(case/"constant/polyMesh/boundary").read_text()
    for patch in ("xMin","xMax","yMin","yMax","zMin","zMax"):
        if not re.search(rf"\b{patch}\s*\{{\s*type\s+cyclic;",boundary):
            fail("3-D native cyclic boundary absent "+patch)
    text=(case/"constant/polyMesh/points").read_text()
    m=re.search(r"(?m)^\s*(\d+)\s*\n\s*\(\s*\n",text)
    if not m:fail("native warped points absent")
    end=re.search(r"(?m)^\s*\)\s*$",text[m.end():])
    if not end:fail("native warped points malformed")
    rows=[tuple(map(float,row)) for row in VEC.findall(text[m.end():m.end()+end.start()])]
    if len(rows)!=(nt+1)**3 or len(rows)!=int(m.group(1)):
        fail("native vertex count")
    moved=0
    for x,y,z in rows:
        iy,iz=round(y*nt),round(z*nt)
        if not (0<=iy<=nt and 0<=iz<=nt):fail("native y/z lattice bounds")
        by,bz=iy/nt,iz/nt
        if max(abs(by-y),abs(bz-z))>2e-10:fail("native y/z vertex mismatch")
        amp=.04*math.sin(math.pi*y)*math.sin(math.pi*z)
        lo,hi=max(0.,x-.04),min(1.,x)
        for _ in range(56):
            mid=(lo+hi)/2
            if mid+amp*math.sin(math.pi*mid)<x:lo=mid
            else:hi=mid
        bx=round((lo+hi)*.5*nt)/nt
        boundary_vertex=bx in (0,1) or by in (0,1) or bz in (0,1)
        delta=0. if boundary_vertex else .04*math.sin(math.pi*bx)*math.sin(math.pi*by)*math.sin(math.pi*bz)
        if abs(x-bx-delta)>2e-10:fail("native interior warped physical vertex mismatch")
        moved+=abs(delta)>1e-12
    if moved!=(nt-1)**3:fail("3-D interior warp missing")
    c=field(case/"0/C",n,"vector","0 1 0 0 0 0 0")
    v=field(case/"0/V",n,"scalar","0 3 0 0 0 0 0")
    if min(v)<=0 or abs(sum(v)-1)>1e-9:fail("native mesh volumes invalid")
    if any(any(not (0<q<1) for q in row) for row in c):
        fail("native cell centre outside cube")
    return c,v,moved

def expected(c,mode,sp,time):
    u=U[mode]
    k2=3*PI2**2
    a=.01*k2+.6
    b=.02*k2+.2
    amp_shape=math.prod(math.cos(PI2*(pos-time*ui)) for pos,ui in zip(c,u))
    if sp=="HT":
        return 1e6*math.exp(-.6*time)+2e5*math.exp(-a*time)*amp_shape
    return (math.exp(-.2*time)*(2e5+1e6*(-math.expm1(-.4*time)))
        +(5e4*math.exp(-b*time)+.4*2e5*(math.exp(-a*time)-math.exp(-b*time))/(b-a))*amp_shape)

def audit(case,mode,nt):
    if mode not in U or nt not in THRESHOLDS:fail("invalid frozen mode / resolution")
    n=nt**3
    c,v,moved=physical_geometry(case,n,nt)
    props=(case/"constant/transportProperties").read_text()
    for required in ("momentumMode prescribed;","activityMode htHto;",
                     "kOx [0 0 -1 0 0 0 0] 0.4;","lambdaT [0 0 -1 0 0 0 0] 0.2;",
                     "DHT [0 2 -1 0 0 0 0] 0.01;","DHTO [0 2 -1 0 0 0 0] 0.02;",
                     "auditActivityFlux true;"):
        if required not in props:fail("configured physics mismatch: "+required)
    if any(x in props for x in ("sourceHT_Bq_s","sourceHTO_Bq_s","sourceStart","sourceEnd","sourceShape")):
        fail("source-free experiment configured an emission")
    ctl=(case/"system/controlDict").read_text()
    for required in ("endTime 0.2;","deltaT 0.005;","writeInterval 40;"):
        if required not in ctl:fail("native clock/steps mismatch")
    schemes=(case/"system/fvSchemes").read_text()
    for required in ("div(phi,HT) Gauss upwind;","div(phi,HTO) Gauss upwind;",
                     "laplacianSchemes { default Gauss linear corrected; }"):
        if required not in schemes:fail("FVM scheme mismatch")
    for f in ("0/U","0/c","0/HT","0/HTO","0.2/U","0.2/HT","0.2/HTO"):
        ftxt=(case/f).read_text()
        for patch in ("xMin","xMax","yMin","yMax","zMin","zMax"):
            if not re.search(rf"\b{patch}\s*\{{\s*type\s+cyclic;",ftxt):
                fail("native field not three-dimensional periodic "+f+" "+patch)
    u0=field(case/"0/U",n,"vector","0 1 -1 0 0 0 0")
    uf=field(case/"0.2/U",n,"vector","0 1 -1 0 0 0 0")
    if max(abs(val-true) for actual in (u0,uf) for row in actual for val,true in zip(row,U[mode]))>1e-12:
        fail("actual native prescribed U does not match frozen mode")
    s=(case/"solve.log").read_text()
    times=[float(t) for t in re.findall(r"(?m)^Time = ("+NUM+r")\s*$",s)]
    h=[tuple(map(float,t)) for t in INV.findall(s)]
    hf=[tuple(map(float,t)) for t in HTF.findall(s)]
    zf=[tuple(map(float,t)) for t in ZTF.findall(s)]
    if any(len(t)!=40 for t in (times,h,hf,zf)) or any(abs(t-DT*i)>2e-12 for i,t in enumerate(times,1)):
        fail("40 independent native time and activity steps required")
    if any(len(re.findall(r"Solving for "+name+r",",s))!=40 for name in ("HT","HTO","c")):
        fail("40 actual HT/HTO/c native matrix solves required")
    if "Solving for p," in s or "Solving for U" in s or "sourceShape=gaussian" in s or "emitted_HT_Bq" in s:
        fail("prescribed source-free simulation performed other physics")
    if "mode=prescribed" not in s or "Time = 0.2" not in s:
        fail("prescribed final native mode not solved")
    ini={sp:field(case/"0"/sp,n,"scalar",DIM) for sp in ("HT","HTO")}
    fin={sp:field(case/"0.2"/sp,n,"scalar",DIM) for sp in ("HT","HTO")}
    err={};maxerr={}
    for sp in ("HT","HTO"):
        if min(ini[sp]+fin[sp])<=0:fail("native species nonpositive")
        ierror=max(abs(actual-expected(pos,mode,sp,0)) for actual,pos in zip(ini[sp],c))
        if ierror>1e-6:fail("invalid native centre 3-D initial wave "+sp)
        residual=[actual-expected(pos,mode,sp,T) for actual,pos in zip(fin[sp],c)]
        err[sp]=math.sqrt(sum(vol*dr**2 for vol,dr in zip(v,residual))/sum(v))
        maxerr[sp]=max(map(abs,residual))
        lim=THRESHOLDS[nt][0 if sp=="HT" else 1]
        if not 0<err[sp]<lim or maxerr[sp]>lim:
            fail(f"{mode} {nt} {sp} analytic RMS={err[sp]:.9g} Linf={maxerr[sp]:.9g} exceeds {lim}")
        if max(fin[sp])-min(fin[sp])<1e3:
            fail("three-dimensional species eigenmode flattened")
        integral=sum(vol*val for vol,val in zip(v,fin[sp]))
        logged=h[-1][0 if sp=="HT" else 1]
        if abs(integral-logged)>1e-5:
            fail("actual native field Bq integral differs from solver log "+sp)
    old=[sum(vi*xi for vi,xi in zip(v,ini[sp])) for sp in ("HT","HTO")]
    maxbalance=0.
    minflux=float("inf")
    for i,((newh,newz),fht,fz) in enumerate(zip(h,hf,zf),1):
        xl,xr,net,loggedht=fht
        yl,yr,netz,loggedz,loggedtotal=fz
        vals=[newh,newz,*fht,*fz]
        if not all(math.isfinite(q) for q in vals) or min(newh,newz)<=0:
            fail("nonfinite/negative actual time-integrated inventory/flux")
        if abs(net)>1e-5 or abs(netz)>1e-5 or abs(xl+xr)>1e-5 or abs(yl+yr)>1e-5:
            fail("three dimensional cyclic net face flux not balanced")
        minflux=min(minflux,abs(xl),abs(xr),abs(yl),abs(yr))
        ht=newh-old[0]+DT*net+DT*.6*newh
        z=newz-old[1]+DT*netz+DT*.2*newz-DT*.4*newh
        total=newh+newz-sum(old)+DT*(net+netz)+DT*.2*(newh+newz)
        maxbalance=max(maxbalance,*map(abs,(ht,z,total,loggedht,loggedz,loggedtotal,
                                             ht-loggedht,z-loggedz,total-loggedtotal)))
        if maxbalance>1e-5:fail(f"step {i} native HT/HTO signed Bq budget failure")
        old=[newh,newz]
    print(f"M1C11_CASE {mode} N={nt} cells={n} realwarped={moved} "
          f"HT_RMS={err['HT']:.10g} HTO_RMS={err['HTO']:.10g} "
          f"HT_Linf={maxerr['HT']:.10g} HTO_Linf={maxerr['HTO']:.10g} "
          f"maxBudget_Bq={maxbalance:.9g} xCyclicMin_Bq_s={minflux:.9g}")
    return {"mode":mode,"N":nt,"cells":n,"HT_RMS":err["HT"],"HTO_RMS":err["HTO"],
            "HT_Linf":maxerr["HT"],"HTO_Linf":maxerr["HTO"],"max_budget_Bq":maxbalance,
            "x_min_abs_cyclic_flux_Bq_s":minflux,"midpoint_nonuniform":True}

def matrix(root):
    out={mode:{nt:audit(root/f"{mode}_{nt}",mode,nt) for nt in (8,12,16)}
         for mode in U}
    for mode in U:
        for sp in ("HT","HTO"):
            ratios=[out[mode][b][sp+"_RMS"]/out[mode][a][sp+"_RMS"] for a,b in ((8,12),(12,16))]
            print(f"M1C11_REFINEMENT {mode} {sp} RMS_8={out[mode][8][sp+'_RMS']:.10g} "
                  f"RMS_12={out[mode][12][sp+'_RMS']:.10g} "
                  f"RMS_16={out[mode][16][sp+'_RMS']:.10g} ratios={ratios}")
            if not all(0<r<.94 for r in ratios):
                fail(f"frozen two-step 3-D continuum RMS refinement for {mode} {sp}: {ratios}")
    # Comparing real native final U-advected and U-zero fields is an independent advection smoke.
    for nt in (8,12,16):
        n=nt**3
        hmove=field(root/f"advect_{nt}/0.2/HT",n,"scalar",DIM)
        hstill=field(root/f"diffuse_{nt}/0.2/HT",n,"scalar",DIM)
        difference=max(abs(a-b) for a,b in zip(hmove,hstill))
        if difference<3000:fail("native 3-D nonzero U had no observable advection")
        print(f"M1C11_NATIVE_ADVECT_MINUS_DIFFUSE N={nt} HT_maxdiff_Bq_m3={difference:.9g}")
    report={"scope":"source-free genuinely 3-D advective-diffusive-reactive analytic eigenmode only",
            "native_cases":out,"source_terms_configured":False,
            "conclusion":"SOURCE_FREE_3D_ANALYTIC_VERIFICATION_SCOPED"}
    (root/"m1c11-report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("M1C11_REPORT_JSON="+json.dumps(report,sort_keys=True,separators=(",",":")))
    print("PASS: M1C-11 six real v2606 source-free warped 3-D periodic native analytical fields, Bq budgets, two refinement ratios and genuine transport")

if __name__=="__main__":
    if len(sys.argv)==5 and sys.argv[1]=="--single":
        audit(Path(sys.argv[2]),sys.argv[3],int(sys.argv[4]))
        print("PASS: M1C-11 source-free native analytic and signed Bq case audit")
    elif len(sys.argv)==3 and sys.argv[1]=="--matrix":matrix(Path(sys.argv[2]))
    else:raise SystemExit("usage: check_3d_analytic.py --single CASE MODE N | --matrix ROOT")
