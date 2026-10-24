#!/usr/bin/env python3
"""Independent fail-closed 5-case native OpenFOAM M1C-8 refinement audit."""
from pathlib import Path
import json
import math
import re
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_gaussian_source import field
from check_open_gaussian_release import INV,EM,HTF,ZTF,FLOW,DIV,records

CASES={
 "space_lo": (6,16,6,.02,20),
 "baseline": (8,20,6,.02,20),
 "space_hi": (12,30,8,.02,20),
 "time_lo": (8,20,6,.04,10),
 "time_hi": (8,20,6,.01,40)
}
RECEPTORS={"R1":(.60,.62,.55),"R2":(.75,.50,.50),"R3":(.88,.62,.55)}
VNUM=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VEC=re.compile(r"\(\s*("+VNUM+r")\s+("+VNUM+r")\s+("+VNUM+r")\s*\)")
# Dynamic-step analogue of previously frozen M1C-7 readback, without importing
# its 960-cell/20-step assumptions.
def fail(s):raise SystemExit("FAIL: M1C-8 "+s)
def rows(s,rg,name,n):
    values=[tuple(map(float,r)) for r in rg.findall(s)]
    if len(values)!=n or any(not math.isfinite(x) for row in values for x in row):
        fail(f"{name}: invalid finite {n}-step native records")
    return values
def get(case,mode):
    case=Path(case)
    nx,ny,nz,dt,steps=CASES[mode]; n=nx*ny*nz
    if not (case/"constant/polyMesh/owner").is_file():
        fail("missing actual native mesh owner")
    check=(case/"mesh-check.log").read_text()
    if "Mesh OK." not in check or not re.search(
       rf"(?m)^\s*cells:\s*{n}\s*$",check):
        fail(f"{mode} not native checkMesh Mesh OK. on {n} cells")
    pts=(case/"constant/polyMesh/points").read_text()
    pointrows=[tuple(map(float,r)) for r in VEC.findall(pts)]
    expected=(nx+1)*(ny+1)*(nz+1)
    if len(pointrows)!=expected:fail(f"{mode} native polyMesh points mismatch")
    moved=0
    for x,y,z in pointrows:
        base_x=round(x*nx)/nx
        base_y=round(y*ny)/ny
        base_z=round(z*nz)/nz
        if abs(y-base_y)>1e-10 or abs(z-base_z)>1e-10:
            fail("physical mesh y/z vertices do not match frozen grid")
        edge=any(abs(q)<1e-10 or abs(q-1)<1e-10 for q in (base_x,base_y,base_z))
        dx=0 if edge else .04*math.sin(math.pi*base_x)*math.sin(math.pi*base_y)*math.sin(math.pi*base_z)
        if abs(x-base_x-dx)>2e-10:fail(f"{mode} frozen warped polyMesh vertex not reproduced")
        if abs(dx)>1e-12:moved+=1
    if moved<50:fail("native 3D interior warp missing")
    centers=field(case/"0/C",n,"vector","0 1 0 0 0 0 0")
    vols=field(case/"0/V",n,"scalar","0 3 0 0 0 0 0")
    if min(vols)<=0 or abs(sum(vols)-1)>1e-8:
        fail("native physical volume mismatch")
    iniU=field(case/"0/U",n,"vector","0 1 -1 0 0 0 0")
    if max(abs(a-b) for (x,y,z),u in zip(centers,iniU)
           for a,b in zip(u,(.05*math.sin(math.pi*y)*math.sin(math.pi*z),0,0)))>1e-12:
        fail("physical native initial PISO velocity mismatch")
    props=(case/"constant/transportProperties").read_text()
    for literal in ("momentumMode pisoLaminar;","auditActivityFlux true;",
        "sourceShape gaussian;","sourceCenter [0 1 0 0 0 0 0] (0.38 0.62 0.55);",
        "sourceSigma [0 1 0 0 0 0 0] 0.13;",
        "sourceHT_Bq_s [0 0 -2 0 0 0 0] 60000;",
        "sourceHTO_Bq_s [0 0 -2 0 0 0 0] 12000;",
        "sourceStart [0 0 1 0 0 0 0] 0.075;",
        "sourceEnd [0 0 1 0 0 0 0] 0.225;",
        "kOx [0 0 -1 0 0 0 0] 0.4;",
        "lambdaT [0 0 -1 0 0 0 0] 0.2;",
        "DHT [0 2 -1 0 0 0 0] 0.01;",
        "DHTO [0 2 -1 0 0 0 0] 0.02;",
        "nu [0 2 -1 0 0 0 0] 0.05;"):
        if literal not in props:fail(f"{mode} changed fixed physical input {literal}")
    ctrl=(case/"system/controlDict").read_text()
    delta=re.search(r"(?m)^\s*deltaT\s+("+VNUM+r")\s*;",ctrl)
    if not delta or abs(float(delta.group(1))-dt)>1e-14:
        fail("configured time step unexpectedly changed")
    if f"writeInterval {steps};" not in ctrl or "endTime 0.4;" not in ctrl:
        fail("wrong final write time/end clock")
    pbc=(case/"0/p").read_text()
    for patch,val in (("xMin",.1),("xMax",0.)):
        pm=re.search(rf"\b{patch}\s*\{{\s*type fixedValue;\s*value uniform ("+VNUM+r");",pbc)
        if not pm or abs(float(pm.group(1))-val)>1e-14:
            fail(f"{mode} pressure boundary {patch} changed")
    ubc=(case/"0/U").read_text()
    for patch in ("yMin","yMax","zMin","zMax"):
        if not re.search(rf"\b{patch}\s*\{{\s*type fixedValue;\s*value uniform \(0 0 0\);",ubc):
            fail(f"{mode} lateral wall no-slip {patch} changed")
    for sp in ("HT","HTO"):
        bc=(case/"0"/sp).read_text()
        if "internalField uniform 0;" not in bc:fail("nonzero initial species")
        for patch in ("xMin","xMax"):
            if not re.search(rf"\b{patch}\s*\{{\s*type inletOutlet;\s*inletValue uniform 0;",bc):
                fail(f"external {sp} activity injected at {patch}")
    s=(case/"solve.log").read_text()
    tt=[float(x) for x in re.findall(r"(?m)^Time = ("+VNUM+r")\s*$",s)]
    if len(tt)!=steps or any(abs(t-(i+1)*dt)>1e-10 for i,t in enumerate(tt)):
        fail(f"{mode} actual time step count/clock !={steps}*{dt}")
    if len(re.findall("Solving for p,",s))!=9*steps or any(
       len(re.findall("Solving for "+sp+",",s))!=steps for sp in ("HT","HTO")):
        fail(f"{mode} native PISO/HT/HTO equations were not actually solved")
    inv=rows(s,INV,"inventory",steps)
    em=rows(s,EM,"emission",steps)
    ht=rows(s,HTF,"HT flux",steps)
    zt=rows(s,ZTF,"HTO flux",steps)
    flow=rows(s,FLOW,"pressure boundary",steps)
    div=rows(s,DIV,"pressure continuity",steps)
    oldh=oldz=0.; netexp=[0.,0.]; qcum=[0.,0.]; maxbudget=maxflow=maxsource=0.
    for i,(t,(h,z),(jh,jz,tau),(fl,fr,fn,rbh),
           (gl,gr,gn,rbz,rbt),(qi,qo),(dv,bnv,ke)) in enumerate(
         zip(tt,inv,em,ht,zt,flow,div)):
        overlap=max(0,min(t,.225)-max(t-dt,.075))
        expectedh,expectedz=60000*overlap,12000*overlap
        maxsource=max(maxsource,abs(jh-expectedh),abs(jz-expectedz),abs(tau-overlap))
        qcum[0]+=expectedh;qcum[1]+=expectedz
        bh=h-oldh+dt*fn+dt*.6*h-expectedh
        bz=z-oldz+dt*gn+dt*.2*z-dt*.4*h-expectedz
        bt=h+z-oldh-oldz+dt*(fn+gn)+dt*.2*(h+z)-expectedh-expectedz
        maxbudget=max(maxbudget,*map(abs,(bh,bz,bt,rbh,rbz,rbt,
                                         bh-rbh,bz-rbz,bt-rbt)))
        maxflow=max(maxflow,abs(qi+qo),abs(dv),abs(bnv))
        if max(abs(fn-fl-fr),abs(gn-gl-gr))>1e-3:fail("nonzero side-wall species flux")
        if qi>-.003 or qo<.003:fail("pressure-inlet flow direction/magnitude invalid")
        if min(h,z)<-1e-8:fail("negative activity inventory")
        netexp[0]+=dt*fn;netexp[1]+=dt*gn
        oldh,oldz=h,z
    if maxsource>1e-9 or abs(qcum[0]-9000)>1e-9 or abs(qcum[1]-1800)>1e-9:
        fail(f"{mode} wrong native exact-window activity Bq source")
    if maxbudget>1e-4 or maxflow>1e-7:
        fail(f"{mode} native budget/continuity failed {maxbudget} {maxflow}")
    out={"name":mode,"cells":n,"nx":nx,"ny":ny,"nz":nz,"dt_s":dt,
         "steps":steps,"warped_vertices":moved,"HT_released_Bq":qcum[0],
         "HTO_released_Bq":qcum[1],"max_Bq_budget_residual":maxbudget,
         "max_m3_s_pressure_residual":maxflow,
         "HT_net_export_Bq":netexp[0],"HTO_net_export_Bq":netexp[1]}
    for sp,j in (("HT",0),("HTO",1)):
        vals=field(case/"0.4"/sp,n,"scalar","0 -3 -1 0 0 0 0")
        if min(vals)<-1e-5:fail("negative final physical cells")
        integrated=sum(v*x for v,x in zip(vols,vals))
        if abs(integrated-inv[-1][j])>1e-4:
            fail(f"{mode} native final activity does not match solver logs {sp}")
        out[f"{sp}_inventory_Bq"]=integrated
        for name,r in RECEPTORS.items():
            weights=[v*math.exp(-sum((a-b)**2 for a,b in zip(p,r))/.02)
                     for v,p in zip(vols,centers)]
            normalizer=sum(weights)
            if normalizer<=0:fail("receptor has zero support")
            conc=sum(w*x for w,x in zip(weights,vals))/normalizer
            nearest=min(math.dist(p,r) for p in centers)
            if not math.isfinite(conc) or conc<0 or nearest>.17:
                fail("receptor concentration bad/physical centre missing")
            out[f"{name}_{sp}_Bq_m3"]=conc
            out[f"{name}_nearest_cell_m"]=nearest
    u=field(case/"0.4/U",n,"vector","0 1 -1 0 0 0 0")
    if max(abs(x) for row in u for x in row)<.001:fail("pressure flow disappeared")
    wall=(case/"wall-seconds.txt").read_text().strip()
    out["solver_wall_s"]=float(wall)
    if not (0<out["solver_wall_s"]<1000):fail("missing bounded native wall time")
    print(f"PASS: M1C-8 {mode}: {n} cells {steps} steps, source=9000/1800 Bq, "
          f"maxBqBudget={maxbudget:.12g}, maxFlow={maxflow:.12g}")
    return out

def analyze(root):
    data={mode:get(root/mode,mode) for mode in CASES}
    columns=[f"R{r}_{sp}_Bq_m3" for r in range(1,4) for sp in ("HT","HTO")]
    columns+=["HT_inventory_Bq","HTO_inventory_Bq","HT_net_export_Bq","HTO_net_export_Bq"]
    comparison={}
    flags=[]
    for axis,sequence in (("space",("space_lo","baseline","space_hi")),
                          ("time",("time_lo","baseline","time_hi"))):
        lo,mid,hi=(data[k] for k in sequence)
        checks={}
        for key in columns:
            a,b,c=(item[key] for item in (lo,mid,hi))
            dlo=abs(b-a);dhi=abs(c-b)
            improved=dhi<dlo and dlo>1e-9
            checks[key]=dict(coarse=a,mid=b,fine=c,change_lo_mid=dlo,
                             change_mid_hi=dhi,ratio=dhi/max(dlo,1e-9),
                             improved=improved,
                             fine_relative_change=dhi/max(abs(c),1e-9))
            if key in ("R1_HT_Bq_m3","R1_HTO_Bq_m3"):
                flags.append(improved)
        comparison[axis]=checks
    convergence="OBSERVED_LIMITED_TREND" if all(flags) else "CONVERGENCE_NOT_DEMONSTRATED"
    report={"physical_setup":"M1C-8 fixed 1m PISO warped duct, not real atmospheric plume",
            "runs":data,"comparisons":comparison,"convergence_assessment":convergence}
    (root/"convergence-report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    for axis in comparison:
        for sp in ("HT","HTO"):
            key="R1_"+sp+"_Bq_m3";m=comparison[axis][key]
            print(f"CONVERGENCE {axis} {sp} R1 L={m['coarse']:.9g} M={m['mid']:.9g} "
                  f"H={m['fine']:.9g} lowDelta={m['change_lo_mid']:.9g} "
                  f"highDelta={m['change_mid_hi']:.9g} ratio={m['ratio']:.8g} "
                  f"improved={m['improved']}")
    print("M1C-8 NUMERICAL_CONVERGENCE="+convergence)
    # Full native numerical audit metrics, source-time integrations and
    # three fixed receptor concentrations survive ephemeral Docker cleanup.
    print("M1C8_REPORT_JSON="+json.dumps(report,sort_keys=True,separators=(",",":")))
    print("PASS: M1C-8 all 5 native mesh/time conservation-and-receptor audits (convergence separately classified)")
    return convergence

if __name__=="__main__":
    if len(sys.argv)==3 and sys.argv[1]=="--matrix":
        analyze(Path(sys.argv[2]))
    elif len(sys.argv)==4 and sys.argv[1]=="--single":
        get(Path(sys.argv[2]),sys.argv[3])
    else:
        raise SystemExit("usage: check_convergence.py --matrix ROOT | --single CASE MODE")
