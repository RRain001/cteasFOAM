#!/usr/bin/env python3
"""M1C-12 independently audit time-window/source-free twin 3D FVM fields."""
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_3d_analytic import (
    field, physical_geometry, expected, INV, HTF, ZTF, U, DIM, THRESHOLDS
)
from check_gaussian_source import EM

START=.073
END=.137
QH=60000.
QZ=12000.
END_TIME=.2
STEPS={"coarse":(20,.01),"normal":(40,.005),"fine":(80,.0025)}
MESHES={("8","coarse"),("8","normal"),("8","fine"),("12","normal"),("16","normal")}
REQUIRED={
 "kOx [0 0 -1 0 0 0 0] 0.4;",
 "lambdaT [0 0 -1 0 0 0 0] 0.2;",
 "DHT [0 2 -1 0 0 0 0] 0.01;",
 "DHTO [0 2 -1 0 0 0 0] 0.02;",
 "momentumMode prescribed;",
 "activityMode htHto;",
 "auditActivityFlux true;"
}
def fail(msg):raise SystemExit("FAIL: M1C-12 "+msg)

def analytic_pulse(t):
    def J(rate):
        if t<=START:return 0.
        return (math.exp(-rate*(t-min(t,END)))-math.exp(-rate*(t-START)))/rate
    jh,jz=J(.6),J(.2)
    return (QH*jh, QZ*jz+QH*(jz-jh))

def overlap(i,dt):
    t=i*dt
    return max(0.,min(t,END)-max(t-dt,START))

def backward_euler_pulse(dt,steps):
    h=z=0.
    history=[]
    for i in range(1,steps+1):
        ov=overlap(i,dt)
        h=(h+QH*ov)/(1+.6*dt)
        z=(z+QZ*ov+dt*.4*h)/(1+.2*dt)
        history.append((h,z))
    return history

def audit(case,source,grid,level):
    if source not in ("control","pulse") or (str(grid),level) not in MESHES:
        fail("frozen source or resolution unsupported")
    nt=int(grid);nc=nt**3
    steps,dt=STEPS[level]
    centers,volumes,moved=physical_geometry(case,nc,nt)
    props=(case/"constant/transportProperties").read_text()
    for required in REQUIRED:
        if required not in props:fail(f"changed native physical coefficient {required}")
    source_props={
       "sourceShape uniform;",
       "sourceHT_Bq_s [0 0 -2 0 0 0 0] 60000;",
       "sourceHTO_Bq_s [0 0 -2 0 0 0 0] 12000;",
       "sourceStart [0 0 1 0 0 0 0] 0.072999999999999995;",
       "sourceEnd [0 0 1 0 0 0 0] 0.13700000000000001;"
    }
    # Float output of Python .17g is exactly pinned by generator; no
    # acceptance-driven fuzzy interpretation of an altered window.
    if source=="pulse":
        for t in source_props:
            if t not in props:fail("source window/rate or uniform shape changed: "+t)
    elif any(x in props for x in ("sourceShape","sourceHT_Bq_s","sourceHTO_Bq_s","sourceStart","sourceEnd")):
        fail("source-free matched control configured a source")
    ctl=(case/"system/controlDict").read_text()
    for token in ("endTime 0.2;",f"deltaT {dt};",f"writeInterval {steps};"):
        if token not in ctl:fail("time integration contract changed "+token)
    scheme=(case/"system/fvSchemes").read_text()
    for token in ("div(phi,HT) Gauss upwind;","div(phi,HTO) Gauss upwind;",
                  "laplacianSchemes { default Gauss linear corrected; }"):
        if token not in scheme:fail("FVM scheme changed "+token)
    for fieldname in ("0/U","0/HT","0/HTO","0.2/U","0.2/HT","0.2/HTO"):
        text=(case/fieldname).read_text()
        for patch in ("xMin","xMax","yMin","yMax","zMin","zMax"):
            if not re.search(rf"\b{patch}\s*\{{\s*type\s+cyclic;",text):
                fail("native periodic face absent "+fieldname+"/"+patch)
    for path in ("0/U","0.2/U"):
        native_u=field(case/path,nc,"vector","0 1 -1 0 0 0 0")
        if max(abs(v-ref) for row in native_u for v,ref in zip(row,U["advect"]))>1e-12:
            fail("genuine 3-D prescribed U not matched")
    log=(case/"solve.log").read_text()
    times=[float(v) for v in re.findall(r"(?m)^Time = ([-+0-9.eE]+)\s*$",log)]
    histories=[tuple(map(float,p)) for p in INV.findall(log)]
    fh=[tuple(map(float,p)) for p in HTF.findall(log)]
    fz=[tuple(map(float,p)) for p in ZTF.findall(log)]
    emitted=[tuple(map(float,p)) for p in EM.findall(log)]
    if any(len(rows)!=steps for rows in (times,histories,fh,fz)):
        fail("native time/inventory/face history incomplete")
    if len(emitted)!=(steps if source=="pulse" else 0):
        fail("emission log count inconsistent with physical source")
    if any(abs(t-i*dt)>2e-12 for i,t in enumerate(times,1)):
        fail("actual native timestep clock incorrect")
    if any(len(re.findall(r"Solving for "+name+r",",log))!=steps for name in ("HT","HTO","c")):
        fail("real native 3-D activity/passive FVM solves absent")
    if "Solving for p," in log or "Solving for U" in log or "sourceShape=gaussian" in log:
        fail("unexpected pressure/velocity or Gaussian source branch")
    ini={k:field(case/"0"/k,nc,"scalar",DIM) for k in ("HT","HTO")}
    end={k:field(case/"0.2"/k,nc,"scalar",DIM) for k in ("HT","HTO")}
    ampa=analytic_pulse(END_TIME) if source=="pulse" else (0.,0.)
    exact_errors={}
    for species,j in (("HT",0),("HTO",1)):
        if min(ini[species]+end[species])<=0:fail("non-positive native activity")
        initerr=max(abs(value-expected(c,"advect",species,0))
                    for value,c in zip(ini[species],centers))
        if initerr>1e-6:fail("native initial 3D analytic profile changed")
        errors=[value-expected(c,"advect",species,END_TIME)-ampa[j]
                for value,c in zip(end[species],centers)]
        r=math.sqrt(sum(v*x*x for v,x in zip(volumes,errors))/sum(volumes))
        maximum=max(abs(x) for x in errors)
        lim=THRESHOLDS[nt][j]
        if not (0<r<lim and maximum<lim):
            fail(f"native {source}_{grid}_{level} continuum {species} RMS={r:.12g} Linf={maximum:.12g} > {lim}")
        if abs(sum(a*v for a,v in zip(end[species],volumes))-histories[-1][j])>1e-5:
            fail("native final HT/HTO field Bq integral mismatches solver log")
        exact_errors[species]={"RMS_Bq_m3":r,"Linf_Bq_m3":maximum}
    previous=[sum(x*v for x,v in zip(ini[k],volumes)) for k in ("HT","HTO")]
    source_total_h=source_total_z=0.
    max_budget=max_source=0.
    for i,((hh,zz),htflux,zflux) in enumerate(zip(histories,fh,fz),1):
        hl,hr,net,report_ht=htflux
        zl,zr,netz,report_z,report_total=zflux
        if not all(math.isfinite(x) for x in (hh,zz,*htflux,*zflux)):
            fail("nonfinite logged activity/flux")
        if abs(net)>1e-5 or abs(netz)>1e-5 or abs(hl+hr)>1e-5 or abs(zl+zr)>1e-5:
            fail("cyclic native signed flux inconsistent")
        ov=overlap(i,dt) if source=="pulse" else 0.
        qh,qz=QH*ov,QZ*ov
        source_total_h+=qh;source_total_z+=qz
        if source=="pulse":
            rh,rz,observed_overlap=emitted[i-1]
            max_source=max(max_source,abs(rh-qh),abs(rz-qz),abs(observed_overlap-ov))
            if max_source>1e-8:fail("native emitted HT/HTO Bq or partial-step source overlap corrupted")
        rht=hh-previous[0]+dt*net+dt*.6*hh-qh
        rz=zz-previous[1]+dt*netz+dt*.2*zz-dt*.4*hh-qz
        rt=hh+zz-sum(previous)+dt*(net+netz)+dt*.2*(hh+zz)-qh-qz
        max_budget=max(max_budget,*map(abs,(rht,rz,rt,report_ht,report_z,report_total,
                                            rht-report_ht,rz-report_z,rt-report_total)))
        if max_budget>1e-5:fail(f"native stepwise HT/HTO signed activity Bq budget step {i}")
        previous=[hh,zz]
    if abs(source_total_h-(3840 if source=="pulse" else 0))>1e-8 or \
       abs(source_total_z-(768 if source=="pulse" else 0))>1e-8:
        fail("integrated activity source window did not equal frozen emitted Bq")
    result={"source":source,"N":nt,"level":level,"dt_s":dt,"steps":steps,
            "native_cells":nc,"native_interior_warped_vertices":moved,
            "source_Bq":{"HT":source_total_h,"HTO":source_total_z},
            "max_budget_Bq":max_budget,"max_source_Bq":max_source,
            "continuum_errors":exact_errors}
    print(f"M1C12_CASE {source} N={nt} dt={dt} cells={nc}"
          f" emittedHT={source_total_h:.12g} emittedHTO={source_total_z:.12g}"
          f" HT_RMS={exact_errors['HT']['RMS_Bq_m3']:.10g}"
          f" HTO_RMS={exact_errors['HTO']['RMS_Bq_m3']:.10g}"
          f" max_Bq_budget={max_budget:.9g}")
    return result

def matrix(root):
    results={}
    temporal={}
    for grid,level in ((8,"coarse"),(8,"normal"),(8,"fine"),(12,"normal"),(16,"normal")):
        cases={}
        for typ in ("control","pulse"):
            cases[typ]=audit(root/f"{typ}_{grid}_{level}",typ,grid,level)
        n=grid**3
        native1=root/f"pulse_{grid}_{level}"
        native0=root/f"control_{grid}_{level}"
        steps,dt=STEPS[level]
        direct_recurrence=backward_euler_pulse(dt,steps)
        expected_discrete=direct_recurrence[-1]
        expected_continuum=analytic_pulse(END_TIME)
        pair_metrics={}
        for species,j in (("HT",0),("HTO",1)):
            pulse=field(native1/"0.2"/species,n,"scalar",DIM)
            base=field(native0/"0.2"/species,n,"scalar",DIM)
            additions=[p-b for p,b in zip(pulse,base)]
            true_discrete=expected_discrete[j]
            max_diff=max(abs(a-true_discrete) for a in additions)
            spread=max(additions)-min(additions)
            if max_diff>1e-5 or spread>1e-5:
                fail(f"matched native uniform source {species} differs from independent BE recurrence: maxerr={max_diff} spread={spread}")
            actual_diff=sum(additions)/n
            continuum_error=actual_diff-expected_continuum[j]
            pair_metrics[species]={
                "native_pulse_minus_control_Bq_m3":actual_diff,
                "native_spatial_range_Bq_m3":spread,
                "exact_BE_increment_Bq_m3":true_discrete,
                "exact_continuum_increment_Bq_m3":expected_continuum[j],
                "signed_time_discretization_error_Bq_m3":continuum_error,
                "max_native_BE_error_Bq_m3":max_diff
            }
        results[f"{grid}_{level}"]={"native":cases,"pulse_pair":pair_metrics}
        if grid==8:temporal[level]=pair_metrics
        print(f"M1C12_PAIR N={grid} dt={dt} "
              f"HT_increment={pair_metrics['HT']['native_pulse_minus_control_Bq_m3']:.12g} "
              f"HTO_increment={pair_metrics['HTO']['native_pulse_minus_control_Bq_m3']:.12g} "
              f"max_native_BE={max(v['max_native_BE_error_Bq_m3'] for v in pair_metrics.values()):.9g}")
    for sp in ("HT","HTO"):
        errors=[abs(temporal[k][sp]["signed_time_discretization_error_Bq_m3"])
                for k in ("coarse","normal","fine")]
        ratios=[errors[i+1]/errors[i] for i in range(2)]
        print(f"M1C12_TIME {sp} errors={errors} ratios={ratios}")
        if not all(a>b>0 and b/a<.8 for a,b in zip(errors,errors[1:])):
            fail("pulse-specific 3-step time convergence for "+sp)
    report={"scope":"native uniform finite-window source vs source-free genuine 3D analytic HT/HTO twins, not a Gaussian spatial source validation",
            "reference_parent_sha":"015969ca27034634d954e727ad2aa16051a4c941",
            "source_window_s":[START,END],
            "source_rates_Bq_s":{"HT":QH,"HTO":QZ},
            "pairs":results,
            "time":temporal,
            "conclusion":"SCOPED_3D_UNIFORM_SOURCE_AND_FINITE_PULSE_VERIFICATION"}
    (root/"m1c12-report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("M1C12_REPORT_JSON="+json.dumps(report,sort_keys=True,separators=(",",":")))
    print("PASS: M1C-12 ten actual native periodic 3D FVM runs, matched source-control BE and exact continuum time response, full signed Bq history")

def backward_euler_pulse(dt,steps):
    h=z=0.
    out=[]
    for i in range(1,steps+1):
        ov=overlap(i,dt)
        h=(h+QH*ov)/(1+.6*dt)
        z=(z+QZ*ov+dt*.4*h)/(1+.2*dt)
        out.append((h,z))
    return out

if __name__=="__main__":
    if len(sys.argv)==5 and sys.argv[1]=="--single":
        audit(Path(sys.argv[2]),sys.argv[3],int(sys.argv[4].split(":")[0]),sys.argv[4].split(":")[1])
        print("PASS: M1C-12 native 3D analytic source-control audit")
    elif len(sys.argv)==3 and sys.argv[1]=="--matrix":
        matrix(Path(sys.argv[2]))
    else:
        raise SystemExit("usage: check_pulsed_analytic.py --single CASE control|pulse Nx:coarse|normal|fine | --matrix ROOT")
