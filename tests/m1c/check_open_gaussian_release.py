#!/usr/bin/env python3
"""Independent physical-cell and time-resolved signed Bq budget audit for M1C-7."""
from pathlib import Path
import math
import re
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"m1b1"))
from check_warped_open_pressure import audit_physical_points
sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_gaussian_source import field as native_field

NUM=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
INV=re.compile(r"CTEAS_M1C activity_HT_Bq=("+NUM+r") activity_HTO_Bq=("+NUM+r")")
EM=re.compile(r"CTEAS_M1C emitted_HT_Bq=("+NUM+r") emitted_HTO_Bq=("+NUM+r") source_overlap_s=("+NUM+r")")
HTF=re.compile(r"CTEAS_M1C flux_HT_xMin_Bq_s=("+NUM+r") flux_HT_xMax_Bq_s=("+NUM+r") flux_HT_net_Bq_s=("+NUM+r") flux_HT_budget_residual_Bq=("+NUM+r")")
ZTF=re.compile(r"CTEAS_M1C flux_HTO_xMin_Bq_s=("+NUM+r") flux_HTO_xMax_Bq_s=("+NUM+r") flux_HTO_net_Bq_s=("+NUM+r") flux_HTO_budget_residual_Bq=("+NUM+r") flux_total_budget_residual_Bq=("+NUM+r")")
FLOW=re.compile(r"CTEAS_M1B1 pressure_inlet_flux=("+NUM+r") pressure_outlet_flux=("+NUM+r")")
DIV=re.compile(r"CTEAS_M1B1 max_div_phi=("+NUM+r") boundary_net_flux=("+NUM+r") kinetic_energy=("+NUM+r")")

def fail(msg):
    raise SystemExit("FAIL: M1C-7 "+msg)

def records(text,pattern,name):
    rows=[tuple(map(float,row)) for row in pattern.findall(text)]
    if len(rows)!=20 or not all(math.isfinite(v) for row in rows for v in row):
        fail(f"{name}: not 20 finite native solver rows")
    return rows

def check(case, log, mode, tag):
    if mode not in ("release","zero_rate") or tag not in ("serial","mpi"):
        fail("mode or tag invalid")
    case=Path(case)
    log=Path(log)
    if not log.is_file():fail("missing actual native solver log")
    s=log.read_text()
    if tag=="mpi" and not re.search(r"nProcs\s*:\s*2",s):
        fail("claimed MPI did not use two processes")
    for p in ("constant/polyMesh/owner","constant/polyMesh/points","0/C","0/V"):
        if not (case/p).is_file():fail(f"native mesh missing {p}")
    quality=(case/"mesh-check.log").read_text()
    if "Mesh OK." not in quality or not re.search(r"(?m)^\s*cells:\s*960\s*$",quality):
        fail("not checkMesh-valid actual 960-cell native mesh")
    if audit_physical_points(case,True)!=665:
        fail("not 665 physically warped 3-D interior vertices")
    props=(case/"constant/transportProperties").read_text()
    qh,qz=((60000,12000) if mode=="release" else (0,0))
    for token in (
        "momentumMode pisoLaminar;","auditActivityFlux true;",
        "sourceShape gaussian;","sourceCenter [0 1 0 0 0 0 0] (0.38 0.62 0.55);",
        "sourceSigma [0 1 0 0 0 0 0] 0.13;",
        f"sourceHT_Bq_s [0 0 -2 0 0 0 0] {qh};",
        f"sourceHTO_Bq_s [0 0 -2 0 0 0 0] {qz};",
        "sourceStart [0 0 1 0 0 0 0] 0.075;",
        "sourceEnd [0 0 1 0 0 0 0] 0.225;",
        "kOx [0 0 -1 0 0 0 0] 0.4;",
        "lambdaT [0 0 -1 0 0 0 0] 0.2;",
        "DHT [0 2 -1 0 0 0 0] 0.01;",
        "DHTO [0 2 -1 0 0 0 0] 0.02;",
        "nu [0 2 -1 0 0 0 0] 0.05;",
    ):
        if token not in props:fail(f"wrong source/flow input {token}")
    if "fluxRequired { default no; p; HT; HTO; }" not in (case/"system/fvSchemes").read_text():
        fail("FVM species fluxRequired missing")
    for species in ("HT","HTO"):
        f=(case/"0"/species).read_text()
        for patch in ("xMin","xMax"):
            if not re.search(rf"\b{patch}\s*\{{\s*type inletOutlet;\s*inletValue uniform 0;",f):
                fail(f"species source case has external boundary inlet {patch}/{species}")
    pressure=(case/"0/p").read_text()
    pm = re.search(r"\bxMin\s*\{\s*type fixedValue;\s*value uniform ("+NUM+r");",pressure)
    if not pm or abs(float(pm.group(1))-.1)>1e-14:
        fail("pressure source not genuinely driven")
    times=[float(t) for t in re.findall(r"(?m)^Time = ("+NUM+r")\s*$",s)]
    if len(times)!=20 or any(abs(t-.02*i)>2e-12 for i,t in enumerate(times,1)):
        fail("wrong 20-step Euler physical clock")
    if (len(re.findall(r"Solving for p,",s))!=180 or
        len(re.findall(r"Solving for HT,",s))!=20 or
        len(re.findall(r"Solving for HTO,",s))!=20 or
        "mode=pisoLaminar" not in s):
        fail("native 180 PISO pressure/20 HT/HTO implicit solves missing")
    inv=records(s,INV,"inventory")
    em=records(s,EM,"emitted")
    ht=records(s,HTF,"htflux")
    zt=records(s,ZTF,"htoflux")
    flow=records(s,FLOW,"pressure flux")
    div=records(s,DIV,"continuity")
    # Recompute actual cell volume normalized Gaussian from native cell centres.
    volumes=native_field(case/"0/V",960,"scalar","0 3 0 0 0 0 0")
    centers=native_field(case/"0/C",960,"vector","0 1 0 0 0 0 0")
    if min(volumes)<=0 or abs(sum(volumes)-1)>1e-8:
        fail("invalid actual native polyMesh cell volumes")
    initialU=native_field(case/"0/U",960,"vector","0 1 -1 0 0 0 0")
    initial_u_error=max(
        abs(actual-predicted)
        for pos,u in zip(centers,initialU)
        for actual,predicted in zip(
            u,(0.05*math.sin(math.pi*pos[1])*math.sin(math.pi*pos[2]),0.,0.)
        )
    )
    if initial_u_error>1e-12:
        fail(f"native startup U not frozen analytic divergence-free profile {initial_u_error}")
    center=(.38,.62,.55); sigma=.13
    weights=[math.exp(-sum((a-b)**2 for a,b in zip(p,center))/(2*sigma*sigma)) for p in centers]
    weighted=sum(v*w for v,w in zip(volumes,weights))
    if weighted<=0 or not math.isfinite(weighted) or max(weights)/min(weights)<1.5:
        fail("invalid/non-local Gaussian source")
    moment=[sum(v*w*c[i] for v,w,c in zip(volumes,weights,centers))/weighted for i in range(3)]
    momenterr=math.sqrt(sum((a-b)**2 for a,b in zip(moment,center)))
    if momenterr>.15:fail("Gaussian source centroid not local")
    oldh=oldz=0.
    releasedh=releasedz=0.
    maxbudget=max_logged=maxsource=maxflow=0.
    for step,(t,(h,z),(jh,jz,source_s),(hl,hr,hn,hres),(zl,zr,zn,zres,tres),
               (qin,qout),(maxdiv,netvol,ke)) in enumerate(
                  zip(times,inv,em,ht,zt,flow,div),1):
        tau=max(0.,min(t,.225)-max(t-.02,.075))
        sourceh,sourcez=qh*tau,qz*tau
        releasedh+=sourceh; releasedz+=sourcez
        maxsource=max(maxsource,abs(jh-sourceh),abs(jz-sourcez),abs(source_s-tau))
        bh=h-oldh+.02*hn+.02*.6*h-sourceh
        bz=z-oldz+.02*zn+.02*.2*z-.02*.4*h-sourcez
        bt=h+z-oldh-oldz+.02*(hn+zn)+.02*.2*(h+z)-sourceh-sourcez
        maxbudget=max(maxbudget,*map(abs,(bh,bz,bt)))
        max_logged=max(max_logged,*map(abs,(hres,zres,tres,hres-bh,zres-bz,tres-bt)))
        maxflow=max(maxflow,abs(qin+qout),abs(netvol),abs(maxdiv))
        if max(abs(hn-hl-hr),abs(zn-zl-zr))>1e-3:
            fail(f"nonzero insulated-wall species flux step {step}")
        if max(abs(bh),abs(bz),abs(bt),abs(hres),abs(zres),abs(tres))>1e-4:
            fail(f"signed activity residual step {step}: {bh},{bz},{bt}")
        if max(abs(qin+qout),abs(netvol),abs(maxdiv))>1e-7:
            fail(f"native PISO volume conservation failed step {step}")
        if qin>-.003 or qout<.003:
            fail(f"physical pressure forward flow insufficient/sign wrong step {step}")
        if min(h,z)<-1e-8:
            fail("negative integrated activity")
        oldh,oldz=h,z
    if maxsource>1e-9:
        fail(f"wrong time overlap emission {maxsource}")
    if abs(releasedh-qh*.15)>1e-9 or abs(releasedz-qz*.15)>1e-9:
        fail("integrated source Bq wrong")
    out={}
    for sp,idx in (("HT",0),("HTO",1)):
        initial=native_field(case/"0"/sp,960,"scalar","0 -3 -1 0 0 0 0")
        final=native_field(case/"0.4"/sp,960,"scalar","0 -3 -1 0 0 0 0")
        if max(map(abs,initial))>1e-12 or min(final)<-1e-5:
            fail(f"wrong HT/HTO initial state or negative final {sp}")
        reconstructed=sum(a*v for a,v in zip(final,volumes))
        if abs(reconstructed-inv[-1][idx])>1e-4:
            fail(f"native activity field vs inventory log {sp}")
        if mode=="release":
            if reconstructed<=1 or max(final)<=1:
                fail(f"localized release did not create physical {sp} activity")
        else:
            if max(map(abs,final))>1e-8 or max(abs(row[idx]) for row in inv)>1e-8:
                fail(f"zero-release created {sp} species")
        out[sp]=final
    U=native_field(case/"0.4/U",960,"vector","0 1 -1 0 0 0 0")
    maxU=max(abs(x) for vec in U for x in vec)
    if maxU < .001:fail("not a resolved pressure driven velocity")
    # Nearest native physical receptor in positive-x direction of source.
    d2=[sum((p[i]-([.60,.62,.55][i]))**2 for i in range(3)) for p in centers]
    receptor=min(range(len(d2)),key=d2.__getitem__)
    if d2[receptor]>.02 or (mode=="release" and
        (out["HT"][receptor]<=1e-9 or out["HTO"][receptor]<=1e-9)):
        fail("downstream native 3-D receptor not sampled/positive")
    if mode=="zero_rate":
        # The configured source window can have positive overlap seconds
        # despite zero Bq/s rates; overlap time is NOT injected activity.
        if max(abs(v) for jh,jz,tau in em for v in (jh,jz))>1e-8:
            fail("zero-release generated nonzero emitted activity Bq")
        for rows in (ht,zt):
            if max(abs(v) for row in rows for v in row)>1e-8:
                fail("zero-release generated nonzero species boundary flux")
    print(f"PASS: M1C-7 {mode}/{tag} cells=960 Bq_HT_injected={releasedh:.12g} "
          f"Bq_HTO_injected={releasedz:.12g} max_budget_Bq={maxbudget:.12g} "
          f"max_logged_Bq={max_logged:.12g} max_flow_m3_s={maxflow:.12g} "
          f"centroid_drift_m={momenterr:.12g} downstream_HT_Bq_m3={out['HT'][receptor]:.12g} "
          f"downstream_HTO_Bq_m3={out['HTO'][receptor]:.12g}")
    return dict(inv=inv,ht=ht,zt=zt,flow=flow,U=U,HT=out["HT"],HTO=out["HTO"])

def parity(a,b):
    for key,lim in (("inv",20),("ht",25),("zt",25),("flow",1e-6)):
        delta=max(abs(x-y) for ra,rb in zip(a[key],b[key]) for x,y in zip(ra,rb))
        if delta>lim:fail(f"MPI parity {key} {delta} > {lim}")
        print(f"M1C-7 MPI {key} maxdelta={delta:.12g}")
    for sp,linf,rms in (("HT",100,30),("HTO",20,6),("U",.003,.003)):
        va=a[sp];vb=b[sp]
        if sp=="U":
            delta=[x-y for xrow,yrow in zip(va,vb) for x,y in zip(xrow,yrow)]
        else: delta=[x-y for x,y in zip(va,vb)]
        err=max(map(abs,delta))
        root=math.sqrt(sum(x*x for x in delta)/len(delta))
        if err>linf or root>rms:fail(f"MPI parity {sp} L∞={err} RMS={root}")
        print(f"M1C-7 MPI {sp} Linf={err:.12g} RMS={root:.12g}")
    print("PASS: M1C-7 source/PISO actual two-rank MPI native field and budget parity")

if __name__=="__main__":
    if len(sys.argv)==5:
        check(*sys.argv[1:])
    elif len(sys.argv)==9:
        a=check(*sys.argv[1:5]);b=check(*sys.argv[5:9]);parity(a,b)
    else:
        raise SystemExit("usage: check_open_gaussian_release.py CASE LOG release|zero_rate serial|mpi [CASE LOG release|zero_rate serial|mpi]")
