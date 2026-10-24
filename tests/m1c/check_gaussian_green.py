#!/usr/bin/env python3
"""Independent M1C-13 periodic Gaussian Green-kernel time convolution; no numpy."""
from pathlib import Path
import sys
import json
import math
import re

sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_3d_analytic import field,physical_geometry,INV,HTF,ZTF,DIM
from check_gaussian_source import EM

SOURCE_CENTER=(.42,.55,.50)
SIGMA=.08
WIND=(1.,.05,.03)
START=.073
END=.137
T=.4
D_HT=.01
D_HTO=.02
OXID=.4
DECAY=.2
RATES=(60000.,12000.)
POS={"R1":(.60,.62,.55),"R2":(.75,.50,.50),"R3":(.88,.62,.55)}
NUM=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
GRIDS={(8,"normal"):80,(12,"normal"):80,(16,"normal"):80,(8,"coarse"):40}
CYCLE=("xMin","xMax","yMin","yMax","zMin","zMax")
DIMS="0 -3 -1 0 0 0 0"

def fail(s):raise SystemExit("FAIL: M1C-13 "+s)

def gl_rule(n):
    """Independent Legendre roots and positive weights on [0,1], no scipy."""
    roots=[]
    weights=[]
    for j in range(1,n+1):
        z=math.cos(math.pi*(j-.25)/(n+.5))
        for _ in range(50):
            p0,p1=1.,z
            for k in range(2,n+1):
                p0,p1=p1,((2*k-1)*z*p1-(k-1)*p0)/k
            pp=n*(z*p1-p0)/(z*z-1)
            nz=z-p1/pp
            if abs(nz-z)<3e-15: z=nz;break
            z=nz
        roots.append((z+1)/2)
        weights.append(1/((1-z*z)*pp*pp))
    pairs=sorted(zip(roots,weights))
    if abs(sum(weights)-1)>1e-12 or min(weights)<=0:fail("independent quadrature invalid")
    return pairs

GL24=gl_rule(24)
GL32=gl_rule(32)

def gaussian_wrapped(x,age,variance,images):
    if not 0<variance<1:fail("nonpositive/huge Gaussian Green variance")
    factor=1/math.sqrt(2*math.pi*variance)
    p=1.
    for coord,c,u in zip(x,SOURCE_CENTER,WIND):
        d=coord-c-u*age
        summation=sum(math.exp(-(d+j)**2/(2*variance))
                      for j in range(-images,images+1))
        p*=factor*summation
    return p

def green(x,quad,images=2):
    ht=hto_direct=hto_converted=0.
    duration=END-START
    for t_weight,wt in quad:
        s=START+duration*t_weight
        age=T-s
        weighted=duration*wt
        variance_ht=SIGMA**2+2*D_HT*age
        variance_z=SIGMA**2+2*D_HTO*age
        ht+=weighted*math.exp(-(OXID+DECAY)*age)*gaussian_wrapped(x,age,variance_ht,images)
        hto_direct+=weighted*math.exp(-DECAY*age)*gaussian_wrapped(x,age,variance_z,images)
        inner=0.
        for frac,wq in quad:
            r=age*frac
            var=SIGMA**2+2*D_HT*r+2*D_HTO*(age-r)
            inner+=wq*math.exp(-OXID*r)*gaussian_wrapped(x,age,var,images)
        hto_converted+=weighted*math.exp(-DECAY*age)*age*inner
    return RATES[0]*ht, RATES[1]*hto_direct+RATES[0]*OXID*hto_converted

def tail_mass():
    probability=1.
    for c in SOURCE_CENTER:
        p=.5*(math.erf((1-c)/(math.sqrt(2)*SIGMA))
             -math.erf(-c/(math.sqrt(2)*SIGMA)))
        probability*=p
    lost=1-probability
    if not 0<lost<2e-7:fail(f"reference wrapped-vs-truncated tail unexpectedly large {lost}")
    return lost

def nearest(centres,p):
    i=min(range(len(centres)),key=lambda j:sum((v-ref)**2 for v,ref in zip(centres[j],p)))
    distance=math.dist(centres[i],p)
    return i,distance

def check(case,nx,level):
    if (nx,level) not in GRIDS:fail("unknown frozen grid/clock")
    n=nx**3;count=GRIDS[nx,level]
    dt=.005 if level=="normal" else .01
    centres,volumes,moved=physical_geometry(case,n,nx)
    props=(case/"constant/transportProperties").read_text()
    for literal in (
        "sourceShape gaussian;","sourceCenter [0 1 0 0 0 0 0] (0.42 0.55 0.5);",
        "sourceSigma [0 1 0 0 0 0 0] 0.08;",
        "sourceHT_Bq_s [0 0 -2 0 0 0 0] 60000;",
        "sourceHTO_Bq_s [0 0 -2 0 0 0 0] 12000;",
        "sourceStart [0 0 1 0 0 0 0] 0.073;",
        "sourceEnd [0 0 1 0 0 0 0] 0.137;",
        "kOx [0 0 -1 0 0 0 0] 0.4;",
        "lambdaT [0 0 -1 0 0 0 0] 0.2;",
        "DHT [0 2 -1 0 0 0 0] 0.01;",
        "DHTO [0 2 -1 0 0 0 0] 0.02;",
        "momentumMode prescribed;","activityMode htHto;","auditActivityFlux true;"
    ):
        if literal not in props:fail("physical Gaussian source/clock config changed "+literal)
    ctrl=(case/"system/controlDict").read_text()
    for lit in ("endTime 0.4;",f"deltaT {dt};",f"writeInterval {count};"):
        if lit not in ctrl:fail("native clock changed "+lit)
    for p in ("0/HT","0/HTO","0/U","0.4/HT","0.4/HTO","0.4/U"):
        t=(case/p).read_text()
        for patch in CYCLE:
            if not re.search(rf"\b{patch}\s*\{{\s*type\s+cyclic;",t):
                fail("native field missing periodic patch "+p+"/"+patch)
    u0=field(case/"0/U",n,"vector","0 1 -1 0 0 0 0")
    u1=field(case/"0.4/U",n,"vector","0 1 -1 0 0 0 0")
    if max(abs(v-ref) for rows in (u0,u1) for row in rows for v,ref in zip(row,WIND))>1e-12:
        fail("native actual prescribed vector wind differs")
    species=["HT","HTO"]
    ini={s:field(case/"0"/s,n,"scalar",DIMS) for s in species}
    fin={s:field(case/"0.4"/s,n,"scalar",DIMS) for s in species}
    for s in species:
        if max(map(abs,ini[s]))>1e-10 or min(fin[s])<-1e-6:fail("nonzero initial or negative final native source field")
    lines=(case/"solve.log").read_text()
    times=[float(v) for v in re.findall(r"(?m)^Time = ("+NUM+r")\s*$",lines)]
    inv=[tuple(map(float,m)) for m in INV.findall(lines)]
    htflux=[tuple(map(float,m)) for m in HTF.findall(lines)]
    zflux=[tuple(map(float,m)) for m in ZTF.findall(lines)]
    em=[tuple(map(float,m)) for m in EM.findall(lines)]
    if any(len(v)!=count for v in (times,inv,htflux,zflux,em)):
        fail("native all-step Gaussian source/flux/inventory logs incomplete")
    if any(abs(t-i*dt)>2e-12 for i,t in enumerate(times,1)):
        fail("native time grid changed")
    if any(len(re.findall(r"Solving for "+s+r",",lines))!=count for s in ("HT","HTO","c")):
        fail("true native species/passive solver did not run each step")
    if "CTEAS_M1C sourceShape=gaussian" not in lines or "Solving for p," in lines:
        fail("real normalized Gaussian source branch missing / spurious pressure solver")
    previous=[0.,0.]
    source_totals=[0.,0.]
    maxbudget=maxsource=maxflux=0.
    nonzero_steps=0
    for i,((hh,zz),(qh,qz,seconds),fht,fz) in enumerate(zip(inv,em,htflux,zflux),1):
        overlap=max(0.,min(i*dt,END)-max((i-1)*dt,START))
        expecth=RATES[0]*overlap
        expectz=RATES[1]*overlap
        source_totals[0]+=expecth;source_totals[1]+=expectz
        maxsource=max(maxsource,abs(qh-expecth),abs(qz-expectz),abs(seconds-overlap))
        if maxsource>1e-8:fail("native Gaussian source emitted HT/HTO Bq or time overlap corrupted")
        if overlap>0:nonzero_steps+=1
        lh,rh,fluxh,reportedh=fht
        lz,rz,fluxz,reportedz,reportedtot=fz
        if not all(math.isfinite(v) for v in (hh,zz,*fht,*fz)):fail("native inventory nonfinite")
        maxflux=max(maxflux,abs(fluxh),abs(fluxz),abs(lh+rh),abs(lz+rz))
        if maxflux>1e-5:fail("native cyclic species activity flux not balanced")
        bh=hh-previous[0]+dt*fluxh+dt*.6*hh-expecth
        bz=zz-previous[1]+dt*fluxz+dt*.2*zz-dt*.4*hh-expectz
        bt=hh+zz-sum(previous)+dt*(fluxh+fluxz)+dt*.2*(hh+zz)-expecth-expectz
        maxbudget=max(maxbudget,*map(abs,(bh,bz,bt,reportedh,reportedz,reportedtot,
                                          bh-reportedh,bz-reportedz,bt-reportedtot)))
        if maxbudget>1e-5:fail(f"step {i} native Gaussian Bq conservation failed")
        previous=[hh,zz]
    if max(abs(source_totals[0]-3840),abs(source_totals[1]-768))>1e-8:
        fail("whole time-window integrated source not exactly prescribed")
    if nonzero_steps<4 or em[0][0]!=0 or em[-1][0]!=0:
        fail("native physical source timing/spatial pulse incomplete")
    for j,s in enumerate(species):
        actual=sum(a*v for a,v in zip(fin[s],volumes))
        if abs(actual-inv[-1][j])>1e-5:
            fail("final native species Bq inventory does not match real physical field "+s)
    sample={}
    for name,nominal in POS.items():
        idx,dist=nearest(centres,nominal)
        maxdist={8:.16,12:.11,16:.09}[nx]
        if dist>=maxdist:fail("actual warped native receptor cell farther than frozen distance")
        actual=tuple(fin[s][idx] for s in species)
        ref24=green(centres[idx],GL24,2)
        ref32=green(centres[idx],GL32,2)
        ref7=green(centres[idx],GL24,3)
        for a,b,c in zip(ref24,ref32,ref7):
            if abs(a-b)>max(.05,1e-5*abs(b)):
                fail("independent 24-versus-32 time quadrature not converged")
            if abs(a-c)>1e-8:
                fail("independent periodic image tails not converged")
        records={}
        for s,observed,reference in zip(species,actual,ref32):
            if not math.isfinite(reference) or reference<=0:fail("Green-function concentration missing")
            error=observed-reference
            rel=abs(error)/reference
            if reference>=10 and rel>={8:.85,12:.65,16:.50}[nx]:
                fail(f"native {nx} {level} {name} {s} Green relative error={rel:.7g} > frozen limit")
            records[s]={"native_Bq_m3":observed,"green_Bq_m3":reference,
                        "signed_error_Bq_m3":error,"abs_relative_error":rel,
                        "quadrature24_vs32_abs_Bq_m3":abs(ref24[0 if s=="HT" else 1]-reference)}
        sample[name]={"native_centre_m":centres[idx],"native_cell_index":idx,
                      "nominal_receptor_m":nominal,"cell_distance_m":dist,"species":records}
        print(f"M1C13_RECEPTOR N={nx} dt={dt} {name} cell={idx} dist={dist:.8g} "
              f"HTnative={actual[0]:.9g} HTgreen={ref32[0]:.9g} "
              f"HTrel={records['HT']['abs_relative_error']:.7g} "
              f"HTOnative={actual[1]:.9g} HTOgreen={ref32[1]:.9g} "
              f"HTOrel={records['HTO']['abs_relative_error']:.7g}")
    print(f"M1C13_CASE nx={nx} dt={dt} cells={n} trueWarped={moved} "
          f"releasedHT={source_totals[0]:.10g} releasedHTO={source_totals[1]:.10g} "
          f"maxBqResidual={maxbudget:.9g} maxSourceError={maxsource:.9g}")
    return {"nx":nx,"level":level,"cells":n,"steps":count,"dt_s":dt,"source_total_Bq":source_totals,
            "max_budget_Bq":maxbudget,"max_source_error_Bq":maxsource,"lost_source_tail":tail_mass(),
            "receivers":sample}

def matrix(root):
    records={}
    for nx,level in ((8,"coarse"),(8,"normal"),(12,"normal"),(16,"normal")):
        records[f"{nx}_{level}"]=check(root/f"{nx}_{level}",nx,level)
    for name in POS:
        for species in ("HT","HTO"):
            seq=[records[f"{nx}_normal"]["receivers"][name]["species"][species]["abs_relative_error"] for nx in (8,12,16)]
            print(f"M1C13_SPATIAL_DIAGNOSTIC {name} {species} relative_errors_8_12_16={seq} "
                  f"fine_better_than_coarse={seq[-1]<seq[0]}")
    result={"scope":"wrapped full-space periodic Gaussian Green approximation vs truncated volume-renormalized native physical Gaussian; synthetic 1m cube",
            "parent_sha":"424ef085aa085c156c51e52eb8f93b07a955d408","tail_outside_cube":tail_mass(),
            "native_cases":records,"conclusion":"SCOPED_3D_LOCAL_GAUSSIAN_GREEN_DIAGNOSTIC"}
    (root/"m1c13-report.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("M1C13_REPORT_JSON="+json.dumps(result,sort_keys=True,separators=(",",":")))
    print("PASS: M1C-13 four true native 3D Gaussian plume runs, independent periodic Green integral, image/time quadrature, receptor diagnostics and signed HT/HTO Bq")

if __name__=="__main__":
    if len(sys.argv)==5 and sys.argv[1]=="--single":
        check(Path(sys.argv[2]),int(sys.argv[3]),sys.argv[4])
        print("PASS: M1C-13 genuine native Gaussian plume field vs independent 3D Green response")
    elif len(sys.argv)==3 and sys.argv[1]=="--matrix":
        matrix(Path(sys.argv[2]))
    else:raise SystemExit("usage: check_gaussian_green.py --single CASE NX coarse|normal | --matrix ROOT")
