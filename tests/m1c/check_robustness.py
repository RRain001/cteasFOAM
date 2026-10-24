#!/usr/bin/env python3
"""Independent M1C-20 multi-axis MPI and time-step source/activity audits."""
import math
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_gaussian_green import (
    field, DIM, INV, HTF, ZTF, EM, START, END, RATES, NUM,
    POS, green, GL32
)
from check_fixed_receptors import grid, weights

def fail(why):
    raise SystemExit("FAIL: M1C20 " + why)

def records(path, count, dt):
    text=(path/'solve.log').read_text()
    times=[float(x) for x in re.findall(r"(?m)^Time = ("+NUM+r")\s*$",text)]
    inventories=[tuple(map(float,x)) for x in INV.findall(text)]
    fluxh=[tuple(map(float,x)) for x in HTF.findall(text)]
    fluxz=[tuple(map(float,x)) for x in ZTF.findall(text)]
    emissions=[tuple(map(float,x)) for x in EM.findall(text)]
    if any(len(v)!=count for v in (times,inventories,fluxh,fluxz,emissions)):
        fail(f"true native 5-way {count}-step records missing at {path}")
    if any(abs(t-i*dt)>2e-11 for i,t in enumerate(times,1)):
        fail("actual clock does not match Δt and end time")
    previous=(0.,0.)
    released=[0.,0.]
    maxbudget=0.
    for i,(inv,fht,fz,em) in enumerate(zip(inventories,fluxh,fluxz,emissions),1):
        h,z=inv
        overlap=max(0.,min(i*dt,END)-max((i-1)*dt,START))
        qh,qz=RATES[0]*overlap,RATES[1]*overlap
        released[0]+=qh
        released[1]+=qz
        if max(abs(em[0]-qh),abs(em[1]-qz),abs(em[2]-overlap))>1e-8:
            fail(f"native actual Gaussian finite source overlap at step {i}")
        lh,rh,fh,bh=fht
        lz,rz,fz,bz,bt=fz
        residualh=h-previous[0]+dt*fh+dt*.6*h-qh
        residualz=z-previous[1]+dt*fz+dt*.2*z-dt*.4*h-qz
        residualt=h+z-sum(previous)+dt*(fh+fz)+dt*.2*(h+z)-qh-qz
        current=max(map(abs,(residualh,residualz,residualt,bh,bz,bt,fh,fz,lh+rh,lz+rz)))
        maxbudget=max(maxbudget,current)
        if not math.isfinite(current) or current>1e-5:
            fail(f"native HT/HTO/total signed activity + cyclic face Bq at step {i}: {current}")
        previous=inv
    if max(abs(released[0]-3840.),abs(released[1]-768.))>1e-8:
        fail("total integrated window source")
    logs=re.findall(r"CTEAS_M1C19 species=(HT|HTO) face_lambda_min=("+NUM+r") face_lambda_max=("+NUM+r") limited_faces=(\d+) activity_min_Bq_m3=("+NUM+r")",text)
    if len(logs)!=2*count:
        fail("native per-species face MULES records incomplete")
    for j,(species,low,high,limited,minimum) in enumerate(logs):
        if species!=("HT" if j%2==0 else "HTO"):
            fail("native HT/HTO face MULES chronology wrong")
        lo,hi,mn=float(low),float(high),float(minimum)
        if not all(math.isfinite(v) for v in (lo,hi,mn)) or not (0<=lo<=hi<=1) or mn < -1e-8 or int(limited)<0:
            fail("invalid native limited face coefficient/cell value")
    return inventories,maxbudget

def mpi_check(serial,parallel,log_path,nprocs,label):
    original=serial/'solve.log'
    for sp in ("HT","HTO"):
        a=field(serial/'0.4'/sp,512,'scalar',DIM)
        b=field(parallel/'0.4'/sp,512,'scalar',DIM)
        delta=max(abs(x-y) for x,y in zip(a,b))
        print(f"M1C20_MPI dir={label} ranks={nprocs} species={sp} max_abs_Bq_m3={delta:.12g}")
        if delta>3e-4:
            fail("native MPI decomposed reconstructed-cell HT/HTO parity")
    # The true MPI log (not reconstructed serial metadata) must contain every step.
    text=log_path.read_text()
    if f"nProcs : {nprocs}" not in text:
        fail("real MPI rank count absent")
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as t:
        p=Path(t);(p/'solve.log').write_text(text)
        mpi,budget=records(p,80,.005)
    baseline,_=records(serial,80,.005)
    for step,(a,b) in enumerate(zip(baseline,mpi),1):
        if max(abs(x-y) for x,y in zip(a,b))>1e-5:
            fail(f"native MPI activity parity failed step={step}")
    print(f"PASS: M1C20 native MPI {label} ranks={nprocs} 80-step signed Bq and reconstructed field parity")

def temporal(case,count,dt,label):
    inv,budget=records(case,count,dt)
    centres,volumes,lookup=grid(case,8)
    results={}
    for sp in ("HT","HTO"):
        actual=field(case/'0.4'/sp,512,'scalar',DIM)
        if min(actual)<-1e-8 or not all(math.isfinite(x) for x in actual):
            fail(f"true final native species values invalid dt={dt}")
        j=0 if sp=="HT" else 1
        if abs(sum(x*y for x,y in zip(volumes,actual))-inv[-1][j])>1e-5:
            fail("field vs logged final integrated Bq mismatch")
        results[sp]=actual
    for probe,position in POS.items():
        _,interpolate=weights(position,centres,lookup,8)
        analytic=green(position,GL32)
        for j,sp in enumerate(("HT","HTO")):
            native=interpolate(results[sp])
            error=abs(native-analytic[j])/analytic[j]
            print(f"M1C20_TEMPORAL dt={dt} probe={probe} species={sp} relative_green_error={error:.12g}")
            if not math.isfinite(error) or error>1.5:
                fail(f"native R1/R2/R3 physical Green receptor threshold dt={dt}")
    print(f"PASS: M1C20 native temporal {label} {count} steps source/Bq/Green/face/positive")

if __name__=="__main__":
    args=sys.argv[1:]
    if len(args)==5 and args[0]=="mpi":
        mpi_check(Path(args[1]),Path(args[2]),Path(args[3]),int(args[4].split(":")[0]),args[4].split(":")[1])
    elif len(args)==4 and args[0]=="time":
        temporal(Path(args[1]),int(args[2]),float(args[3]),str(args[3]))
    else:
        raise SystemExit("usage: check_robustness.py mpi SERIAL PARALLEL LOG NUMPROCS:LABEL | time CASE STEPS DT")
