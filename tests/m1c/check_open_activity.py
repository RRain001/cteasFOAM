#!/usr/bin/env python3
"""Fail-closed independent parsing of M1C-3 real 3-D HT/HTO inlet transport."""
from pathlib import Path
import math
import re
import sys

# Reuse the previously validated independent geometry verifier, which
# reconstructs all 1323 vertices and checks 665 actual interior warps.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"m1b1"))
from check_warped_open_pressure import audit_physical_points

NUM=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
INV=re.compile(r"CTEAS_M1C activity_HT_Bq=("+NUM+r") activity_HTO_Bq=("+NUM+r")")
HTF=re.compile(r"CTEAS_M1C flux_HT_xMin_Bq_s=("+NUM+r") flux_HT_xMax_Bq_s=("+NUM+r") flux_HT_net_Bq_s=("+NUM+r") flux_HT_budget_residual_Bq=("+NUM+r")")
HTOF=re.compile(r"CTEAS_M1C flux_HTO_xMin_Bq_s=("+NUM+r") flux_HTO_xMax_Bq_s=("+NUM+r") flux_HTO_net_Bq_s=("+NUM+r") flux_HTO_budget_residual_Bq=("+NUM+r") flux_total_budget_residual_Bq=("+NUM+r")")
FLOW=re.compile(r"CTEAS_M1B1 pressure_inlet_flux=("+NUM+r") pressure_outlet_flux=("+NUM+r")")
DIV=re.compile(r"CTEAS_M1B1 max_div_phi=("+NUM+r") boundary_net_flux=("+NUM+r") kinetic_energy=("+NUM+r")")
SCA=re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
VEC=re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
VNUM=re.compile(r"\(\s*("+NUM+r")\s+("+NUM+r")\s+("+NUM+r")\s*\)")

def fail(msg):
    raise SystemExit("FAIL: M1C-3 "+msg)

def field(path,kind="scalar",dims=None):
    s=Path(path).read_text()
    if dims is not None and not re.search(r"dimensions\s+\["+re.escape(dims)+r"\]",s):
        fail(f"wrong dimensions {path}")
    pat = VEC if kind=="vector" else SCA
    m=pat.search(s)
    if m:
        if int(m.group(1))!=960: fail(f"wrong actual 960-cell field {path}")
        vals = ([tuple(map(float,v)) for v in VNUM.findall(m.group(2))]
                if kind=="vector" else [float(v) for v in m.group(2).split()])
    else:
        if kind=="vector":
            m=re.search(r"internalField\s+uniform\s+\(([^()]*)\)\s*;",s)
            vals=[tuple(map(float,m.group(1).split()))]*960 if m else []
        else:
            m=re.search(r"internalField\s+uniform\s+("+NUM+r")\s*;",s)
            vals=[float(m.group(1))]*960 if m else []
    if len(vals)!=960 or not all(math.isfinite(x) for v in vals for x in (v if kind=="vector" else (v,))):
        fail(f"unreadable/nonfinite native field {path}")
    return vals

def audit(case,logfile,mode):
    case=Path(case)
    s=Path(logfile).read_text()
    rows=[[float(x) for x in r] for r in INV.findall(s)]
    ht=[[float(x) for x in r] for r in HTF.findall(s)]
    hto=[[float(x) for x in r] for r in HTOF.findall(s)]
    flow=[[float(x) for x in r] for r in FLOW.findall(s)]
    div=[[float(x) for x in r] for r in DIV.findall(s)]
    if any(len(r)!=20 for r in (rows,ht,hto,flow,div)):
        fail(f"incomplete 20-step flux/inventory/flow log mode={mode}")
    if (len(re.findall(r"(?m)^Time = ",s))!=20 or
        len(re.findall(r"Solving for p,",s))!=180 or
        len(re.findall(r"Solving for HT,",s))!=20 or
        len(re.findall(r"Solving for HTO,",s))!=20 or
        "mode=pisoLaminar" not in s):
        fail(f"not 20 actual PISO HT/HTO matrix solves {mode}")
    if not all(math.isfinite(v) for arr in (rows,ht,hto,flow,div) for r in arr for v in r):
        fail(f"nonfinite log {mode}")
    volumes=field(case/"0/V",dims="0 3 0 0 0 0 0")
    if min(volumes)<=0 or abs(sum(volumes)-1)>1e-8:
        fail(f"invalid actual polyMesh native volumes {mode}")
    init={}
    final={}
    active="xMax" if mode=="reverse" else "xMin"
    for name,value in (("HT",1000000),("HTO",200000)):
        stxt=(case/"0"/name).read_text()
        for patch, val in (("xMin",value if active=="xMin" else 0),
                           ("xMax",value if active=="xMax" else 0)):
            m=re.search(rf"\b{patch}\s*\{{\s*type\s+inletOutlet;\s*inletValue\s+uniform\s+("+NUM+r")\s*;",stxt)
            if not m or abs(float(m.group(1))-val)>1e-12:
                fail(f"wrong real species inflow BC {name}/{patch}/{mode}")
        for wall in ("yMin","yMax","zMin","zMax"):
            if not re.search(rf"\b{wall}\s*\{{\s*type\s+zeroGradient;",stxt):
                fail(f"wrong wall insulation {name}/{wall}/{mode}")
        init[name]=field(case/"0"/name,dims="0 -3 -1 0 0 0 0")
        final[name]=field(case/"0.4"/name,dims="0 -3 -1 0 0 0 0")
        if max(map(abs,init[name]))>1e-12:
            fail(f"not initially empty HT/HTO {name}/{mode}")
        if min(final[name]) < -1e-5 or max(final[name]) > 1.5e6:
            fail(f"physically negative/unbounded activity {name}/{mode}")
        native=sum(v*vol for v,vol in zip(final[name],volumes))
        logged=rows[-1][0 if name=="HT" else 1]
        if abs(native-logged)>1e-4:
            fail(f"native final volume integral and reported mass disagree {name}/{mode}: {native} {logged}")
    props=(case/"constant/transportProperties").read_text()
    for entry in ("activityMode htHto;","auditActivityFlux true;",
                  "kOx [0 0 -1 0 0 0 0] 0.4;",
                  "lambdaT [0 0 -1 0 0 0 0] 0.2;",
                  "DHT [0 2 -1 0 0 0 0] 0.01;",
                  "DHTO [0 2 -1 0 0 0 0] 0.02;",
                  "momentumMode pisoLaminar;"):
        if entry not in props: fail(f"wrong physical coefficient {entry}")
    if "fluxRequired { default no; p; HT; HTO; }" not in (case/"system/fvSchemes").read_text():
        fail("missing real species FVM matrix fluxRequired")
    pressure=(case/"0"/"p").read_text()
    for patch,val in (("xMin",-.1 if mode=="reverse" else 0 if mode=="zero" else .1),("xMax",0)):
        m=re.search(rf"\b{patch}\s*\{{\s*type\s+fixedValue;\s*value\s+uniform\s*("+NUM+r")\s*;",pressure)
        if not m or abs(float(m.group(1))-val)>1e-12:
            fail(f"wrong fixed pressure BC {mode}/{patch}")
    prev_h=prev_z=0.
    maxres=0.
    for i,((h,z),(hmin,hmax,hnet,hr),(zmin,zmax,znet,zr,rt),(qin,qout),(divphi,netvol,ke)) in enumerate(zip(rows,ht,hto,flow,div)):
        rh=h-prev_h+.02*hnet+.02*.6*h
        rz=z-prev_z+.02*znet+.02*.2*z-.02*.4*h
        rtot=h+z-prev_h-prev_z+.02*(hnet+znet)+.02*.2*(h+z)
        maxres=max(maxres,*map(abs,(rh,rz,rtot,hr,zr,rt,rh-hr,rz-zr,rtot-rt)))
        if max(abs(hnet-hmin-hmax),abs(znet-zmin-zmax))>1e-3:
            fail(f"nonzero insulating wall activity flux {mode} step={i}")
        if max(map(abs,(rh,rz,rtot,hr,zr,rt)))>1e-4:
            fail(f"species/total signed discrete Bq budget {mode} step={i} residual={rh} {rz} {rtot}")
        if max(abs(qin+qout),abs(netvol),abs(divphi))>1e-7:
            fail(f"volume flux or flow divergence nonconservative {mode} step={i}")
        if min(h,z) < -1e-8:
            fail(f"negative species inventory {mode} step={i}")
        prev_h,prev_z=h,z
    qin,qout=flow[-1]
    hl,hr,hn,_=ht[-1]
    zl,zr,zn,_,_=hto[-1]
    if mode=="zero":
        if (max(abs(v) for arr in (rows,ht,hto,flow) for row in arr for v in row)>1e-8
            or max(abs(c) for name in ("HT","HTO") for c in final[name])>1e-8
            or max(abs(x) for v in field(case/"0.4/U",kind="vector") for x in v)>1e-8):
            fail("zero-gradient control generated activity without flow")
    else:
        inlet_h,inlet_z=(hr,zr) if mode=="reverse" else (hl,zl)
        if not (inlet_h < -1000 and inlet_z < -1000 and min(rows[-1])>0):
            fail(f"wrong nonzero active species inlet flux {mode}: {inlet_h} {inlet_z}")
        if mode=="reverse" and not (qin>.003 and qout<-.003):
            fail("reversed pressure/volume flow sign")
        if mode=="warped" and not (qin<-.003 and qout>.003):
            fail("forward pressure/volume flow sign")
    print(f"M1C-3 mode={mode} cells=960 HT_final_Bq={rows[-1][0]:.12g} HTO_final_Bq={rows[-1][1]:.12g} "
          f"HT_in_Bq_s={(hr if mode=='reverse' else hl):.12g} HTO_in_Bq_s={(zr if mode=='reverse' else zl):.12g} "
          f"max_budget_Bq={maxres:.12g} final_qin_m3_s={qin:.12g} final_qout_m3_s={qout:.12g}")

if len(sys.argv)!=8:
    raise SystemExit("usage: check_open_activity.py CHECK_LOG FWD_CASE FWD_LOG REV_CASE REV_LOG ZERO_CASE ZERO_LOG")
quality=Path(sys.argv[1]).read_text()
if "Mesh OK." not in quality or not re.search(r"(?m)^\s*cells:\s*960\s*$",quality):
    fail("native checkMesh 960 cells")
for case in (sys.argv[2],sys.argv[4],sys.argv[6]):
    if audit_physical_points(case,True)!=665:
        fail("exact native mesh warped vertex reconstruction failed")
    pts=(Path(case)/"constant/polyMesh/points").read_bytes()
    if case!=sys.argv[2] and pts!=(Path(sys.argv[2])/"constant/polyMesh/points").read_bytes():
        fail("not exactly the same native warped geometry across three cases")
for case,log,mode in ((sys.argv[2],sys.argv[3],"warped"),
                      (sys.argv[4],sys.argv[5],"reverse"),
                      (sys.argv[6],sys.argv[7],"zero")):
    audit(case,log,mode)
print("PASS: M1C-3 official OpenFOAM v2606 warped pressure-open HT/HTO signed Bq conservation")
