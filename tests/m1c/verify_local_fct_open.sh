#!/usr/bin/env bash
# Sourced by run_local_fct.sh, matching its compiled source.
for mode in base face; do
  c="$work/open_$mode"
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$c" zero
  blockMesh -case "$c" > "$c/blockMesh.log" 2>&1 || exit 61
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$c" zero
  checkMesh -case "$c" -allGeometry -allTopology > "$c/mesh-check.log" 2>&1 || exit 62
  grep -q 'Mesh OK.' "$c/mesh-check.log" || exit 63
  python3 /workspace/tests/m1c/make_open_activity.py "$c" zero
  python3 - "$c" "$mode" <<'PYCASE'
from pathlib import Path
import sys
case=Path(sys.argv[1]);mode=sys.argv[2]
p=case/'constant/transportProperties'
src=p.read_text();assert src.count('momentumMode pisoLaminar;')==1
p.write_text(src.replace('momentumMode pisoLaminar;','momentumMode prescribed;')+
             'enforceNonnegativeActivity true;\nauditActivityExtrema true;\n'+
             ('activityConvectionBlend faceMules;\n' if mode=='face' else ''))
u=case/'0/U';s=u.read_text()
before='internalField uniform (0 0 0);'
assert s.count(before)==1
s=s.replace(before,'internalField uniform (0.2 0 0);')
for patch in ('xMin','xMax'):
    before=patch+' { type pressureInletOutletVelocity; value uniform (0 0 0); }'
    after=patch+' { type fixedValue; value uniform (0.2 0 0); }'
    assert s.count(before)==1
    s=s.replace(before,after)
u.write_text(s)
if mode=='face':
    f=case/'system/fvSchemes';s=f.read_text();needle='div(phi,HTO) Gauss upwind;'
    assert s.count(needle)==1
    f.write_text(s.replace(needle,needle+
       ' div(phi,HTHigh) Gauss limitedLinear 1; div(phi,HTOHigh) Gauss limitedLinear 1;'))
print('M1C19_OPEN_NATIVE_INPUT mode='+mode+' mesh=960 cells xMin/xMax actual nonperiodic')
PYCASE
  timeout 180s "$exe" -case "$c" > "$c/solve.log" 2>&1 || { tail -70 "$c/solve.log"; exit 64; }
done
python3 - "$work/open_base" "$work/open_face" <<'PYVERIFY'
import re,sys,math
from pathlib import Path
sys.path.insert(0,'/workspace/tests/m1c')
from check_gaussian_green import field,DIM,INV,HTF,ZTF
a,b=map(Path,sys.argv[1:])
num=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
for name,c in (('base',a),('face',b)):
 s=(c/'solve.log').read_text()
 steps=re.findall(r'^Time = ',s,re.M)
 budgets=HTF.findall(s);other=ZTF.findall(s);inv=INV.findall(s)
 if len(steps)!=20 or len(budgets)!=20 or len(other)!=20 or len(inv)!=20:
  raise SystemExit('FAIL: M1C19 open genuine 20-step native HT/HTO + flux budget incomplete '+name)
 for row in budgets:
  if abs(float(row[3]))>1e-5:raise SystemExit('FAIL: M1C19 open HT signed budget '+name)
 for row in other:
  if max(abs(float(x)) for x in row[3:])>1e-5:raise SystemExit('FAIL: M1C19 open HTO/total signed budget '+name)
 boundary=[tuple(map(float,v)) for v in re.findall(r'flux_HT_xMin_Bq_s=('+num+r') flux_HT_xMax_Bq_s=('+num+r')',s)]
 if len(boundary)!=20 or not any(x < -1e4 for x,y in boundary):
  raise SystemExit('FAIL: M1C19 missing real nonzero inflow HT boundary flux '+name)
 if name=='face':
  if len(re.findall(r'CTEAS_M1C19 species=HT face_lambda_min=',s))!=20 or len(re.findall(r'CTEAS_M1C19 species=HTO face_lambda_min=',s))!=20:
   raise SystemExit('FAIL: M1C19 missing actual open-domain face-limiter log')
 for sp in ('HT','HTO'):
  v=field(c/'0.4'/sp,960,'scalar',DIM)
  if not all(map(math.isfinite,v)) or min(v)<-1e-8:
   raise SystemExit('FAIL: M1C19 true open finite-volume native negative/nonfinite field '+name+' '+sp)
 print(f'M1C19_OPEN_NATIVE mode={name} max_boundary_inflow_HT_Bq_s={-min(x for x,y in boundary):.12g}')
for sp in ('HT','HTO'):
 v1=field(a/'0.4'/sp,960,'scalar',DIM);v2=field(b/'0.4'/sp,960,'scalar',DIM)
 d=max(abs(x-y) for x,y in zip(v1,v2))
 print(f'M1C19_OPEN_FACE_DELTA {sp} max_abs_Bq_m3={d:.12g}')
print('PASS: M1C-19 genuine open inletOutlet wall/patch face-FCT solver, 20-step signed HT/HTO boundary Bq and nonnegative cells')
PYVERIFY
