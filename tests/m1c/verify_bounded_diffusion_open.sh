#!/usr/bin/env bash
for direction in forward reverse;do
  c="$work/open_$direction"
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$c" zero
  blockMesh -case "$c" > "$c/blockMesh.log" 2>&1 || exit 91
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$c" zero
  checkMesh -case "$c" -allGeometry -allTopology > "$c/mesh-check.log" 2>&1 || exit 92
  grep -q "Mesh OK." "$c/mesh-check.log" || exit 93
  python3 /workspace/tests/m1c/make_open_activity.py "$c" "$([[ "$direction" == reverse ]] && echo reverse || echo zero)"
  python3 - "$c" "$direction" <<'PYOPEN'
from pathlib import Path
import sys
case=Path(sys.argv[1]);direction=sys.argv[2]
p=case/"constant/transportProperties";s=p.read_text()
assert s.count("momentumMode pisoLaminar;")==1
s=s.replace("momentumMode pisoLaminar;","momentumMode prescribed;")
p.write_text(s+"activityConvectionBlend faceMules;\nactivityDiffusionCorrection faceMules;\nenforceNonnegativeActivity true;\nauditActivityExtrema true;\n")
u=case/"0/U";s=u.read_text()
old="internalField uniform (0 0 0);";assert s.count(old)==1
flow="-0.2" if direction=="reverse" else "0.2"
s=s.replace(old,f"internalField uniform ({flow} 0 0);")
for patch in ("xMin","xMax"):
    old=patch+" { type pressureInletOutletVelocity; value uniform (0 0 0); }"
    assert s.count(old)==1
    s=s.replace(old,patch+f" {{ type fixedValue; value uniform ({flow} 0 0); }}")
u.write_text(s)
schemes=case/"system/fvSchemes";s=schemes.read_text()
old="div(phi,HTO) Gauss upwind;";assert s.count(old)==1
s=s.replace(old,old+" div(phi,HTHigh) Gauss limitedLinear 1; div(phi,HTOHigh) Gauss limitedLinear 1;")
old="laplacianSchemes { default Gauss linear corrected; }";assert s.count(old)==1
s=s.replace(old,"laplacianSchemes { default Gauss linear corrected; laplacian(DHT,HTLowDiff) Gauss linear uncorrected; laplacian(DHTO,HTOLowDiff) Gauss linear uncorrected; }")
old="snGradSchemes { default corrected; }";assert s.count(old)==1
s=s.replace(old,"snGradSchemes { default corrected; snGrad(HTHighDiff) corrected; snGrad(HTLowDiff) uncorrected; snGrad(HTOHighDiff) corrected; snGrad(HTOLowDiff) uncorrected; }")
schemes.write_text(s)
print(f"M1C21_OPEN_SETUP direction={direction} U_x={flow} xMin/xMax real inletOutlet")
PYOPEN
  timeout 240s "$exe" -case "$c" > "$c/solve.log" 2>&1 || { echo "M1C21_OPEN_FAILED $direction";tail -70 "$c/solve.log";exit 94; }
  python3 - "$c" "$direction" <<'PYAUDIT'
import sys,math,re
from pathlib import Path
sys.path.insert(0,'/workspace/tests/m1c')
from check_gaussian_green import field,DIM,HTF,ZTF,INV,NUM
c=Path(sys.argv[1]);d=sys.argv[2];s=(c/"solve.log").read_text()
ht=HTF.findall(s);z=ZTF.findall(s);inventory=INV.findall(s)
if len(ht)!=20 or len(z)!=20 or len(inventory)!=20:raise SystemExit("FAIL: M1C21 open actual native signed HT/HTO flux chronology")
for j,(h,v) in enumerate(zip(ht,z),1):
    if max(abs(float(h[3])),abs(float(v[3])),abs(float(v[4])))>1e-5:
        raise SystemExit(f"FAIL: M1C21 open signed Bq balance at step {j}")
    left,right=map(float,h[:2])
    if d=="forward" and not (left<-1e4 and abs(left)>abs(right)):
        raise SystemExit(f"FAIL: M1C21 xMin incoming physical boundary at step {j}")
    if d=="reverse" and not (right<-1e4 and abs(right)>abs(left)):
        raise SystemExit(f"FAIL: M1C21 xMax incoming physical boundary at step {j}")
for sp in ("HT","HTO"):
    a=field(c/"0.4"/sp,960,'scalar',DIM)
    if not all(map(math.isfinite,a)) or min(a)<-1e-8:raise SystemExit("FAIL: open native species negative/nonfinite "+sp)
    correction=re.findall(r"CTEAS_M1C21_DIFFUSION species="+sp+r" max_abs_high_minus_low_Bq_s=("+NUM+r") nonzero_internal_faces=(\d+)",s)
    if len(correction)!=20 or not any(float(flux)>1e-12 and int(n)>0 for flux,n in correction):
        raise SystemExit("FAIL: open actual nonorthogonal correction absent "+sp)
print(f"PASS: M1C21 genuinely nonperiodic {d} inletOutlet positive species and 20-step signed Bq budgets with nonzero corrected diffusion")
PYAUDIT
done
