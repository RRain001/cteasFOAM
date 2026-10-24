#!/usr/bin/env bash
# M1C-20 sourced by frozen M1C-19 runner after baseline serial/MPI/open checks.
# Requires work/exe and already solved serial 8_face, open_face.
set -euo pipefail
for shape in "y:2:1_2_1" "z:2:1_1_2" "xy:4:2_2_1"; do
  IFS=: read -r direction nRanks dims <<< "$shape"
  IFS=_ read -r nx ny nz <<< "$dims"
  caseDir="$work/m1c20_mpi_$direction"
  cp -a "$work/8_face" "$caseDir"
  cat > "$caseDir/system/decomposeParDict" <<DICT
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $nRanks;
method simple;
simpleCoeffs { n ($nx $ny $nz); delta 0.001; }
distributed no;
roots ();
DICT
  decomposePar -case "$caseDir" -time 0 > "$work/decomp_$direction.log" 2>&1 || { tail -65 "$work/decomp_$direction.log";exit 81; }
  for ((rank=0;rank<nRanks;rank++)); do
    test -f "$caseDir/processor$rank/constant/polyMesh/owner" || exit 82
  done
  OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 240s \
    mpirun --oversubscribe -np "$nRanks" "$exe" -parallel -case "$caseDir" \
    > "$work/m1c20-mpi-$direction.log" 2>&1 || { tail -100 "$work/m1c20-mpi-$direction.log";exit 83; }
  reconstructPar -case "$caseDir" -latestTime > "$work/reconstruct_$direction.log" 2>&1 || { tail -65 "$work/reconstruct_$direction.log";exit 84; }
  python3 /workspace/tests/m1c/check_robustness.py mpi \
    "$work/8_face" "$caseDir" "$work/m1c20-mpi-$direction.log" "$nRanks:$direction"
done

# Same physical source/geometry/reference, exactly two distinct CFL-level dt.
for setting in "fine:0.0025:160" "coarse:0.01:40"; do
  IFS=: read -r name dt steps <<< "$setting"
  c="$work/m1c20_time_$name"
  cp -a "$work/8_face" "$c"
  python3 - "$c/system/controlDict" "$dt" "$steps" <<'PYCLOCK'
from pathlib import Path
import sys
p=Path(sys.argv[1]);s=p.read_text()
assert s.count("deltaT 0.005;")==1
assert s.count("writeInterval 80;")==1
p.write_text(s.replace("deltaT 0.005;",f"deltaT {sys.argv[2]};").replace("writeInterval 80;",f"writeInterval {sys.argv[3]};"))
PYCLOCK
  if [[ "$name" == fine ]]; then
    # Preserve the genuine M1C-19 corrected-diffusion counterexample as
    # a NEGATIVE TEST. The solver must terminate instead of clipping HT.
    if timeout 120s "$exe" -case "$c" > "$work/m1c20-corrected-diffusion-reject.log" 2>&1; then
      echo "FAIL: M1C20 native corrected nonorthogonal diffusion unexpectedly accepted fine dt"; exit 87
    fi
    grep -q "CTEAS_M1C19_INVALID_LOW species=HT" "$work/m1c20-corrected-diffusion-reject.log" || {
      tail -80 "$work/m1c20-corrected-diffusion-reject.log"; exit 88;
    }
    echo "M1C20_EXPECTED_NEGATIVE corrected_nonorthogonal_diffusion dt=0.0025 true_low_HT_negative_no_clipping"
  fi
  # Independent remedial numerical experiment: no claim that the original
  # corrected non-orthogonal diffusion passes this smaller native time step.
  # Replace *only in these two exploratory cases* the deferred corrected
  # diffusion by uncorrected two-point implicit diffusion (positive stencil).
  # Original M1C-19 full and default correction remain frozen/untouched.
  python3 - "$c/system/fvSchemes" <<'PYDIFF'
from pathlib import Path
import sys
p=Path(sys.argv[1]);s=p.read_text()
old="laplacianSchemes { default Gauss linear corrected; }"
assert s.count(old)==1
p.write_text(s.replace(old,"laplacianSchemes { default Gauss linear uncorrected; }"))
PYDIFF
  timeout 240s "$exe" -case "$c" > "$c/solve.log" 2>&1 || { echo "M1C20_TIME_UNCORRECTED_NATIVE_FAILED case=$name"; tail -90 "$c/solve.log"; exit 85; }
  python3 /workspace/tests/m1c/check_robustness.py time "$c" "$steps" "$dt"
done

# Truly reversed inletOutlet flow, at xMax; original xMin positive-flow experiment remains unchanged.
c="$work/m1c20_open_reverse"
cp -a "$work/open_face" "$c"
python3 - "$c" <<'PYREVERSE'
from pathlib import Path
import sys
case=Path(sys.argv[1])
u=case/"0/U";s=u.read_text()
for fragment in ("internalField uniform (0.2 0 0);",
                 "xMin { type fixedValue; value uniform (0.2 0 0); }",
                 "xMax { type fixedValue; value uniform (0.2 0 0); }"):
    assert s.count(fragment)==1,fragment
u.write_text(s.replace("(0.2 0 0)","(-0.2 0 0)"))
for name,val in (("HT","1000000"),("HTO","200000")):
    p=case/"0"/name;s=p.read_text()
    left=f"xMin {{ type inletOutlet; inletValue uniform {val}; value uniform 0; }}"
    right="xMax { type inletOutlet; inletValue uniform 0; value uniform 0; }"
    assert s.count(left)==1 and s.count(right)==1
    s=s.replace(left,f"xMin {{ type inletOutlet; inletValue uniform 0; value uniform 0; }}")
    s=s.replace(right,f"xMax {{ type inletOutlet; inletValue uniform {val}; value uniform 0; }}")
    p.write_text(s)
print("M1C20_NATIVE_REVERSED_OPEN active_inlet=xMax U=(-0.2,0,0)")
PYREVERSE
timeout 180s "$exe" -case "$c" > "$c/solve.log" 2>&1 || { tail -85 "$c/solve.log";exit 86; }
python3 - "$c" <<'PYOPEN'
from pathlib import Path
import math,re,sys
sys.path.insert(0,"/workspace/tests/m1c")
from check_gaussian_green import field,DIM,INV,HTF,ZTF
case=Path(sys.argv[1]);s=(case/"solve.log").read_text()
steps=re.findall(r"(?m)^Time = ",s)
h=HTF.findall(s);z=ZTF.findall(s);inv=INV.findall(s)
if any(len(x)!=20 for x in (steps,h,z,inv)):
    raise SystemExit("FAIL: M1C20 reversed open true native 20-step flux/inventory evidence")
boundaries=[]
for j,(hh,zz) in enumerate(zip(h,z),1):
    left,right,net,res=map(float,hh)
    leftz,rightz,netz,resz,rest=map(float,zz)
    if max(abs(res),abs(resz),abs(rest))>1e-5:
        raise SystemExit(f"FAIL: M1C20 reversed native open Bq budget at {j}")
    if not all(math.isfinite(x) for x in (left,right,net,leftz,rightz,netz,res,resz,rest)):
        raise SystemExit("FAIL: nonfinite reversed physical flux")
    boundaries.append((left,right))
if not any(right<-1e4 and abs(right)>abs(left) for left,right in boundaries):
    raise SystemExit("FAIL: xMax expected real negative signed incoming flux absent")
for species in ("HT","HTO"):
    arr=field(case/"0.4"/species,960,"scalar",DIM)
    if not all(map(math.isfinite,arr)) or min(arr)<-1e-8:
        raise SystemExit("FAIL: reversed real OpenFOAM positive species fields")
if len(re.findall(r"CTEAS_M1C19 species=HT face_lambda_min=",s))!=20 or len(re.findall(r"CTEAS_M1C19 species=HTO face_lambda_min=",s))!=20:
    raise SystemExit("FAIL: reversed open per-step MULES limiter absent")
print(f"M1C20_REVERSED_OPEN max_HT_xMax_inflow_Bq_s={-min(x[1] for x in boundaries):.12g}")
print("PASS: M1C20 native reversed physical open signed HT/HTO/total Bq and positivity")
PYOPEN
echo "PASS: M1C20 diagnostic cases: multi-axis 2/4 rank, reversed open, and separate uncorrected-diffusion 40/160-step positive Bq checks; corrected-diffusion fine-dt NEGATIVE confirmed"
