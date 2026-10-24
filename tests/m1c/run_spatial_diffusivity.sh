#!/usr/bin/env bash
set -eo pipefail
set +e; source /usr/lib/openfoam/openfoam2606/etc/bashrc; set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 11
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c22-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then echo "M1C22_FAILED rc=$rc"; for f in "$work"/build.log "$work"/*/solve.log; do if [[ -f "$f" ]];then echo "=== $f ===";tail -35 "$f";fi;done;fi;rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere" && wclean && wmake > "$work/build.log" 2>&1) || { tail -65 "$work/build.log";exit 20; }
exe="$FOAM_USER_APPBIN/cteasAtmosphereFoam"
test -x "$exe"
for n in 8 16;do
 base="$work/base$n"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$base" zero "$n"x"$n"x"$n"
 python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$base" advect
 python3 /workspace/tests/m1c/make_gaussian_green.py "$base" normal
 python3 - "$base" <<'PYCONFIG'
from pathlib import Path
import sys
c=Path(sys.argv[1])
p=c/"constant/transportProperties"
p.write_text(p.read_text()+"activityConvectionBlend faceMules;\nactivityDiffusionCorrection faceMules;\nenforceNonnegativeActivity true;\nauditActivityExtrema true;\n")
f=c/"system/fvSchemes";s=f.read_text()
old="div(phi,HTO) Gauss upwind;";assert s.count(old)==1
s=s.replace(old,old+" div(phi,HTHigh) Gauss limitedLinear 1; div(phi,HTOHigh) Gauss limitedLinear 1;")
old="laplacianSchemes { default Gauss linear corrected; }";assert s.count(old)==1
s=s.replace(old,"laplacianSchemes { default Gauss linear corrected; laplacian(DHT,HTLowDiff) Gauss linear uncorrected; laplacian(DHTO,HTOLowDiff) Gauss linear uncorrected; }")
old="snGradSchemes { default corrected; }";assert s.count(old)==1
s=s.replace(old,"snGradSchemes { default corrected; snGrad(HTHighDiff) corrected; snGrad(HTLowDiff) uncorrected; snGrad(HTOHighDiff) corrected; snGrad(HTOLowDiff) uncorrected; }")
f.write_text(s)
PYCONFIG
 blockMesh -case "$base" > "$base/blockMesh.log" 2>&1 || exit 21
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$base" zero "$n"x"$n"x"$n"
 checkMesh -case "$base" -allGeometry -allTopology > "$base/meshcheck.log" 2>&1 || exit 22
 grep -q 'Mesh OK.' "$base/meshcheck.log" || exit 23
 postProcess -case "$base" -func writeCellCentres -time 0 > "$base/centres.log" 2>&1 || exit 24
 postProcess -case "$base" -func writeCellVolumes -time 0 > "$base/volumes.log" 2>&1 || exit 25
 for mode in base uniform variable; do
   c="$work/$mode$n"
   if [[ "$mode" != base ]];then
      cp -a "$base" "$c"
      python3 /workspace/tests/m1c/make_spatial_diffusivity.py "$c" "$((n*n*n))" "$([[ "$mode" == uniform ]] && echo constant || echo variable)"
      printf '%s\n' 'activityDiffusivityMode field;' >> "$c/constant/transportProperties"
   fi
   timeout 240s "$exe" -case "$c" > "$c/solve.log" 2>&1 || { echo "M1C22_NATIVE_FAILED n=$n mode=$mode";tail -55 "$c/solve.log";exit 30; }
   python3 /workspace/tests/m1c/check_spatial_diffusivity.py "$c" "$n" "$mode"
 done
 python3 /workspace/tests/m1c/check_spatial_diffusivity.py --compare "$work/base$n" "$work/uniform$n" "$work/variable$n" "$n"
done

# True native parallel processing for one spatially variable coefficient case.
c="$work/variable8mpi"
cp -a "$work/variable8" "$c"
rm -rf "$c/0.4"
cat > "$c/system/decomposeParDict" <<'DICT'
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
distributed no;
roots ();
DICT
decomposePar -case "$c" -time 0 > "$work/decompose.log" 2>&1 || { tail -70 "$work/decompose.log";exit 41; }
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 240s mpirun --oversubscribe -np 2 "$exe" -parallel -case "$c" > "$c/solve.log" 2>&1 || { tail -60 "$c/solve.log";exit 42; }
reconstructPar -case "$c" -latestTime > "$work/reconstruct.log" 2>&1 || { tail -60 "$work/reconstruct.log";exit 43; }
python3 /workspace/tests/m1c/check_robustness.py mpi "$work/variable8" "$c" "$c/solve.log" "2:x"
python3 /workspace/tests/m1c/check_spatial_diffusivity.py "$c" 8 variable
echo 'PASS: M1C22 actual two-rank variable coefficient processor-face Bq/HT/HTO parity'

# Native negative input controls. The physically invalid field and missing
# field must fail BEFORE a time step, not fall back to scalar D.
for invalid in missing negative invalidMode incompatibleMode; do
 c="$work/invalid-$invalid"
 cp -a "$work/variable8" "$c"
 rm -rf "$c/0.4"
 python3 - "$c" "$invalid" <<'PYBAD'
from pathlib import Path
import re,sys
c=Path(sys.argv[1]);mode=sys.argv[2]
if mode=="missing":(c/"0/DHTSpatial").unlink()
if mode=="negative":
 p=c/"0/DHTSpatial";s=p.read_text()
 s,n=re.subn(r"(internalField nonuniform List<scalar>\s+512\s*\(\s*)([-+0-9.eE]+)",r"\g<1>-0.01",s,count=1)
 assert n==1
 p.write_text(s)
if mode in ("invalidMode","incompatibleMode"):
 p=c/"constant/transportProperties";s=p.read_text()
 if mode=="invalidMode":
  assert s.count("activityDiffusivityMode field;")==1
  s=s.replace("activityDiffusivityMode field;","activityDiffusivityMode bogus;")
 else:
  assert s.count("activityDiffusionCorrection faceMules;")==1
  s=s.replace("activityDiffusionCorrection faceMules;","activityDiffusionCorrection none;")
 p.write_text(s)
PYBAD
 case "$invalid" in
 missing) marker='DHTSpatial';;
 negative) marker='M1C22_NONPOSITIVE_DIFFUSIVITY';;
 invalidMode) marker='M1C22_UNSUPPORTED_DIFFUSIVITY_MODE';;
 incompatibleMode) marker='M1C22_FIELD_REQUIRES_BOUNDED_HT_HTO_MULES';;
 esac
 if timeout 45s "$exe" -case "$c" > "$work/invalid-$invalid.log" 2>&1;then echo "FAIL: invalid spatial D $invalid accepted";exit 44;fi
 grep -q "$marker" "$work/invalid-$invalid.log" || { tail -40 "$work/invalid-$invalid.log";exit 45; }
 if grep -q '^Time = ' "$work/invalid-$invalid.log";then echo 'FAIL: invalid spatial input advanced clock';exit 46;fi
 echo "PASS: M1C22 native fail-closed $invalid"
done

echo 'PASS: M1C22 actual native spatial coefficient input / Bq conservation / uniform equivalence scoped periodic mesh'
