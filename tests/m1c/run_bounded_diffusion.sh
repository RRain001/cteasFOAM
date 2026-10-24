#!/usr/bin/env bash
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 19
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c21-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then echo "M1C21_FAILED rc=$rc"; for f in "$work"/build.log "$work"/*/solve.log; do if [[ -f "$f" ]];then echo ==== "$f";tail -45 "$f";fi;done;fi; rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere" && wclean && wmake > "$work/build.log" 2>&1) || { tail -70 "$work/build.log";exit 20; }
exe="$FOAM_USER_APPBIN/cteasAtmosphereFoam"
test -x "$exe" || exit 21
for setting in "fine:8:0.0025:160" "normal:8:0.005:80" "coarse:8:0.01:40" "spatial:16:0.005:80"; do
 IFS=: read -r name n dt steps <<< "$setting"
 c="$work/$name"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$c" zero "$n"x"$n"x"$n"
 python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$c" advect
 python3 /workspace/tests/m1c/make_gaussian_green.py "$c" normal
 python3 - "$c" "$dt" "$steps" <<'PYCONFIG'
from pathlib import Path
import sys
c=Path(sys.argv[1]);dt=sys.argv[2];steps=sys.argv[3]
p=c/"constant/transportProperties"
p.write_text(p.read_text()+"activityConvectionBlend faceMules;\nactivityDiffusionCorrection faceMules;\n"+"enforceNonnegativeActivity true;\nauditActivityExtrema true;\n")
f=c/"system/fvSchemes";s=f.read_text()
old="div(phi,HTO) Gauss upwind;"
assert s.count(old)==1
s=s.replace(old,old+" div(phi,HTHigh) Gauss limitedLinear 1; div(phi,HTOHigh) Gauss limitedLinear 1;")
old="laplacianSchemes { default Gauss linear corrected; }"
assert s.count(old)==1
s=s.replace(old, "laplacianSchemes { default Gauss linear corrected; laplacian(DHT,HTLowDiff) Gauss linear uncorrected; laplacian(DHTO,HTOLowDiff) Gauss linear uncorrected; }")
old="snGradSchemes { default corrected; }"
assert s.count(old)==1
s=s.replace(old,"snGradSchemes { default corrected; snGrad(HTHighDiff) corrected; snGrad(HTLowDiff) uncorrected; snGrad(HTOHighDiff) corrected; snGrad(HTOLowDiff) uncorrected; }")
f.write_text(s)
control=c/"system/controlDict";s=control.read_text()
for old,new in (("deltaT 0.005;",f"deltaT {dt};"),("writeInterval 80;",f"writeInterval {steps};")):
 assert s.count(old)==1,(old,s)
 s=s.replace(old,new)
control.write_text(s)
PYCONFIG
 blockMesh -case "$c" > "$c/blockMesh.log" 2>&1 || exit 22
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$c" zero "$n"x"$n"x"$n"
 checkMesh -case "$c" -allGeometry -allTopology > "$c/mesh-check.log" 2>&1 || exit 23
 grep -q "Mesh OK." "$c/mesh-check.log" || exit 24
 postProcess -case "$c" -func writeCellCentres -time 0 > "$c/centres.log" 2>&1 || exit 25
 postProcess -case "$c" -func writeCellVolumes -time 0 > "$c/volumes.log" 2>&1 || exit 26
 timeout 360s "$exe" -case "$c" > "$c/solve.log" 2>&1 || { echo "M1C21_NATIVE_FAILED $name";tail -55 "$c/solve.log";exit 30; }
 python3 /workspace/tests/m1c/check_bounded_diffusion.py "$c" "$n" "$steps" "$dt"
done
source /workspace/tests/m1c/verify_bounded_diffusion_mpi.sh
source /workspace/tests/m1c/verify_bounded_diffusion_open.sh
source /workspace/tests/m1c/verify_bounded_diffusion_failclosed.sh
echo "PASS: M1C21 native warped periodic fine-clock, true 2/4 rank and physical two-way open correction tests; independent science NOT_RUN"
