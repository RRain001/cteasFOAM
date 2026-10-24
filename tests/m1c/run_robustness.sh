#!/usr/bin/env bash
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c19-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then echo "FAIL: M1C19 rc=$rc"; for f in "$work"/build.log "$work"/*/solve.log; do if [[ -f "$f" ]]; then echo "== $f =="; tail -65 "$f";fi;done;fi; rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere" && wclean && wmake > "$work/build.log" 2>&1) || { tail -95 "$work/build.log"; exit 20; }
exe="$FOAM_USER_APPBIN/cteasAtmosphereFoam"
test -x "$exe" || exit 21
for n in 8 16 24; do
  for mode in base face; do
    c="$work/${n}_$mode"
    python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$c" zero "${n}x${n}x${n}"
    python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$c" advect
    python3 /workspace/tests/m1c/make_gaussian_green.py "$c" normal
    python3 - "$c" "$mode" <<'PYCASE'
from pathlib import Path
import sys
case=Path(sys.argv[1]);mode=sys.argv[2]
p=case/'constant/transportProperties'
p.write_text(p.read_text()+'auditActivityExtrema true;\nenforceNonnegativeActivity true;\n'+('activityConvectionBlend faceMules;\n' if mode=='face' else ''))
if mode=='face':
 f=case/'system/fvSchemes';txt=f.read_text();needle='div(phi,HTO) Gauss upwind;'
 assert txt.count(needle)==1
 txt=txt.replace(needle,needle+' div(phi,HTHigh) Gauss limitedLinear 1; div(phi,HTOHigh) Gauss limitedLinear 1;')
 f.write_text(txt)
print(f'M1C19_REAL_NATIVE_CASE n={sys.argv[1]} mode={mode}')
PYCASE
    blockMesh -case "$c" > "$c/blockMesh.log" 2>&1 || { tail -40 "$c/blockMesh.log"; exit 22; }
    python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$c" zero "${n}x${n}x${n}"
    checkMesh -case "$c" -allGeometry -allTopology > "$c/mesh-check.log" 2>&1 || { tail -40 "$c/mesh-check.log"; exit 23; }
    grep -q 'Mesh OK.' "$c/mesh-check.log" || exit 24
    postProcess -case "$c" -func writeCellCentres -time 0 > "$c/centres.log" 2>&1 || exit 25
    postProcess -case "$c" -func writeCellVolumes -time 0 > "$c/volumes.log" 2>&1 || exit 26
    timeout 600s "$exe" -case "$c" > "$c/solve.log" 2>&1 || { tail -80 "$c/solve.log"; exit 31; }
  done
  python3 /workspace/tests/m1c/check_local_fct.py --compare "$work" "$n"
done
source /workspace/tests/m1c/verify_local_fct_mpi.sh
source /workspace/tests/m1c/verify_local_fct_open.sh
source /workspace/tests/m1c/verify_robustness.sh
source /workspace/tests/m1c/verify_local_fct_negative.sh
# Native negative control: invalid requested high-order scheme must fail.
case_dir="$work/8_face"
sed -i 's/div(phi,HTHigh) Gauss limitedLinear 1;/div(phi,HTHigh) Gauss unrecognizedHighScheme19;/' "$case_dir/system/fvSchemes"
if timeout 60s "$exe" -case "$case_dir" > "$work/invalid-scheme.log" 2>&1; then echo 'FAIL: M1C19 invalid high-order scheme accepted';exit 40;fi
grep -q unrecognizedHighScheme19 "$work/invalid-scheme.log" || exit 41
echo 'PASS: M1C-19 bad native high-order scheme rejected'
echo 'PASS: M1C-19 experimental facewise MULES true native 3-D periodic, 2-rank MPI, open-boundary positive Bq flux tests; independent external science NOT_RUN'
