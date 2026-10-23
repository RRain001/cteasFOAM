#!/usr/bin/env bash
# Actual unsteady, fully 3D, nonlinear periodic Beltrami pressure-PISO PDE.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B1 ABC setup_rc=$setup_rc OpenFOAM=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for exe in wmake wclean blockMesh python3 timeout; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-abc3d-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "M1B1 ABC real 3D regression FAILED rc=$rc";
 for file in "$work"/compile.log "$work"/mesh*.log "$work"/solve*.log; do
   if [[ -f "$file" ]]; then echo "=== $file ==="; tail -150 "$file"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(
 cd "$work/atmosphere"
 wclean
 wmake > "$work/compile.log" 2>&1 || { cat "$work/compile.log"; exit 21; }
)
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/compile.log"; then
 echo "FAIL: hidden ABC C++ compiler error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for n in 8 16; do
 python3 /workspace/tests/m1b1/make_abc_case.py "$work/case$n" "$n"
 blockMesh -case "$work/case$n" > "$work/mesh$n.log" 2>&1 || {
   cat "$work/mesh$n.log"; exit 23;
 }
 echo "[M1B-1d] solve real ABC fully 3D n=$n^3, 20 physical steps"
 timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" \
   -case "$work/case$n" > "$work/solve$n.log" 2>&1 || {
     cat "$work/solve$n.log"; exit 24;
   }
 grep -q 'CTEAS_M1B1 max_div_phi=' "$work/solve$n.log" || {
   echo "FAIL: no actual solved pressure-corrected 3D face-divergence logs"; exit 25;
 }
done
python3 /workspace/tests/m1b1/check_abc.py \
 "$work/case8" "$work/solve8.log" "$work/case16" "$work/solve16.log"
echo "PASS: real nonlinear three-component 3D ABC analytic Navier-Stokes regression"
