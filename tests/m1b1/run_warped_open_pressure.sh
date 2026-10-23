#!/usr/bin/env bash
# M1B-1g: actual 3D warped nonperiodic pressure-driven duct and zero-dp control.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B-1g OpenFOAM setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for exe in wmake wclean blockMesh checkMesh postProcess python3 timeout; do
  command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-open-warp-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "M1B-1g failed rc=$rc; retaining failure evidence below";
 for f in "$work"/build.log "$work"/mesh*.log "$work"/check*.log "$work"/centres*.log "$work"/solve*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -120 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(
 cd "$work/atmosphere"
 wclean
 wmake > "$work/build.log" 2>&1 || { cat "$work/build.log"; exit 21; }
)
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then
 echo "FAIL: hidden compiler error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for mode in flat warped zero; do
 case_dir="$work/$mode"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" "$mode"
 blockMesh -case "$case_dir" > "$work/mesh$mode.log" 2>&1 || {
  cat "$work/mesh$mode.log"; exit 23;
 }
 if [[ "$mode" != flat ]]; then
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" "$mode" || exit 24
 fi
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$work/check$mode.log" 2>&1 || {
  tail -120 "$work/check$mode.log"; exit 25;
 }
 grep -q "Mesh OK." "$work/check$mode.log" || { tail -120 "$work/check$mode.log"; exit 25; }
 postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$work/centres$mode.log" 2>&1 || {
  tail -120 "$work/centres$mode.log"; exit 26;
 }
 test -s "$case_dir/0/C" || { echo "FAIL: native centres absent"; exit 26; }
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$case_dir" "$mode" || exit 27
 timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$work/solve$mode.log" 2>&1 || {
  tail -130 "$work/solve$mode.log"; exit 28;
 }
 echo "M1B-1g mode=$mode actual_steps=$(grep -c '^CTEAS_M1B1 nonorthogonal_correctors=2$' "$work/solve$mode.log") pressure_solves=$(grep -c 'Solving for p,' "$work/solve$mode.log")"
done
python3 /workspace/tests/m1b1/check_warped_open_pressure.py \
 "$work/flat" "$work/checkflat.log" "$work/solveflat.log" \
 "$work/warped" "$work/checkwarped.log" "$work/solvewarped.log" \
 "$work/zero" "$work/checkzero.log" "$work/solvezero.log"
echo "PASS: M1B-1g 3D distorted/nonperiodic fixed pressure channel scientific software gate"
