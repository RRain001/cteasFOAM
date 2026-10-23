#!/usr/bin/env bash
# Real nonorthogonal 3D ABC geometry, centroids and nonlinear pressure flow.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B-1f2 setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for exe in wmake wclean blockMesh checkMesh postProcess python3 timeout; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-warped3d-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "M1B-1f2 nonorthogonal regression failed rc=$rc";
 for f in "$work"/compile.log "$work"/mesh*.log "$work"/warp*.log "$work"/check*.log "$work"/centres*.log "$work"/seed*.log "$work"/solve*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -140 "$f"; fi
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
 echo "FAIL: hidden C++ error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for n in 8 16; do
 case_dir="$work/case$n"
 python3 /workspace/tests/m1b1/make_abc_case.py "$case_dir" "$n"
 blockMesh -case "$case_dir" > "$work/mesh$n.log" 2>&1 || {
   cat "$work/mesh$n.log"; exit 23;
 }
 echo "[M1B-1f] warp actual OpenFOAM polyMesh n=$n amplitude=.06"
 python3 /workspace/tests/m1b1/warp_abc_mesh.py warp "$case_dir" "$n" 0.06 \
   > "$work/warp$n.log" 2>&1 || { cat "$work/warp$n.log"; exit 24; }
 cat "$work/warp$n.log"
 checkMesh -case "$case_dir" -allGeometry -allTopology \
   > "$work/check$n.log" 2>&1 || {
     cat "$work/check$n.log"; exit 25;
   }
 grep -q 'Mesh OK.' "$work/check$n.log" || {
   echo "FAIL: checkMesh not accepted"; tail -120 "$work/check$n.log"; exit 25;
 }
 echo "[M1B-1f] write actual cell centres via OpenFOAM"
 postProcess -case "$case_dir" -func writeCellCentres -time 0 \
   > "$work/centres$n.log" 2>&1 || {
     cat "$work/centres$n.log"; exit 26;
   }
 test -s "$case_dir/0/C" || {
   echo "FAIL: writeCellCentres did not emit 0/C"; exit 26;
 }
 python3 /workspace/tests/m1b1/warp_abc_mesh.py seed "$case_dir" "$n" \
   > "$work/seed$n.log" 2>&1 || { cat "$work/seed$n.log"; exit 27; }
 cat "$work/seed$n.log"
 sed -i '/^nPressureCorrectors 3;$/a nNonOrthogonalCorrectors 2;' \
   "$case_dir/constant/transportProperties"
 grep -q '^nNonOrthogonalCorrectors 2;$' \
   "$case_dir/constant/transportProperties" || exit 27
 echo "[M1B-1f2] real three-dimensional nonorthogonal inner PISO n=$n"
 timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" \
   > "$work/solve$n.log" 2>&1 || {
     cat "$work/solve$n.log"; exit 28;
   }
 # Verify actual pressure matrix executions, not merely config text.
 # 20 timesteps x 3 outer PISO x (2+1) inner nonorth = 180 solves.
 actual_20=$(grep -c '^CTEAS_M1B1 nonorthogonal_correctors=2$' \
   "$work/solve$n.log" || true)
 actual_180=$(grep -c 'Solving for p,' "$work/solve$n.log" || true)
 echo "M1B-1f2 n=$n actual_nonorth_records=$actual_20 actual_p_solutions=$actual_180"
 [[ "$actual_20" -eq 20 && "$actual_180" -eq 180 ]] || {
   echo "FAIL: did not execute 180 actual pressure solves"; exit 29;
 }
done
python3 /workspace/tests/m1b1/check_warped_abc.py \
 "$work/case8" "$work/check8.log" "$work/solve8.log" \
 "$work/case16" "$work/check16.log" "$work/solve16.log" \
 > "$work/nonorth_metrics.log"
cat "$work/nonorth_metrics.log"
python3 /workspace/tests/m1b1/check_nonorth_pressure_consistency.py \
 "$work/nonorth_metrics.log"
echo "PASS: actual inner nonorth pressure correctors and warped ABC scientific gates"
