#!/usr/bin/env bash
# Real official v2606 geometry/ABC sensitivity at two interior warp amplitudes.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
rc=$?
set -e
echo "M1B-1f4 setup_rc=$rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for cmd in wmake wclean blockMesh checkMesh postProcess python3 timeout; do
  command -v "$cmd" >/dev/null || { echo "FAIL: missing $cmd"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-warpsensitivity-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "M1B-1f4 failed rc=$rc";
 for f in "$work"/compile.log "$work"/mesh*.log "$work"/warp*.log "$work"/check*.log "$work"/centres*.log "$work"/seed*.log "$work"/solve*.log "$work"/metrics*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -130 "$f"; fi
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
 echo "FAIL: hidden wmake compiler error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22

for tag in 04 08; do
  if [[ "$tag" == 04 ]]; then amplitude=0.04; else amplitude=0.08; fi
  for n in 8 16; do
    key="$tag"_"$n"
    case_dir="$work/case$key"
    python3 /workspace/tests/m1b1/make_abc_case.py "$case_dir" "$n"
    blockMesh -case "$case_dir" > "$work/mesh$key.log" 2>&1 || {
      cat "$work/mesh$key.log"; exit 23;
    }
    echo "[M1B-1f4] actual polyMesh warp n=$n amplitude=$amplitude"
    python3 /workspace/tests/m1b1/warp_abc_mesh.py warp "$case_dir" "$n" "$amplitude" \
      > "$work/warp$key.log" 2>&1 || { cat "$work/warp$key.log"; exit 24; }
    checkMesh -case "$case_dir" -allGeometry -allTopology \
      > "$work/check$key.log" 2>&1 || { tail -120 "$work/check$key.log"; exit 25; }
    grep -q "Mesh OK." "$work/check$key.log" || exit 25
    postProcess -case "$case_dir" -func writeCellCentres -time 0 \
      > "$work/centres$key.log" 2>&1 || { tail -120 "$work/centres$key.log"; exit 26; }
    test -s "$case_dir/0/C" || exit 26
    python3 /workspace/tests/m1b1/warp_abc_mesh.py seed "$case_dir" "$n" \
      > "$work/seed$key.log" 2>&1 || { cat "$work/seed$key.log"; exit 27; }
    sed -i '/^nPressureCorrectors 3;$/a nNonOrthogonalCorrectors 2;' \
      "$case_dir/constant/transportProperties"
    grep -q '^nNonOrthogonalCorrectors 2;$' \
      "$case_dir/constant/transportProperties" || exit 27
    timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" \
      > "$work/solve$key.log" 2>&1 || { tail -150 "$work/solve$key.log"; exit 28; }
    steps=$(grep -c '^CTEAS_M1B1 nonorthogonal_correctors=2$' "$work/solve$key.log" || true)
    solves=$(grep -c 'Solving for p,' "$work/solve$key.log" || true)
    echo "M1B-1f4 amplitude=$amplitude n=$n actual_steps=$steps p_solves=$solves"
    [[ "$steps" -eq 20 && "$solves" -eq 180 ]] || {
      echo "FAIL: actual inner nonorth pressure solve count"; exit 29;
    }
  done
  python3 /workspace/tests/m1b1/check_warped_abc.py \
    "$work/case$tag"_8 "$work/check$tag"_8.log "$work/solve$tag"_8.log \
    "$work/case$tag"_16 "$work/check$tag"_16.log "$work/solve$tag"_16.log \
    > "$work/metrics$tag.log" || {
      cat "$work/metrics$tag.log"; exit 30;
    }
  cat "$work/metrics$tag.log"
done
python3 /workspace/tests/m1b1/check_warp_amplitude_trend.py \
 "$work/metrics04.log" "$work/metrics08.log"
echo "PASS: true warped 3D mesh-amplitude pressure-flow scientific sensitivity"
