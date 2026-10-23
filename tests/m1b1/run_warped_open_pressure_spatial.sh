#!/usr/bin/env bash
# M1B-1g-3: actual three-level 3D warped duct spatial refinement.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
rc=$?
set -e
echo "M1B-1g-3 official OpenFOAM setup_rc=$rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for cmd in python3 wmake wclean blockMesh checkMesh postProcess timeout; do
 command -v "$cmd" >/dev/null || { echo "FAIL: missing $cmd"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-duct-h-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1B-1g-3 spatial study rc=$rc";
 for f in "$work"/build.log "$work"/*/check.log "$work"/*/mesh.log "$work"/*/centre.log "$work"/*/solve.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -60 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || {
 cat "$work/build.log"; exit 21;
}
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then
 echo "FAIL: hidden compiler error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for dim in 8x20x6 12x30x9 16x40x12; do
 case_dir="$work/$dim"
 echo "M1B-1g-3 start resolution=$dim"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" warped "$dim"
 blockMesh -case "$case_dir" > "$case_dir/mesh.log" 2>&1 || {
  tail -90 "$case_dir/mesh.log"; exit 23;
 }
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" warped "$dim" || exit 24
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/check.log" 2>&1 || {
  tail -100 "$case_dir/check.log"; exit 25;
 }
 grep -q 'Mesh OK.' "$case_dir/check.log" || exit 25
 postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$case_dir/centre.log" 2>&1 || {
  tail -100 "$case_dir/centre.log"; exit 26;
 }
 test -s "$case_dir/0/C" || exit 26
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$case_dir" warped "$dim" || exit 27
 python3 - "$case_dir/system/controlDict" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1])
s=p.read_text()
if s.count("deltaT 0.02;")!=1:
    raise SystemExit("FAIL: missing base timestep")
s=s.replace("deltaT 0.02;", "deltaT 0.01;")
p.write_text(s)
PY
 timeout --signal=TERM --kill-after=10s 600s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || {
  echo "FAIL: spatial solve exited nonzero at $dim"
  tail -100 "$case_dir/solve.log"; exit 31;
 }
 echo "M1B-1g-3 completed resolution=$dim"
done
python3 /workspace/tests/m1b1/check_warped_open_pressure_spatial.py \
 "$work/8x20x6" "$work/8x20x6/check.log" "$work/8x20x6/solve.log" \
 "$work/12x30x9" "$work/12x30x9/check.log" "$work/12x30x9/solve.log" \
 "$work/16x40x12" "$work/16x40x12/check.log" "$work/16x40x12/solve.log"
echo "PASS: M1B-1g-3 official-v2606 3D physical warped duct spatial convergence"
