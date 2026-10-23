#!/usr/bin/env bash
# M1B-1g-6: real 3-D pressure-driven finite-volume nonuniform c budget.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B-1g-6 setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for app in python3 blockMesh checkMesh postProcess wmake wclean timeout; do
 command -v "$app" >/dev/null || { echo "FAIL: missing actual OpenFOAM executable $app"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-open-tracer-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1B-1g-6 rc=$rc";
 for f in "$work"/build.log "$work"/mesh.log "$work"/check.log "$work"/centres.log "$work"/solve*.log "$work"/negative.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -100 "$f"; fi
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
 echo "FAIL: hidden compiler error"; exit 22;
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22

for mode in warped reverse zero; do
 case_dir="$work/$mode"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" "$mode"
 if [[ "$mode" == warped ]]; then
  blockMesh -case "$case_dir" > "$work/mesh.log" 2>&1 || { tail -90 "$work/mesh.log"; exit 23; }
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" "$mode"
  checkMesh -case "$case_dir" -allGeometry -allTopology > "$work/check.log" 2>&1 || {
   tail -90 "$work/check.log"; exit 24;
  }
  grep -q "Mesh OK." "$work/check.log" || exit 24
  postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$work/centres.log" 2>&1 || {
   tail -90 "$work/centres.log"; exit 25;
  }
 else
  cp -a "$work/warped/constant/polyMesh" "$case_dir/constant/"
  cp -a "$work/warped/0/C" "$case_dir/0/"
 fi
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$case_dir" "$mode"
 python3 /workspace/tests/m1b1/make_open_pressure_tracer.py "$case_dir" "$mode"
 timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$work/solve-$mode.log" 2>&1 || {
  tail -100 "$work/solve-$mode.log"; exit 31;
 }
 echo "M1B-1g-6 completed mode=$mode steps=$(grep -c '^CTEAS_M1B1 nonorthogonal_correctors=2$' "$work/solve-$mode.log") pressure_solves=$(grep -c 'Solving for p,' "$work/solve-$mode.log") tracer_flux_records=$(grep -c '^CTEAS_M1A scalar_flux_xMin=' "$work/solve-$mode.log")"
done
args=( "$work/check.log" "$work/warped" "$work/solve-warped.log"
       "$work/reverse" "$work/solve-reverse.log"
       "$work/zero" "$work/solve-zero.log" )
python3 /workspace/tests/m1b1/check_open_pressure_tracer.py "${args[@]}"
cp "$work/solve-warped.log" "$work/solve-corrupt.log"
python3 - "$work/solve-corrupt.log" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);data=p.read_text()
pat=r"scalar_budget_residual=([-+0-9.eE]+)"
changed,count=re.subn(pat,"scalar_budget_residual=0.01",data,count=1)
if count!=1: raise SystemExit("FAIL: unable to inject a real tracer-budget corruption")
p.write_text(changed)
PY
if python3 /workspace/tests/m1b1/check_open_pressure_tracer.py \
 "$work/check.log" "$work/warped" "$work/solve-corrupt.log" \
 "$work/reverse" "$work/solve-reverse.log" \
 "$work/zero" "$work/solve-zero.log" > "$work/negative.log" 2>&1; then
 echo "FAIL: deliberately corrupted scalar balance went undetected"; exit 35
fi
grep -q 'FAIL: scalar budget residual' "$work/negative.log" || {
 cat "$work/negative.log"; echo "FAIL: wrong negative control rejection"; exit 36;
}
echo "PASS: deliberate +0.01 m3 scalar-budget corruption rejected"
echo "PASS: M1B-1g-6 official OpenFOAM v2606 nonuniform tracer conservation and reverse-flow"
