#!/usr/bin/env bash
# M1B-1g-7: genuine v2606 pressure-driven tracer D=0 vs D=0.01.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B-1g-7 setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for app in python3 blockMesh checkMesh postProcess wmake wclean timeout; do
 command -v "$app" >/dev/null || { echo "FAIL: missing OpenFOAM $app"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-open-diffusion-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1B-1g-7 rc=$rc";
 for f in "$work"/build.log "$work"/mesh.log "$work"/check.log "$work"/centres.log "$work"/solve*.log "$work"/negative.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -85 "$f"; fi
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

for mode in warped reverse zero; do
 for diffusion in zero positive; do
  if [[ "$diffusion" == positive ]]; then D=0.01; else D=0; fi
  case_dir="$work/$mode-$diffusion"
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" "$mode"
  if [[ "$mode" == warped && "$diffusion" == zero ]]; then
   blockMesh -case "$case_dir" > "$work/mesh.log" 2>&1 || { tail -90 "$work/mesh.log"; exit 23; }
   python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" "$mode"
   checkMesh -case "$case_dir" -allGeometry -allTopology > "$work/check.log" 2>&1 || {
    tail -90 "$work/check.log"; exit 24;
   }
   grep -q 'Mesh OK.' "$work/check.log" || exit 24
   postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$work/centres.log" 2>&1 || {
    tail -90 "$work/centres.log"; exit 25;
   }
  else
   cp -a "$work/warped-zero/constant/polyMesh" "$case_dir/constant/"
   cp -a "$work/warped-zero/0/C" "$case_dir/0/"
  fi
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$case_dir" "$mode"
  python3 /workspace/tests/m1b1/make_open_pressure_tracer.py "$case_dir" "$mode" "$D"
  timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$work/solve-$mode-$diffusion.log" 2>&1 || {
   tail -100 "$work/solve-$mode-$diffusion.log"; exit 31;
  }
  echo "M1B-1g-7 mode=$mode D=$D steps=$(grep -c '^CTEAS_M1B1 nonorthogonal_correctors=2$' "$work/solve-$mode-$diffusion.log") pressure_solves=$(grep -c 'Solving for p,' "$work/solve-$mode-$diffusion.log") scalar_flux_records=$(grep -c '^CTEAS_M1A scalar_flux_xMin=' "$work/solve-$mode-$diffusion.log")"
 done
done

for diffusion in zero positive; do
 args=( "$work/check.log" )
 for mode in warped reverse zero; do
  args+=( "$work/$mode-$diffusion" "$work/solve-$mode-$diffusion.log" )
 done
 python3 /workspace/tests/m1b1/check_open_pressure_tracer.py "${args[@]}"
done
args=()
for mode in warped reverse zero; do
 args+=( "$work/$mode-zero" "$work/solve-$mode-zero.log"
         "$work/$mode-positive" "$work/solve-$mode-positive.log" )
done
python3 /workspace/tests/m1b1/check_open_pressure_tracer_diffusion.py "${args[@]}"

# Deliberately falsify the positive-D input dictionary: checker must reject.
python3 - "$work/warped-positive/constant/transportProperties" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); s=p.read_text()
a="D [0 2 -1 0 0 0 0] 0.01;"
if s.count(a)!=1: raise SystemExit("FAIL: unable to inject D corruption")
p.write_text(s.replace(a,"D [0 2 -1 0 0 0 0] 0;"))
PY
if python3 /workspace/tests/m1b1/check_open_pressure_tracer_diffusion.py "${args[@]}" > "$work/negative.log" 2>&1; then
 echo 'FAIL: corrupted D input went undetected'; exit 35
fi
grep -q 'FAIL: wrong declared diffusivity' "$work/negative.log" || {
 cat "$work/negative.log"; echo 'FAIL: wrong negative control rejection'; exit 36;
}
echo 'PASS: deliberate positive-D dictionary corruption rejected'
echo 'PASS: M1B-1g-7 official v2606 3-D positive molecular diffusion and signed budgets'
