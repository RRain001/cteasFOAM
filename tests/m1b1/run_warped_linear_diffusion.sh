#!/usr/bin/env bash
# M1B-1g-8: analytic affine c=1-x on genuine 3-D nonorthogonal mesh.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B-1g-8 setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for app in python3 blockMesh checkMesh postProcess wmake wclean timeout; do
 command -v "$app" >/dev/null || { echo "FAIL: missing $app"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-affine-diffusion-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1B-1g-8 rc=$rc"
 for f in "$work"/build.log "$work"/mesh.log "$work"/check.log "$work"/centres.log "$work"/solve.log "$work"/negative.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -90 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { tail -90 "$work/build.log"; exit 21; }
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then exit 22; fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
case_dir="$work/affine"
python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" zero
blockMesh -case "$case_dir" > "$work/mesh.log" 2>&1 || { tail -90 "$work/mesh.log"; exit 23; }
python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" zero
checkMesh -case "$case_dir" -allGeometry -allTopology > "$work/check.log" 2>&1 || { tail -90 "$work/check.log"; exit 24; }
grep -q 'Mesh OK.' "$work/check.log" || exit 24
postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$work/centres.log" 2>&1 || { tail -90 "$work/centres.log"; exit 25; }
python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$case_dir" zero
python3 /workspace/tests/m1b1/make_warped_linear_diffusion.py "$case_dir"
timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$work/solve.log" 2>&1 || { tail -90 "$work/solve.log"; exit 31; }
python3 /workspace/tests/m1b1/check_warped_linear_diffusion.py "$case_dir" "$work/check.log" "$work/solve.log"
# Fault injection on the actual written scalar output: must be rejected.
python3 - "$case_dir/0.4/c" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1])
s=p.read_text()
pat=r"(internalField\s+nonuniform\s+List<scalar>\s+960\s*\(\s*)([-+0-9.eE]+)"
m=re.search(pat,s)
if not m: raise SystemExit("FAIL: native final scalar not found for corruption")
s=s[:m.start(2)]+f"{float(m.group(2)) + 0.05:.17g}"+s[m.end(2):]
p.write_text(s)
PY
if python3 /workspace/tests/m1b1/check_warped_linear_diffusion.py "$case_dir" "$work/check.log" "$work/solve.log" > "$work/negative.log" 2>&1; then
 echo 'FAIL: corrupted final scalar field went undetected'; exit 35
fi
grep -q 'FAIL: analytical c=1-x mismatch' "$work/negative.log" || { cat "$work/negative.log"; exit 36; }
echo 'PASS: deliberate +0.05 final scalar corruption rejected'
echo 'PASS: M1B-1g-8 official OpenFOAM v2606 analytic affine diffusion on warped mesh'
