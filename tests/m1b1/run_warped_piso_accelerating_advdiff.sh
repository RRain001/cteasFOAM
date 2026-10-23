#!/usr/bin/env bash
# M1B-1g-11: actual PISO momentum + scalar on warped 3D periodic-x meshes.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B-1g-11 setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for app in python3 blockMesh checkMesh postProcess wmake wclean timeout; do
 command -v "$app" >/dev/null || { echo "FAIL: missing $app"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-piso-accelerating-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1B-1g-11 rc=$rc"
 for f in "$work"/build.log "$work"/coarse.mesh.log "$work"/coarse.check.log "$work"/coarse.solve.log "$work"/fine.mesh.log "$work"/fine.check.log "$work"/fine.solve.log "$work"/negative-c.log "$work"/negative-u.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -85 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { tail -90 "$work/build.log"; exit 21; }
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then exit 22; fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for entry in coarse:8x20x6:960 fine:12x30x9:3240; do
 IFS=: read -r label resolution count <<< "$entry"
 case_dir="$work/$label"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" zero "$resolution"
 python3 /workspace/tests/m1b1/make_warped_piso_accelerating_advdiff.py prepare "$case_dir"
 blockMesh -case "$case_dir" > "$work/$label.mesh.log" 2>&1 || { tail -90 "$work/$label.mesh.log"; exit 23; }
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" zero "$resolution"
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$work/$label.check.log" 2>&1 || { tail -90 "$work/$label.check.log"; exit 24; }
 grep -q 'Mesh OK.' "$work/$label.check.log" || exit 24
 postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$work/$label.centres.log" 2>&1 || { tail -90 "$work/$label.centres.log"; exit 25; }
 python3 /workspace/tests/m1b1/make_warped_piso_accelerating_advdiff.py seed "$case_dir" "$count"
 timeout 240s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$work/$label.solve.log" 2>&1 || { tail -90 "$work/$label.solve.log"; exit 31; }
done
checker=/workspace/tests/m1b1/check_warped_piso_accelerating_advdiff.py
python3 "$checker" "$work/coarse" "$work/coarse.check.log" "$work/coarse.solve.log" "$work/fine" "$work/fine.check.log" "$work/fine.solve.log"
cp "$work/fine/0.4/c" "$work/fine/c-pristine"
python3 - "$work/fine/0.4/c" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]); s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+3240\s*\(\s*)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: fine-grid native c field absent")
s=s[:m.start(2)]+f"{float(m.group(2))+0.08:.17g}"+s[m.end(2):]
p.write_text(s)
PY
if python3 "$checker" "$work/coarse" "$work/coarse.check.log" "$work/coarse.solve.log" "$work/fine" "$work/fine.check.log" "$work/fine.solve.log" > "$work/negative-c.log" 2>&1; then
 echo 'FAIL: corrupted native c not detected'; exit 35
fi
grep -q 'FAIL: PISO scalar analytical error' "$work/negative-c.log" || { cat "$work/negative-c.log"; exit 36; }
cp "$work/fine/c-pristine" "$work/fine/0.4/c"
cp "$work/fine/0.4/U" "$work/fine/U-pristine"
python3 - "$work/fine/0.4/U" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]); s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<vector>\s+3240\s*\(\s*\(\s*)([-+0-9.eE]+)",s)
if m:
    s=s[:m.start(2)]+f"{float(m.group(2))+0.03:.17g}"+s[m.end(2):]
else:
    m=re.search(r"(internalField\s+uniform\s+\(\s*)([-+0-9.eE]+)",s)
    if not m: raise SystemExit("FAIL: fine-grid native U field absent")
    s=s[:m.start(2)]+f"{float(m.group(2))+0.03:.17g}"+s[m.end(2):]
p.write_text(s)
PY
if python3 "$checker" "$work/coarse" "$work/coarse.check.log" "$work/coarse.solve.log" "$work/fine" "$work/fine.check.log" "$work/fine.solve.log" > "$work/negative-u.log" 2>&1; then
 echo 'FAIL: corrupted native U not detected'; exit 37
fi
grep -q 'FAIL: PISO-computed final U error' "$work/negative-u.log" || { cat "$work/negative-u.log"; exit 38; }
echo 'PASS: native final scalar +0.08 and U +0.03 tampering independently rejected'
echo 'PASS: M1B-1g-11 official v2606 warped 3D accelerating PISO advection-diffusion'
