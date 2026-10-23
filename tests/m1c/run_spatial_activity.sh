#!/usr/bin/env bash
# M1C-1: two-species nonzero-U, positive-D, conversion/decay on real warped 3-D meshes.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1C-1 setup_rc=$setup_rc version=${WM_PROJECT_VERSION:-unset}"
[[ "${WM_PROJECT_VERSION:-}" == v2606 || "${WM_PROJECT_VERSION:-}" == 2606 ]] || exit 18
for app in python3 blockMesh checkMesh postProcess wmake wclean timeout; do
 command -v "$app" >/dev/null || { echo "FAIL: missing $app"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c-spatial-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1C-1 rc=$rc work=$work"
 for f in "$work"/build.log "$work"/*.mesh.log "$work"/*.check.log "$work"/*.volumes.log "$work"/*.solve.log "$work"/negative-*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -90 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { tail -100 "$work/build.log"; exit 21; }
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then exit 22; fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for entry in coarse:8x20x6:960 fine:12x30x9:3240; do
 IFS=: read -r label resolution count <<< "$entry"
 case_dir="$work/$label"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" zero "$resolution"
 python3 /workspace/tests/m1b1/make_warped_periodic_advection_diffusion.py prepare "$case_dir"
 python3 /workspace/tests/m1c/make_spatial_activity.py prepare "$case_dir"
 blockMesh -case "$case_dir" > "$work/$label.mesh.log" 2>&1 || { tail -90 "$work/$label.mesh.log"; exit 23; }
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" zero "$resolution"
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$work/$label.check.log" 2>&1 || { tail -90 "$work/$label.check.log"; exit 24; }
 grep -q 'Mesh OK.' "$work/$label.check.log" || exit 24
 postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$work/$label.centres.log" 2>&1 || { tail -90 "$work/$label.centres.log"; exit 25; }
 postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$work/$label.volumes.log" 2>&1 || { tail -90 "$work/$label.volumes.log"; exit 26; }
 test -f "$case_dir/0/V" || { echo 'FAIL: native cell volume field missing'; ls "$case_dir/0"; exit 26; }
 python3 /workspace/tests/m1c/make_spatial_activity.py seed "$case_dir" "$count"
 timeout 240s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$work/$label.solve.log" 2>&1 || { tail -90 "$work/$label.solve.log"; exit 31; }
done
checker=/workspace/tests/m1c/check_spatial_activity.py
python3 "$checker" "$work/coarse" "$work/coarse.check.log" "$work/coarse.solve.log" \
 "$work/fine" "$work/fine.check.log" "$work/fine.solve.log"
cp "$work/fine/0.4/HT" "$work/pristine-HT"
python3 - "$work/fine/0.4/HT" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]); s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+3240\s*\(\s*)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: real final HT field missing")
s=s[:m.start(2)]+f"{float(m.group(2))+100000:.17g}"+s[m.end(2):]
p.write_text(s)
PY
if python3 "$checker" "$work/coarse" "$work/coarse.check.log" "$work/coarse.solve.log" \
 "$work/fine" "$work/fine.check.log" "$work/fine.solve.log" > "$work/negative-field.log" 2>&1; then
 echo 'FAIL: corrupted final HT passed'; exit 35
fi
grep -q 'FAIL: M1C-1' "$work/negative-field.log" || { cat "$work/negative-field.log"; exit 36; }
echo 'PASS: deliberate final HT corruption rejected'
cp "$work/pristine-HT" "$work/fine/0.4/HT"
python3 - "$work/fine.solve.log" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]); s=p.read_text()
m=re.search(r"(flux_HT_net_Bq_s=)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: native HT flux record missing")
s=s[:m.start(2)]+f"{float(m.group(2))+100000:.17g}"+s[m.end(2):]
p.write_text(s)
PY
if python3 "$checker" "$work/coarse" "$work/coarse.check.log" "$work/coarse.solve.log" \
 "$work/fine" "$work/fine.check.log" "$work/fine.solve.log" > "$work/negative-flux.log" 2>&1; then
 echo 'FAIL: corrupted HT flux passed'; exit 37
fi
grep -q 'FAIL: M1C-1' "$work/negative-flux.log" || { cat "$work/negative-flux.log"; exit 38; }
echo 'PASS: deliberate HT boundary flux corruption rejected'
echo 'PASS: M1C-1 official v2606 3-D warped periodic HT/HTO advection-diffusion-reaction'
