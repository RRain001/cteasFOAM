#!/usr/bin/env bash
# M1B-1g-12: unmodified original icoFoam vs cteas PISO, same warped periodic input.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
rc=$?
set -e
echo "M1B-1g-12 setup_rc=$rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for app in python3 blockMesh checkMesh postProcess wmake wclean icoFoam timeout; do
 command -v "$app" >/dev/null || { echo "FAIL: missing $app"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-upstream-periodic-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1B-1g-12 rc=$rc"
 for f in "$work"/*.log; do if [[ -f "$f" ]]; then echo "=== $f ==="; tail -75 "$f"; fi; done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { tail -90 "$work/build.log"; exit 21; }
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for item in coarse:8x20x6:960 fine:12x30x9:3240; do
 IFS=: read -r label dims count <<< "$item"
 ct="$work/ct-$label"; ref="$work/ico-$label"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$ct" zero "$dims"
 python3 /workspace/tests/m1b1/make_warped_piso_accelerating_advdiff.py prepare "$ct"
 python3 /workspace/tests/m1b1/make_warped_periodic_upstream_piso.py "$ct"
 blockMesh -case "$ct" > "$work/$label.mesh.log" 2>&1 || { tail -75 "$work/$label.mesh.log"; exit 23; }
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$ct" zero "$dims"
 checkMesh -case "$ct" -allGeometry -allTopology > "$work/$label.check.log" 2>&1 || { tail -75 "$work/$label.check.log"; exit 24; }
 grep -q 'Mesh OK.' "$work/$label.check.log" || exit 24
 postProcess -case "$ct" -func writeCellCentres -time 0 > "$work/$label.centres.log" 2>&1 || exit 25
 python3 /workspace/tests/m1b1/make_warped_piso_accelerating_advdiff.py seed "$ct" "$count"
 cp -a "$ct" "$ref"
 python3 - "$ref/system/controlDict" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); s=p.read_text()
old="application cteasAtmosphereFoam;"
if s.count(old)!=1: raise SystemExit("FAIL: unexpected controlDict application")
p.write_text(s.replace(old,"application icoFoam;"))
PY
 timeout 240s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$ct" > "$work/$label.native.log" 2>&1 || { tail -100 "$work/$label.native.log"; exit 31; }
 timeout 240s icoFoam -case "$ref" > "$work/$label.original.log" 2>&1 || { tail -100 "$work/$label.original.log"; exit 32; }
done
checker=/workspace/tests/m1b1/check_warped_periodic_upstream_piso.py
args=("$work/ct-coarse" "$work/ico-coarse" "$work/coarse.check.log" "$work/coarse.native.log" "$work/coarse.original.log" "$work/ct-fine" "$work/ico-fine" "$work/fine.check.log" "$work/fine.native.log" "$work/fine.original.log")
python3 "$checker" "${args[@]}"
cp "$work/ico-fine/0.4/U" "$work/original-U-pristine"
python3 - "$work/ico-fine/0.4/U" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<vector>\s+3240\s*\(\s*\(\s*)([-+0-9.eE]+)",s)
if not m:
 m=re.search(r"(internalField\s+uniform\s+\(\s*)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: missing original final U")
s=s[:m.start(2)]+f"{float(m.group(2))+0.02:.17g}"+s[m.end(2):]
p.write_text(s)
PY
if python3 "$checker" "${args[@]}" > "$work/negative-u.log" 2>&1; then
 echo "FAIL: original U tampering accepted"; exit 35
fi
grep -q 'FAIL: original U analytical uniformity' "$work/negative-u.log" || { cat "$work/negative-u.log"; exit 36; }
cp "$work/original-U-pristine" "$work/ico-fine/0.4/U"
cp "$work/ct-fine/0.4/c" "$work/native-c-pristine"
python3 - "$work/ct-fine/0.4/c" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+3240\s*\(\s*)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: missing native final scalar")
s=s[:m.start(2)]+f"{float(m.group(2))+0.08:.17g}"+s[m.end(2):]
p.write_text(s)
PY
if python3 "$checker" "${args[@]}" > "$work/negative-c.log" 2>&1; then
 echo "FAIL: native c tampering accepted"; exit 37
fi
grep -q 'FAIL: periodic PISO scalar analytical error' "$work/negative-c.log" || { cat "$work/negative-c.log"; exit 38; }
echo 'PASS: independent upstream U and native scalar negative controls rejected'
echo 'PASS: M1B-1g-12 official v2606 original unforced periodic PISO parity'
