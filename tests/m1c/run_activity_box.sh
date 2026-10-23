#!/usr/bin/env bash
# M1C-0 real official OpenFOAM v2606 finite-volume HT/HTO reaction box.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1C setup_rc=$setup_rc version=${WM_PROJECT_VERSION:-unset}"
[[ "${WM_PROJECT_VERSION:-}" == v2606 || "${WM_PROJECT_VERSION:-}" == 2606 ]] || exit 18
for cmd in wmake wclean blockMesh checkMesh python3 timeout; do
 command -v "$cmd" >/dev/null || { echo "FAIL: missing $cmd"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1C rc=$rc work=$work"
 for f in "$work"/build.log "$work"/*/solve.log "$work"/*/mesh-check.log "$work"/negative-*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -65 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { tail -100 "$work/build.log"; exit 21; }
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then exit 22; fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
maker=/workspace/tests/m1c/make_activity_box.py
checker=/workspace/tests/m1c/check_activity_box.py
for entry in conversion:0.4:0 decay_conversion:0.4:0.2 decay_only:0:0.2 static:0:0; do
 IFS=: read -r label k lam <<< "$entry"
 case_dir="$work/$label"
 python3 "$maker" "$case_dir" "$k" "$lam" valid
 blockMesh -case "$case_dir" > "$case_dir/mesh.log" 2>&1 || { tail -80 "$case_dir/mesh.log"; exit 23; }
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || { tail -80 "$case_dir/mesh-check.log"; exit 24; }
 grep -q 'Mesh OK.' "$case_dir/mesh-check.log" || exit 24
 timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || { tail -100 "$case_dir/solve.log"; exit 31; }
 python3 "$checker" "$case_dir" "$case_dir/solve.log" "$k" "$lam"
done
cp "$work/decay_conversion/1/HT" "$work/pristine-HT"
python3 - "$work/decay_conversion/1/HT" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]); s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+24\s*\(\s*|internalField\s+uniform\s+)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: no final HT field to tamper")
s=s[:m.start(2)]+f"{float(m.group(2))+100:.17g}"+s[m.end(2):]
p.write_text(s)
PY
if python3 "$checker" "$work/decay_conversion" "$work/decay_conversion/solve.log" 0.4 0.2 > "$work/negative-tamper.log" 2>&1; then
 echo "FAIL: corrupted final HT passed"; exit 35
fi
grep -q 'FAIL: M1C' "$work/negative-tamper.log" || { cat "$work/negative-tamper.log"; exit 36; }
echo "PASS: deliberately corrupted HT field rejected"
for entry in negative-k negative-lambda bad-dims; do
 case_dir="$work/$entry"
 python3 "$maker" "$case_dir" 0.4 0.2 "$entry"
 blockMesh -case "$case_dir" > "$case_dir/mesh.log" 2>&1 || exit 23
 if timeout 30s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1; then
  echo "FAIL: invalid case $entry accepted"; exit 37
 fi
 if [[ "$entry" == bad-dims ]]; then
  grep -q 'M1C HT/HTO fields require Bq/m3 dimensions' "$case_dir/solve.log" || { tail -50 "$case_dir/solve.log"; exit 38; }
 else
  grep -q 'M1C reaction rates kOx and lambdaT must be non-negative' "$case_dir/solve.log" || { tail -50 "$case_dir/solve.log"; exit 38; }
 fi
 echo "PASS: invalid input $entry explicitly rejected"
done
echo "PASS: M1C-0 HT/HTO official v2606 reaction-box verification"
