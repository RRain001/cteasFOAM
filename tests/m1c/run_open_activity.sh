#!/usr/bin/env bash
# M1C-3 real OpenFOAM v2606 warped duct open-boundary HT/HTO budget.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1C-3 setup_rc=$setup_rc version=${WM_PROJECT_VERSION:-unset}"
[[ "${WM_PROJECT_VERSION:-}" == v2606 || "${WM_PROJECT_VERSION:-}" == 2606 ]] || exit 18
for exe in python3 wmake wclean blockMesh checkMesh postProcess timeout; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c-open-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "M1C-3 FAIL rc=$rc evidence at $work";
 for f in "$work"/build.log "$work"/*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -80 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || {
 tail -100 "$work/build.log"; exit 21;
}
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for mode in warped reverse zero; do
 case_dir="$work/$mode"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" "$mode"
 if [[ "$mode" == warped ]]; then
  blockMesh -case "$case_dir" > "$work/mesh.log" 2>&1 || {
   tail -90 "$work/mesh.log"; exit 23; }
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" "$mode"
  checkMesh -case "$case_dir" -allGeometry -allTopology > "$work/check.log" 2>&1 || {
   tail -90 "$work/check.log"; exit 24; }
  grep -q "Mesh OK." "$work/check.log" || exit 24
  postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$work/centres.log" 2>&1 || {
   tail -90 "$work/centres.log"; exit 25; }
  postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$work/volumes.log" 2>&1 || {
   tail -90 "$work/volumes.log"; exit 26; }
 else
  cp -a "$work/warped/constant/polyMesh" "$case_dir/constant/"
  cp -a "$work/warped/0/C" "$case_dir/0/"
  cp -a "$work/warped/0/V" "$case_dir/0/"
 fi
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$case_dir" "$mode"
 python3 /workspace/tests/m1c/make_open_activity.py "$case_dir" "$mode"
 timeout 240s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$work/solve-$mode.log" 2>&1 || {
  tail -110 "$work/solve-$mode.log"; exit 31; }
 echo "M1C-3 mode=$mode steps=$(grep -c '^Time = ' "$work/solve-$mode.log") psolves=$(grep -c 'Solving for p,' "$work/solve-$mode.log") HTsolves=$(grep -c 'Solving for HT,' "$work/solve-$mode.log") HTOsolves=$(grep -c 'Solving for HTO,' "$work/solve-$mode.log")"
done
checker=/workspace/tests/m1c/check_open_activity.py
args=( "$work/check.log" "$work/warped" "$work/solve-warped.log"
       "$work/reverse" "$work/solve-reverse.log"
       "$work/zero" "$work/solve-zero.log" )
python3 "$checker" "${args[@]}"
cp "$work/solve-warped.log" "$work/corrupt-flux.log"
python3 - "$work/corrupt-flux.log" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]); s=p.read_text()
s,n=re.subn(r"(flux_HT_net_Bq_s=)([-+0-9.eE]+)",r"\g<1>100000",s,count=1)
if n!=1: raise SystemExit("FAIL: no actual HT flux to tamper")
p.write_text(s)
PY
if python3 "$checker" "$work/check.log" "$work/warped" "$work/corrupt-flux.log" \
 "$work/reverse" "$work/solve-reverse.log" "$work/zero" "$work/solve-zero.log" > "$work/negative-flux.log" 2>&1; then
 echo "FAIL: corrupted real HT flux undetected"; exit 35
fi
grep -q "FAIL: M1C-3" "$work/negative-flux.log" || {
 cat "$work/negative-flux.log"; exit 36; }
echo "PASS: deliberate HT signed flux corruption rejected"
cp "$work/warped/0.4/HT" "$work/real-HT"
python3 - "$work/warped/0.4/HT" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+960\s*\(\s*)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: actual final native nonuniform HT not found")
s=s[:m.start(2)]+str(float(m.group(2))+100000)+s[m.end(2):]
p.write_text(s)
PY
if python3 "$checker" "${args[@]}" > "$work/negative-field.log" 2>&1; then
 echo "FAIL: corrupted final HT field undetected"; exit 37
fi
grep -q "FAIL: M1C-3" "$work/negative-field.log" || {
 cat "$work/negative-field.log"; exit 38; }
echo "PASS: deliberate final HT field corruption rejected"
echo "PASS: M1C-3 open-boundary HT/HTO real PISO and signed activity budget"
