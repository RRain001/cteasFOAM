#!/usr/bin/env bash
# M1C-9 extends immutable M1C-8 5-run suite with three independent directional refinements.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1C-8 official_setup_rc=$setup_rc"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for exe in python3 wmake wclean blockMesh checkMesh postProcess timeout; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c9-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then
 echo "FAIL: M1C-8 rc=$rc work=$work"
 for f in "$work"/build.log "$work"/*/solve.log "$work"/*/mesh-check.log; do
  if [[ -f "$f" ]]; then echo "== $f =="; tail -70 "$f"; fi
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
for mode in space_lo baseline space_hi time_lo time_hi space_x space_y space_z; do
 case "$mode" in
  space_lo) resolution=6x16x6;;
  space_hi) resolution=12x30x8;;
  space_x) resolution=12x20x6;;
  space_y) resolution=8x30x6;;
  space_z) resolution=8x20x8;;
  *) resolution=8x20x6;;
 esac
 case_dir="$work/$mode"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" warped "$resolution"
 blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || {
  tail -70 "$case_dir/blockMesh.log"; exit 23;
 }
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" warped "$resolution"
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || {
  tail -70 "$case_dir/mesh-check.log"; exit 24;
 }
 grep -q "Mesh OK." "$case_dir/mesh-check.log" || exit 25
 postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$case_dir/centres.log" 2>&1 || exit 26
 postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$case_dir/volumes.log" 2>&1 || exit 27
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$case_dir" warped "$resolution"
 python3 /workspace/tests/m1c/make_open_activity.py "$case_dir" warped
 if [[ "$mode" == space_x || "$mode" == space_y || "$mode" == space_z ]]; then
  python3 /workspace/tests/m1c/make_directional_case.py "$case_dir" "$mode"
 else
  python3 /workspace/tests/m1c/make_convergence_case.py "$case_dir" "$mode"
 fi
 start_ns=$(date +%s%N)
 timeout 300s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || {
  tail -90 "$case_dir/solve.log"; exit 31;
 }
 end_ns=$(date +%s%N)
 python3 - "$start_ns" "$end_ns" "$case_dir/wall-seconds.txt" <<'PYCLOCK'
import sys
from pathlib import Path
Path(sys.argv[3]).write_text(str((int(sys.argv[2])-int(sys.argv[1]))/1e9)+"\n")
PYCLOCK
 if [[ "$mode" == space_x || "$mode" == space_y || "$mode" == space_z ]]; then
  python3 /workspace/tests/m1c/check_receiver_diagnostics.py --single "$case_dir" "$mode"
 else
  python3 /workspace/tests/m1c/check_convergence.py --single "$case_dir" "$mode"
 fi
done
python3 /workspace/tests/m1c/check_convergence.py --matrix "$work"
python3 /workspace/tests/m1c/check_receiver_diagnostics.py --matrix "$work"
baseline="$work/baseline"
cp "$baseline/0.4/HT" "$work/pristine-HT"
python3 - "$baseline/0.4/HT" <<'PYTAMPER'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+960\s*\(\s*)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: baseline native HT cells missing")
p.write_text(s[:m.start(2)]+str(float(m.group(2))+1000)+s[m.end(2):])
PYTAMPER
if python3 /workspace/tests/m1c/check_convergence.py --single "$baseline" baseline > "$work/tamper-field.log" 2>&1; then
 echo "FAIL: M1C-8 corrupted physical HT accepted"; exit 40
fi
grep -q "FAIL: M1C-8" "$work/tamper-field.log" || {
 cat "$work/tamper-field.log"; exit 41;
}
mv "$work/pristine-HT" "$baseline/0.4/HT"
echo "PASS: M1C-8 tampered 960-cell native physical HT rejected"
cp "$baseline/solve.log" "$work/pristine-solve.log"
python3 - "$baseline/solve.log" <<'PYEM'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
s,n=re.subn(r"(emitted_HT_Bq=)([-+0-9.eE]+)",
            lambda m:m.group(1)+str(float(m.group(2))+1000),s,count=1)
if n!=1:raise SystemExit("FAIL: no actual OpenFOAM emitted HT to tamper")
p.write_text(s)
PYEM
if python3 /workspace/tests/m1c/check_convergence.py --single "$baseline" baseline > "$work/tamper-emission.log" 2>&1; then
 echo "FAIL: M1C-8 corrupted native emitted HT accepted"; exit 42
fi
grep -q "FAIL: M1C-8" "$work/tamper-emission.log" || {
 cat "$work/tamper-emission.log"; exit 43;
}
mv "$work/pristine-solve.log" "$baseline/solve.log"
echo "PASS: M1C-8 tampered native emitted HT Bq rejected"
echo "PASS: M1C-9 eight true native v2606 PISO cases, observation-operator/directional mesh diagnostics and strict negative controls"
