#!/usr/bin/env bash
# M1C-12: source pulse matched to native 3-D source-free twins.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for exe in python3 wmake wclean blockMesh checkMesh postProcess timeout; do
 command -v "$exe" >/dev/null || { echo "FAIL: M1C12 missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c12-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then
 echo "FAIL: M1C-12 rc=$rc work=$work"
 for f in "$work"/build.log "$work"/*/mesh-check.log "$work"/*/solve.log "$work"/negative-*.log; do
  if [[ -f "$f" ]]; then echo "== $f =="; tail -45 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || {
 tail -110 "$work/build.log"; exit 21;
}
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for cfg in 8:coarse 8:normal 8:fine 12:normal 16:normal; do
 IFS=: read -r nx level <<< "$cfg"
 for source in control pulse; do
  case_dir="$work/$source"_"$nx"_"$level"
  resolution="$nx"x"$nx"x"$nx"
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" zero "$resolution"
  python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$case_dir" advect
  python3 /workspace/tests/m1c/make_pulsed_analytic.py "$case_dir" "$source" "$level"
  blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || { tail -80 "$case_dir/blockMesh.log"; exit 23; }
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" zero "$resolution"
  checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || { tail -90 "$case_dir/mesh-check.log"; exit 24; }
  grep -q 'Mesh OK.' "$case_dir/mesh-check.log" || exit 24
  postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$case_dir/centres.log" 2>&1 || { tail -70 "$case_dir/centres.log"; exit 25; }
  postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$case_dir/volumes.log" 2>&1 || { tail -70 "$case_dir/volumes.log"; exit 26; }
  python3 /workspace/tests/m1c/make_3d_analytic.py seed "$case_dir" "$((nx*nx*nx))"
  timeout 300s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || { tail -100 "$case_dir/solve.log"; exit 31; }
  python3 /workspace/tests/m1c/check_pulsed_analytic.py --single "$case_dir" "$source" "$nx:$level"
 done
done
checker=/workspace/tests/m1c/check_pulsed_analytic.py
python3 "$checker" --matrix "$work"
native="$work/pulse_16_normal"
cp "$native/0.2/HT" "$work/pristine-finalHT"
python3 - "$native/0.2/HT" <<'PYFIELD'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+4096\s*\(\s*)([-+0-9.eE]+)",s)
if m is None:raise SystemExit("FAIL: native final HT missing")
p.write_text(s[:m.start(2)]+f"{float(m.group(2))+100000:.17g}"+s[m.end(2):])
PYFIELD
if python3 "$checker" --single "$native" pulse 16:normal > "$work/negative-field.log" 2>&1; then
 echo "FAIL: M1C-12 corrupted native HT accepted"; exit 40
fi
grep -q 'FAIL: M1C-12 native pulse_16_normal continuum HT' "$work/negative-field.log" || { cat "$work/negative-field.log"; exit 41; }
mv "$work/pristine-finalHT" "$native/0.2/HT"
echo 'PASS: M1C-12 genuine native 4096-cell HT field corruption rejected'
cp "$native/solve.log" "$work/pristine-solve"
python3 - "$native/solve.log" <<'PYEM'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
s,n=re.subn(r"(emitted_HT_Bq=)([-+0-9.eE]+)",lambda m:m.group(1)+str(float(m.group(2))+1000),s,count=1)
if n!=1:raise SystemExit("FAIL: native emitted HT missing")
p.write_text(s)
PYEM
if python3 "$checker" --single "$native" pulse 16:normal > "$work/negative-emitted.log" 2>&1; then
 echo "FAIL: M1C-12 tampered emitted activity accepted"; exit 42
fi
grep -q 'FAIL: M1C-12 native emitted HT/HTO Bq or partial-step source overlap corrupted' "$work/negative-emitted.log" || { cat "$work/negative-emitted.log"; exit 43; }
mv "$work/pristine-solve" "$native/solve.log"
echo 'PASS: M1C-12 native emitted Bq log tampering rejected'
props="$native/constant/transportProperties"
cp "$props" "$work/pristine-transport"
python3 - "$props" <<'PYPROP'
from pathlib import Path
import sys
p=Path(sys.argv[1]);s=p.read_text()
old='sourceEnd [0 0 1 0 0 0 0] 0.13700000000000001;'
if s.count(old)!=1:raise SystemExit("FAIL: missing native source window")
p.write_text(s.replace(old,'sourceEnd [0 0 1 0 0 0 0] 0.14700000000000002;'))
PYPROP
if python3 "$checker" --single "$native" pulse 16:normal > "$work/negative-source-window.log" 2>&1; then
 echo "FAIL: M1C-12 altered native sourceEnd accepted"; exit 44
fi
grep -q 'FAIL: M1C-12 source window/rate or uniform shape changed:' "$work/negative-source-window.log" || { cat "$work/negative-source-window.log"; exit 45; }
mv "$work/pristine-transport" "$props"
echo 'PASS: M1C-12 real native pulse window modification rejected'
echo 'PASS: M1C-12 10 matched source/control OpenFOAM v2606 cases, true 3D continuum, BE time response, two nuclides and three negative controls'
