#!/usr/bin/env bash
# M1C-11: independent 3-D source-free trigonometric eigenmode / native FVM test.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for exe in python3 wmake wclean blockMesh checkMesh postProcess timeout; do
 command -v "$exe" >/dev/null || { echo "FAIL: M1C-11 missing native tool $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c11-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then
 echo "FAIL: M1C-11 rc=$rc work=$work"
 for f in "$work"/build.log "$work"/*/mesh-check.log "$work"/*/solve.log "$work"/negative-*.log; do
  if [[ -f "$f" ]]; then echo "== $f =="; tail -60 "$f"; fi
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
for mode in advect diffuse; do
 for nx in 8 12 16; do
  case_dir="$work/$mode"_"$nx"
  resolution="$nx"x"$nx"x"$nx"
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" zero "$resolution"
  python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$case_dir" "$mode"
  blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || {
   tail -90 "$case_dir/blockMesh.log"; exit 23;
  }
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" zero "$resolution"
  checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || {
   tail -90 "$case_dir/mesh-check.log"; exit 24;
  }
  grep -q "Mesh OK." "$case_dir/mesh-check.log" || exit 24
  postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$case_dir/centres.log" 2>&1 || {
   tail -90 "$case_dir/centres.log"; exit 25;
  }
  postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$case_dir/volumes.log" 2>&1 || {
   tail -90 "$case_dir/volumes.log"; exit 26;
  }
  python3 /workspace/tests/m1c/make_3d_analytic.py seed "$case_dir" "$((nx*nx*nx))"
  timeout 300s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || {
   tail -90 "$case_dir/solve.log"; exit 31;
  }
  python3 /workspace/tests/m1c/check_3d_analytic.py --single "$case_dir" "$mode" "$nx"
 done
done
python3 /workspace/tests/m1c/check_3d_analytic.py --matrix "$work"
checker=/workspace/tests/m1c/check_3d_analytic.py
native="$work/advect_16"
cp "$native/0.2/HT" "$work/pristine-final-ht"
python3 - "$native/0.2/HT" <<'PYFIELD'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+4096\s*\(\s*)([-+0-9.eE]+)",s)
if m is None:raise SystemExit("FAIL: missing real native HT values for negative control")
p.write_text(s[:m.start(2)]+f"{float(m.group(2))+100000:.17g}"+s[m.end(2):])
PYFIELD
if python3 "$checker" --single "$native" advect 16 > "$work/negative-field.log" 2>&1; then
 echo "FAIL: M1C-11 corrupted native HT field incorrectly accepted"; exit 40
fi
grep -Eq 'FAIL: M1C-11 (actual native field Bq integral differs|advect 16 HT analytic RMS=)' "$work/negative-field.log" || {
 cat "$work/negative-field.log"; exit 41;
}
mv "$work/pristine-final-ht" "$native/0.2/HT"
echo "PASS: M1C-11 genuine native 4096-cell HT field corruption rejected by independent analytic/field-Bq audit"
cp "$native/solve.log" "$work/pristine-solve"
python3 - "$native/solve.log" <<'PYFLUX'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
s,num=re.subn(r"(flux_HT_net_Bq_s=)([-+0-9.eE]+)",
    lambda m:m.group(1)+str(float(m.group(2))+100000),s,count=1)
if num!=1:raise SystemExit("FAIL: missing native signed HT flux for negative control")
p.write_text(s)
PYFLUX
if python3 "$checker" --single "$native" advect 16 > "$work/negative-flux.log" 2>&1; then
 echo "FAIL: M1C-11 corrupted native HT flux accepted"; exit 42
fi
grep -q 'FAIL: M1C-11 three dimensional cyclic net face flux not balanced' "$work/negative-flux.log" || {
 cat "$work/negative-flux.log"; exit 43;
}
mv "$work/pristine-solve" "$native/solve.log"
echo "PASS: M1C-11 genuine native HT signed-face flux corruption rejected"
cp "$native/0/U" "$work/pristine-U"
python3 - "$native/0/U" <<'PYU'
from pathlib import Path
import sys
p=Path(sys.argv[1]);s=p.read_text()
part="internalField uniform (0.2 0.1 0.07);"
if s.count(part)!=1:raise SystemExit("FAIL: no genuine native 3D U to tamper")
p.write_text(s.replace(part,"internalField uniform (0.2 0.2 0.07);"))
PYU
if python3 "$checker" --single "$native" advect 16 > "$work/negative-U.log" 2>&1; then
 echo "FAIL: M1C-11 altered native U accepted"; exit 44
fi
grep -q 'FAIL: M1C-11 actual native prescribed U does not match frozen mode' "$work/negative-U.log" || {
 cat "$work/negative-U.log"; exit 45;
}
mv "$work/pristine-U" "$native/0/U"
echo "PASS: M1C-11 genuine native 3-D prescribed U corruption rejected"
echo "PASS: M1C-11 six native source-free 3-D eigenmode cases, analytic convergence, signed Bq, three independent negative controls"
