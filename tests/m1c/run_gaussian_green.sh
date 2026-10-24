#!/usr/bin/env bash
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for e in python3 wmake wclean blockMesh checkMesh postProcess timeout; do command -v "$e" >/dev/null || exit 18; done
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c13-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then echo "FAIL: M1C-13 rc=$rc work=$work"; for f in "$work"/build.log "$work"/*/mesh-check.log "$work"/*/solve.log "$work"/negative-*.log; do if [[ -f "$f" ]]; then echo "== $f =="; tail -42 "$f"; fi; done; fi; rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { tail -100 "$work/build.log"; exit 21; }
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for cfg in 8:coarse 8:normal 12:normal 16:normal; do
 IFS=: read -r nx level <<< "$cfg"
 case_dir="$work/$nx"_"$level"
 resolution="$nx"x"$nx"x"$nx"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" zero "$resolution"
 python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$case_dir" advect
 python3 /workspace/tests/m1c/make_gaussian_green.py "$case_dir" "$level"
 blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || { tail -80 "$case_dir/blockMesh.log"; exit 23; }
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" zero "$resolution"
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || { tail -80 "$case_dir/mesh-check.log"; exit 24; }
 grep -q 'Mesh OK.' "$case_dir/mesh-check.log" || exit 24
 postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$case_dir/centres.log" 2>&1 || exit 25
 postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$case_dir/volumes.log" 2>&1 || exit 26
 timeout 300s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || { tail -100 "$case_dir/solve.log"; exit 31; }
 python3 /workspace/tests/m1c/check_gaussian_green.py --single "$case_dir" "$nx" "$level"
done
checker=/workspace/tests/m1c/check_gaussian_green.py
python3 "$checker" --matrix "$work"
native="$work/16_normal"
cp "$native/0.4/HT" "$work/pristine-HT"
python3 - "$native" <<'PYFIELD'
from pathlib import Path
import sys,re
sys.path.insert(0,"/workspace/tests/m1c")
from check_gaussian_green import field,nearest,POS
case=Path(sys.argv[1])
c=field(case/"0/C",4096,"vector","0 1 0 0 0 0 0")
i,_=nearest(c,POS["R2"])
p=case/"0.4/HT";s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+4096\s*\()(.*?)\)\s*;",s,re.S)
if not m:raise SystemExit("FAIL: missing real native 4096 HT values")
values=m.group(2).split()
if len(values)!=4096:raise SystemExit("FAIL: missing actual real native cell field")
values[i]=f"{float(values[i])+100000:.17g}"
p.write_text(s[:m.start(2)]+"\n"+"\n".join(values)+"\n"+s[m.end(2):])
print(f"M1C13_TAMPER actual R2-nearest native HT cell index={i}")
PYFIELD
if python3 "$checker" --single "$native" 16 normal > "$work/negative-field.log" 2>&1; then echo "FAIL: M1C-13 tampered native field accepted"; exit 40; fi
grep -q 'FAIL: M1C-13 final native species Bq inventory does not match real physical field HT' "$work/negative-field.log" || { cat "$work/negative-field.log"; exit 41; }
mv "$work/pristine-HT" "$native/0.4/HT"
echo "PASS: M1C-13 actual R2 native receptor-cell corruption rejected on independent inventory"
cp "$native/solve.log" "$work/pristine-solve"
python3 - "$native/solve.log" <<'PYEM'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
s,n=re.subn(r"(emitted_HT_Bq=)([-+0-9.eE]+)",lambda m:m.group(1)+str(float(m.group(2))+1000),s,count=1)
if n!=1:raise SystemExit("FAIL: genuine native emitted HT missing")
p.write_text(s)
PYEM
if python3 "$checker" --single "$native" 16 normal > "$work/negative-emitted.log" 2>&1; then echo "FAIL: M1C-13 tampered Gaussian emission accepted"; exit 42; fi
grep -q 'FAIL: M1C-13 native Gaussian source emitted HT/HTO Bq or time overlap corrupted' "$work/negative-emitted.log" || { cat "$work/negative-emitted.log"; exit 43; }
mv "$work/pristine-solve" "$native/solve.log"
echo "PASS: M1C-13 actual native emitted Gaussian HT Bq tampering rejected"
props="$native/constant/transportProperties"
cp "$props" "$work/pristine-props"
python3 - "$props" <<'PYPROP'
from pathlib import Path
import sys
p=Path(sys.argv[1]);s=p.read_text()
old="sourceSigma [0 1 0 0 0 0 0] 0.08;"
if s.count(old)!=1:raise SystemExit("FAIL: actual native Gaussian sigma missing")
p.write_text(s.replace(old,"sourceSigma [0 1 0 0 0 0 0] 0.12;"))
PYPROP
if python3 "$checker" --single "$native" 16 normal > "$work/negative-sigma.log" 2>&1; then echo "FAIL: M1C-13 altered Gaussian width accepted"; exit 44; fi
grep -q 'FAIL: M1C-13 physical Gaussian source/clock config changed sourceSigma' "$work/negative-sigma.log" || { cat "$work/negative-sigma.log"; exit 45; }
mv "$work/pristine-props" "$props"
echo "PASS: M1C-13 actual native Gaussian sigma tampering rejected"
echo "PASS: M1C-13 four native v2606 spatial Gaussian plume runs, periodic Green semi-analytic references, Bq and three negative controls"
