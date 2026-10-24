#!/usr/bin/env bash
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c14-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then echo "FAIL: M1C-14 exit=$rc"; for f in "$work"/build.log "$work"/*/solve.log "$work"/negative-*.log; do if [[ -f "$f" ]]; then echo "=== $f ==="; tail -32 "$f"; fi; done; fi; rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere" && wclean && wmake > "$work/build.log" 2>&1) || { tail -80 "$work/build.log"; exit 21; }
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam"
for n in 8 12 16 24; do
 case_dir="$work/$n"_normal
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" zero "$n"x"$n"x"$n"
 python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$case_dir" advect
 python3 /workspace/tests/m1c/make_gaussian_green.py "$case_dir" normal
 blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || { tail -70 "$case_dir/blockMesh.log"; exit 23; }
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" zero "$n"x"$n"x"$n"
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || { tail -90 "$case_dir/mesh-check.log"; exit 24; }
 grep -q 'Mesh OK.' "$case_dir/mesh-check.log"
 postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$case_dir/centres.log" 2>&1 || { tail -50 "$case_dir/centres.log"; exit 25; }
 postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$case_dir/volumes.log" 2>&1 || { tail -50 "$case_dir/volumes.log"; exit 26; }
 timeout 420s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || { tail -70 "$case_dir/solve.log"; exit 31; }
 python3 /workspace/tests/m1c/check_fixed_receptors.py --single "$case_dir" "$n"
done
python3 /workspace/tests/m1c/check_fixed_receptors.py --matrix "$work"
native="$work/16_normal"
checker=/workspace/tests/m1c/check_fixed_receptors.py
cp "$native/0.4/HT" "$work/pristine-final-HT"
python3 - "$native" <<'PYHT'
from pathlib import Path
import re,sys
sys.path.insert(0,"/workspace/tests/m1c")
from check_fixed_receptors import grid,weights,POS
p=Path(sys.argv[1])
c,_,lookup=grid(p,16)
anchors,_=weights(POS["R2"],c,lookup,16)
ix=max(anchors,key=lambda pair:pair[1])[0]
f=p/"0.4/HT"
s=f.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+4096\s*\()(.*?)\)\s*;",s,re.S)
assert m is not None
v=m.group(2).split()
assert len(v)==4096
v[ix]=str(float(v[ix])+100000)
f.write_text(s[:m.start(2)]+"\n"+"\n".join(v)+"\n"+s[m.end(2):])
print("M1C14_TAMPER real HT nonzero-weight anchor",ix)
PYHT
if python3 "$checker" --single "$native" 16 > "$work/negative-HT.log" 2>&1; then echo "FAIL: native HT tamper accepted"; exit 41; fi
grep -q 'FAIL: M1C-14 native final source field integrated Bq differs from solver log HT' "$work/negative-HT.log" || { cat "$work/negative-HT.log"; exit 42; }
mv "$work/pristine-final-HT" "$native/0.4/HT"
echo "PASS: M1C14 tampered native nonzero-weight HT anchor rejected"
cp "$native/0/C" "$work/pristine-centres"
python3 - "$native/0/C" <<'PYC'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<vector>\s+4096\s*\()(.*?)\)\s*;",s,re.S)
assert m is not None
a=re.search(r"\(([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\)",m.group(2))
assert a
replacement=f"({a.group(1)} {float(a.group(2))+0.1:.16g} {a.group(3)})"
p.write_text(s[:m.start(2)]+m.group(2)[:a.start()]+replacement+m.group(2)[a.end():]+s[m.end(2):])
PYC
if python3 "$checker" --single "$native" 16 > "$work/negative-centres.log" 2>&1; then echo "FAIL: native coordinates tamper accepted"; exit 43; fi
grep -q 'FAIL: M1C-14 native physical cell-centre y/z lattice geometry corrupted' "$work/negative-centres.log" || { cat "$work/negative-centres.log"; exit 44; }
mv "$work/pristine-centres" "$native/0/C"
echo "PASS: M1C14 tampered native physical cell-centre geometry rejected"
cp "$native/constant/transportProperties" "$work/pristine-props"
sed -i 's/sourceSigma \[0 1 0 0 0 0 0\] 0.08;/sourceSigma [0 1 0 0 0 0 0] 0.12;/' "$native/constant/transportProperties"
if python3 "$checker" --single "$native" 16 > "$work/negative-source.log" 2>&1; then echo "FAIL: source sigma tamper accepted"; exit 45; fi
grep -q 'FAIL: M1C-14 native physical source/wind/reactive parameter modified sourceSigma' "$work/negative-source.log" || { cat "$work/negative-source.log"; exit 46; }
mv "$work/pristine-props" "$native/constant/transportProperties"
echo "PASS: M1C14 tampered real Gaussian source width rejected"
echo "PASS: M1C-14 four native 3D meshes, 24-cube refinement, fixed receptor affine interpolation and error split, strict negative controls"
