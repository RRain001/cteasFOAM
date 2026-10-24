#!/usr/bin/env bash
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c18-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then echo "FAIL: M1C18 rc=$rc"; for f in "$work"/build.log "$work"/*/solve.log; do if [[ -f "$f" ]]; then echo "== $f =="; tail -65 "$f";fi;done;fi; rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere" && wclean && wmake > "$work/build.log" 2>&1) || { tail -95 "$work/build.log"; exit 20; }
exe="$FOAM_USER_APPBIN/cteasAtmosphereFoam"
test -x "$exe" || exit 21
for n in 8 16 24; do
  for mode in base blend; do
    c="$work/${n}_$mode"
    python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$c" zero "${n}x${n}x${n}"
    python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$c" advect
    python3 /workspace/tests/m1c/make_gaussian_green.py "$c" normal
    python3 - "$c" "$mode" <<'PYCASE'
from pathlib import Path
import sys
case=Path(sys.argv[1]);mode=sys.argv[2]
p=case/'constant/transportProperties'
p.write_text(p.read_text()+'auditActivityExtrema true;\nenforceNonnegativeActivity true;\n'+('activityConvectionBlend globalPositive;\n' if mode=='blend' else ''))
if mode=='blend':
 f=case/'system/fvSchemes';txt=f.read_text();needle='div(phi,HTO) Gauss upwind;'
 assert txt.count(needle)==1
 txt=txt.replace(needle,needle+' div(phi,HTHigh) Gauss limitedLinear 1; div(phi,HTOHigh) Gauss limitedLinear 1;')
 f.write_text(txt)
print(f'M1C18_REAL_NATIVE_CASE n={sys.argv[1]} mode={mode}')
PYCASE
    blockMesh -case "$c" > "$c/blockMesh.log" 2>&1 || { tail -40 "$c/blockMesh.log"; exit 22; }
    python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$c" zero "${n}x${n}x${n}"
    checkMesh -case "$c" -allGeometry -allTopology > "$c/mesh-check.log" 2>&1 || { tail -40 "$c/mesh-check.log"; exit 23; }
    grep -q 'Mesh OK.' "$c/mesh-check.log" || exit 24
    postProcess -case "$c" -func writeCellCentres -time 0 > "$c/centres.log" 2>&1 || exit 25
    postProcess -case "$c" -func writeCellVolumes -time 0 > "$c/volumes.log" 2>&1 || exit 26
    timeout 600s "$exe" -case "$c" > "$c/solve.log" 2>&1 || { tail -80 "$c/solve.log"; exit 31; }
  done
  python3 /workspace/tests/m1c/check_flux_blend.py --compare "$work" "$n"
done
# Reject invalid named high-order advection without pretending a normal solution.
case_dir="$work/8_blend"
sed -i 's/div(phi,HTHigh) Gauss limitedLinear 1;/div(phi,HTHigh) Gauss unrecognizedHighOrderScheme;/' "$case_dir/system/fvSchemes"
if timeout 60s "$exe" -case "$case_dir" > "$work/invalid-scheme.log" 2>&1; then
  echo 'FAIL: M1C18 invalid FV high-order scheme accepted'; exit 40
fi
grep -q 'unrecognizedHighOrderScheme' "$work/invalid-scheme.log" || { tail -40 "$work/invalid-scheme.log"; exit 41; }
echo 'PASS: M1C-18 invalid named native high-order scheme rejected'
# Unsupported non-periodic geometry should fail before any timestep.
case_dir="$work/8_base"
sed -i 's/enforceNonnegativeActivity true;/enforceNonnegativeActivity true; activityConvectionBlend globalPositive;/' "$case_dir/constant/transportProperties"
# Existing 8^3 periodic field remains; test invalid activity mode as negative control.
sed -i 's/activityConvectionBlend globalPositive;/activityConvectionBlend nonexistent;/' "$case_dir/constant/transportProperties"
if timeout 60s "$exe" -case "$case_dir" > "$work/invalid-mode.log" 2>&1; then
  echo 'FAIL: M1C18 unsupported mode accepted'; exit 42
fi
grep -q 'Unsupported activityConvectionBlend' "$work/invalid-mode.log" || { tail -40 "$work/invalid-mode.log"; exit 43; }
# Confirm that an actual two-rank decomposed case rejects unsupported MPI
# through the native M1C-18 gate (not a failed mesh/openfoam startup).
parallel="$work/parallel-blend"
cp -a "$work/8_blend" "$parallel"
cat > "$parallel/system/decomposeParDict" <<'DICT'
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
distributed no;
roots ();
DICT
decomposePar -case "$parallel" -time 0 > "$work/decompose.log" 2>&1 || { tail -60 "$work/decompose.log"; exit 44; }
for rank in 0 1; do test -f "$parallel/processor$rank/constant/polyMesh/owner" || exit 45; done
if OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 60s mpirun --oversubscribe -np 2 \
 "$exe" -parallel -case "$parallel" > "$work/mpi-reject.log" 2>&1; then
 echo 'FAIL: M1C18 unsupported MPI accepted'; exit 46
fi
grep -q 'M1C18_GLOBAL_BLEND_UNSUPPORTED requires serial prescribed flow' "$work/mpi-reject.log" || { tail -50 "$work/mpi-reject.log"; exit 47; }
echo 'PASS: M1C-18 genuine decomposed 2-rank MPI input rejected by experimental serial-only gate'
# Verify a forged alpha cannot be mistaken for native physical evidence.
cp "$work/8_blend/solve.log" "$work/native-pristine-alpha.log"
python3 - "$work/8_blend/solve.log" <<'PYTAMPER'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
s,n=re.subn(r'(CTEAS_M1C18 species=HT global_alpha=)([-+0-9.eE]+)',r'\g<1>-0.25',s,count=1)
assert n==1
p.write_text(s)
PYTAMPER
if python3 /workspace/tests/m1c/check_flux_blend.py --compare "$work" 8 > "$work/tampered-alpha-check.log" 2>&1; then
 echo 'FAIL: M1C18 negative forged global alpha accepted'; exit 48
fi
grep -q 'FAIL: M1C18 invalid native alpha' "$work/tampered-alpha-check.log" || { tail -35 "$work/tampered-alpha-check.log"; exit 49; }
mv "$work/native-pristine-alpha.log" "$work/8_blend/solve.log"
echo 'PASS: M1C-18 nonphysical negative global-alpha log rejected by independent checker'
# Nonperiodic geometry: genuine native patch/wall 3-D case, not a string-only mock.
open="$work/open-nonperiodic"
python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$open" zero 8x8x8
blockMesh -case "$open" > "$work/open-blockMesh.log" 2>&1 || exit 50
python3 /workspace/tests/m1c/make_open_activity.py "$open" zero
python3 - "$open/constant/transportProperties" <<'PYOPEN'
from pathlib import Path
import sys
p=Path(sys.argv[1]);s=p.read_text();assert s.count('momentumMode pisoLaminar;')==1
p.write_text(s.replace('momentumMode pisoLaminar;','momentumMode prescribed;')+'activityConvectionBlend globalPositive;\n')
PYOPEN
if timeout 60s "$exe" -case "$open" > "$work/open-reject.log" 2>&1; then
 echo 'FAIL: M1C18 nonperiodic real native mesh accepted'; exit 51
fi
grep -q 'M1C18_GLOBAL_BLEND_UNSUPPORTED requires cyclic-only boundary' "$work/open-reject.log" || { tail -55 "$work/open-reject.log"; exit 52; }
echo 'PASS: M1C-18 genuine nonperiodic mesh rejected by native cyclic-only gate'
echo 'PASS: M1C-18 six real 3D native Gaussian cases, physically scoped positivity/Bq/Green checks and five negative controls; local FCT/MULES remains NOT IMPLEMENTED'
