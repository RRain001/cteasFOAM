#!/usr/bin/env bash
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c15-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then echo "FAIL: M1C-15 rc=$rc"; for f in "$work"/build.log "$work"/*/solve.log "$work"/negative-*.log; do if [[ -f "$f" ]]; then echo "== $f =="; tail -35 "$f"; fi; done; fi; rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { tail -70 "$work/build.log"; exit 21; }
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam"
for n in 8 16 24; do
 for method in cellCentre cellTet4; do
  case_dir="$work/$n"_"$method"
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" zero "$n"x"$n"x"$n"
  python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$case_dir" advect
  python3 /workspace/tests/m1c/make_gaussian_green.py "$case_dir" normal
  python3 - "$case_dir" "$method" <<'PYPROP'
from pathlib import Path
import sys
p=Path(sys.argv[1])/"constant/transportProperties"
s=p.read_text()
if sys.argv[2]=="cellTet4":s+="sourceSpatialQuadrature cellTet4;\n"
s+="writeSourceProfile true;\n"
p.write_text(s)
PYPROP
  blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || { tail -55 "$case_dir/blockMesh.log"; exit 23; }
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" zero "$n"x"$n"x"$n"
  checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || { tail -70 "$case_dir/mesh-check.log"; exit 24; }
  grep -q 'Mesh OK.' "$case_dir/mesh-check.log"
  postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$case_dir/centres.log" 2>&1 || exit 25
  postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$case_dir/volumes.log" 2>&1 || exit 26
  timeout 420s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || { tail -80 "$case_dir/solve.log"; exit 31; }
  python3 /workspace/tests/m1c/check_source_quadrature.py --single "$case_dir" "$n" "$method"
 done
done
python3 /workspace/tests/m1c/check_source_quadrature.py --matrix "$work"
# Genuine two-rank source volume quadrature parity on identical native 8^3 mesh.
parallel="$work/parallel8"
cp -a "$work/8_cellTet4" "$parallel"
cat > "$parallel/system/decomposeParDict" <<'DIC'
FoamFile
{
 version 2.0;
 format ascii;
 class dictionary;
 location "system";
 object decomposeParDict;
}
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
distributed no;
roots ();
DIC
decomposePar -case "$parallel" -time 0 > "$parallel/decompose.log" 2>&1 || { tail -65 "$parallel/decompose.log"; exit 32; }
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 300s mpirun --oversubscribe -np 2 "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$parallel" > "$parallel/solve.log" 2>&1 || { tail -80 "$parallel/solve.log";exit 33; }
grep -Eq 'nProcs[[:space:]]*:[[:space:]]*2' "$parallel/solve.log" || exit 34
reconstructPar -case "$parallel" -latestTime > "$parallel/reconstruct.log" 2>&1 || { tail -65 "$parallel/reconstruct.log";exit 35; }
python3 - "$work/8_cellTet4" "$parallel" <<'PYMPI'
import sys,math
from pathlib import Path
sys.path.insert(0,"/workspace/tests/m1c")
from check_gaussian_green import field,INV
p0,p1=map(Path,sys.argv[1:])
for sp in ("HT","HTO"):
    a=field(p0/"0.4"/sp,512,"scalar","0 -3 -1 0 0 0 0")
    b=field(p1/"0.4"/sp,512,"scalar","0 -3 -1 0 0 0 0")
    dif=max(abs(x-y) for x,y in zip(a,b))
    if dif>1e-5:raise SystemExit(f"FAIL: M1C15 native cellTet4 two-rank {sp} differs {dif}")
    print(f"M1C15_MPI_FIELD {sp} max_abs_diff_Bq_m3={dif:.10g}")
logs=[INV.findall((p/"solve.log").read_text()) for p in (p0,p1)]
if not all(len(x)==80 for x in logs):raise SystemExit("FAIL: M1C15 MPI actual HT/HTO 80-step inventory missing")
diff=max(abs(float(a)-float(b)) for l,r in zip(*logs) for a,b in zip(l,r))
if diff>1e-5:raise SystemExit(f"FAIL: M1C15 native 2-rank species Bq parity {diff}")
print(f"PASS: M1C15 real native 2-rank cellTet4 HT/HTO source/activity parity maxBq={diff:.10g}")
PYMPI

native="$work/16_cellTet4"
cp "$native/0/sourceProfile" "$work/pristine-profile"
python3 - "$native/0/sourceProfile" <<'PYFIELD'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+4096\s*\(\s*)([-+0-9.eE]+)",s)
if not m:raise SystemExit("FAIL: native sourceProfile absent")
p.write_text(s[:m.start(2)]+f"{float(m.group(2))+10:.17g}"+s[m.end(2):])
PYFIELD
if python3 /workspace/tests/m1c/check_source_quadrature.py --single "$native" 16 cellTet4 > "$work/negative-profile.log" 2>&1; then echo "FAIL: source kernel corruption accepted";exit 40;fi
grep -q 'FAIL: M1C-15 source profile native volume-weighted mean differs from one' "$work/negative-profile.log" || { cat "$work/negative-profile.log"; exit 41; }
mv "$work/pristine-profile" "$native/0/sourceProfile"
echo "PASS: M1C15 tampered genuine native sourceProfile normalization rejected"
cp "$native/constant/transportProperties" "$work/pristine-props"
sed -i 's/sourceSpatialQuadrature cellTet4;/sourceSpatialQuadrature typo;/' "$native/constant/transportProperties"
if "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$native" > "$work/negative-unknown.log" 2>&1; then echo "FAIL: unknown quadrature accepted";exit 42;fi
grep -q 'unsupported sourceSpatialQuadrature' "$work/negative-unknown.log" || { tail -30 "$work/negative-unknown.log";exit 43; }
mv "$work/pristine-props" "$native/constant/transportProperties"
echo "PASS: M1C15 unknown source quadrature fail-closed"
cp "$native/constant/transportProperties" "$work/pristine-props"
sed -i 's/sourceShape gaussian;/sourceShape uniform;/' "$native/constant/transportProperties"
if "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$native" > "$work/negative-orphan.log" 2>&1; then echo "FAIL: orphan mode accepted";exit 44;fi
grep -q 'require a configured Gaussian activity source' "$work/negative-orphan.log" || { tail -30 "$work/negative-orphan.log";exit 45; }
mv "$work/pristine-props" "$native/constant/transportProperties"
echo "PASS: M1C15 unsupported uniform-source quadrature fail-closed"
cp "$native/solve.log" "$work/pristine-solve"
python3 - "$native/solve.log" <<'PYLOG'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
s,n=re.subn(r"(emitted_HT_Bq=)([-+0-9.eE]+)",lambda m:m.group(1)+str(float(m.group(2))+1000),s,count=1)
if n!=1:raise SystemExit("FAIL: original emitted native Bq missing")
p.write_text(s)
PYLOG
if python3 /workspace/tests/m1c/check_source_quadrature.py --single "$native" 16 cellTet4 > "$work/negative-emitted.log" 2>&1; then echo "FAIL: tampered Bq log accepted";exit 46;fi
grep -q 'FAIL: M1C-14 genuine Gaussian finite-window emission mismatch' "$work/negative-emitted.log" || { cat "$work/negative-emitted.log";exit 47; }
mv "$work/pristine-solve" "$native/solve.log"
echo "PASS: M1C15 signed native emitted HT activity tampering rejected"
echo "PASS: M1C-15 opt-in conservative Gaussian 3D cell tetra quadrature native source, receptors, moments, 6 cases, 4 negative controls"
