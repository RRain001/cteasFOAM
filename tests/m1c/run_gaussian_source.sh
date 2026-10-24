#!/usr/bin/env bash
# M1C-6 real native Gaussian source, 24/192 cells and two-rank parity.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1C-6 setup_rc=$setup_rc version=${WM_PROJECT_VERSION:-unset}"
[[ "${WM_PROJECT_VERSION:-}" == v2606 || "${WM_PROJECT_VERSION:-}" == 2606 ]] || exit 18
for exe in python3 wmake wclean blockMesh checkMesh postProcess timeout decomposePar mpirun reconstructPar; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c6-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1C-6 rc=$rc work=$work"
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
maker=/workspace/tests/m1c/make_gaussian_source.py
checker=/workspace/tests/m1c/check_gaussian_source.py
for mode in coarse fine parallel; do
 case_dir="$work/$mode"
 if [[ "$mode" == parallel ]]; then size=coarse; else size=$mode; fi
 python3 "$maker" "$case_dir" "$size"
 blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || exit 23
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || {
 tail -80 "$case_dir/mesh-check.log"; exit 24;
 }
 grep -q 'Mesh OK.' "$case_dir/mesh-check.log" || exit 24
 postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$case_dir/centres.log" 2>&1 || exit 25
 postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$case_dir/volumes.log" 2>&1 || exit 26
 if [[ "$mode" == parallel ]]; then
 cat > "$case_dir/system/decomposeParDict" <<'DICT'
FoamFile
{
 version 2.0;
 format ascii;
 class dictionary;
 object decomposeParDict;
}
numberOfSubdomains 2;
method simple;
simpleCoeffs
{
 n (2 1 1);
 delta 0.001;
}
distributed no;
roots ();
DICT
 decomposePar -case "$case_dir" > "$case_dir/decompose.log" 2>&1 || {
 tail -80 "$case_dir/decompose.log"; exit 27;
 }
 for rank in 0 1; do
  test -f "$case_dir/processor$rank/constant/polyMesh/owner" || exit 27
 done
 OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 180s mpirun --oversubscribe -np 2 \
 "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$case_dir" > "$case_dir/solve.log" 2>&1 || {
 tail -90 "$case_dir/solve.log"; exit 31;
 }
 grep -Eq 'nProcs[[:space:]]*:[[:space:]]*2' "$case_dir/solve.log" || exit 32
 reconstructPar -case "$case_dir" -latestTime > "$case_dir/reconstruct.log" 2>&1 || exit 33
 for rank in 0 1; do
  for sp in HT HTO; do test -f "$case_dir/processor$rank/1/$sp" || exit 34; done
 done
 else
 timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || {
 tail -90 "$case_dir/solve.log"; exit 31;
 }
 fi
 python3 "$checker" "$case_dir" "$case_dir/solve.log" "$size"
done
python3 - "$work" <<'PYMPI'
from pathlib import Path
import sys
sys.path.insert(0, "/workspace/tests/m1c")
from check_gaussian_source import check,fail
root=Path(sys.argv[1])
a,ah,az=check(root/"coarse",root/"coarse/solve.log","coarse")
b,bh,bz=check(root/"parallel",root/"parallel/solve.log","coarse")
err=max(abs(x-y) for x,y in zip(ah+az,bh+bz))
inv=max(abs(x-y) for rowa,rowb in zip(a,b) for x,y in zip(rowa,rowb))
if err>1e-6 or inv>1e-6:fail(f"genuine two-rank spatial source mismatch field={err} inventory={inv}")
if (root/"coarse/constant/polyMesh/points").read_bytes() != (root/"parallel/constant/polyMesh/points").read_bytes():
 fail("MPI native geometry was changed")
print(f"PASS: M1C-6 genuine two-rank Gaussian source parity field_max={err:.12g} Bq/m3 inventory_max={inv:.12g} Bq")
PYMPI
python3 - "$work/fine/1/HT" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+192\s*\(\s*)([-+0-9.eE]+)",s)
if not m:raise SystemExit("FAIL: real 192-cell Gaussian HT result missing")
p.write_text(s[:m.start(2)]+str(float(m.group(2))+1000)+s[m.end(2):])
PY
if python3 "$checker" "$work/fine" "$work/fine/solve.log" fine > "$work/negative-field.log" 2>&1; then
 echo "FAIL: corrupted native local HT field accepted"; exit 35
fi
grep -q 'FAIL: M1C-6' "$work/negative-field.log" || exit 36
echo "PASS: M1C-6 tampered native Gaussian HT rejected"
for mode in unknown_shape missing_center negative_sigma bad_center_dims bad_sigma_dims orphan_shape no_support; do
 case_dir="$work/$mode"
 python3 "$maker" "$case_dir" "$mode"
 blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || exit 40
 if timeout 30s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1; then
  echo "FAIL: M1C-6 bad $mode accepted"; exit 41
 fi
 if grep -q '^Time = ' "$case_dir/solve.log"; then
  echo "FAIL: M1C-6 bad $mode emitted after advancing"; exit 42
 fi
 case "$mode" in
  unknown_shape) grep -q 'M1C unsupported sourceShape' "$case_dir/solve.log" || exit 43;;
  missing_center) grep -q 'M1C gaussian source requires sourceCenter and sourceSigma' "$case_dir/solve.log" || exit 43;;
  negative_sigma) grep -q 'M1C gaussian source requires finite center and positive sigma' "$case_dir/solve.log" || exit 43;;
  orphan_shape) grep -q 'M1C sourceShape requires configured Bq/s emission rates' "$case_dir/solve.log" || exit 43;;
  no_support) grep -q 'M1C gaussian source has zero/nonfinite global support' "$case_dir/solve.log" || exit 43;;
  *) grep -q 'FOAM FATAL' "$case_dir/solve.log" || exit 43;;
 esac
 echo "PASS: M1C-6 invalid $mode rejected before first time step"
done
echo "PASS: M1C-6 native normalized Gaussian HT/HTO source and MPI"
