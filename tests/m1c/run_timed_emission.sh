#!/usr/bin/env bash
# M1C-5 official v2606 native finite-duration HT/HTO Bq/s release.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1C-5 setup_rc=$setup_rc version=${WM_PROJECT_VERSION:-unset}"
[[ "${WM_PROJECT_VERSION:-}" == v2606 || "${WM_PROJECT_VERSION:-}" == 2606 ]] || exit 18
for exe in python3 wmake wclean blockMesh checkMesh postProcess timeout; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c5-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1C-5 rc=$rc work=$work"
 for f in "$work"/build.log "$work"/*/solve.log "$work"/negative-*.log; do
  if [[ -f "$f" ]]; then echo "== $f =="; tail -65 "$f"; fi
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
maker=/workspace/tests/m1c/make_timed_emission.py
checker=/workspace/tests/m1c/check_timed_emission.py
for mode in pulse ht_only zero_rate; do
 case_dir="$work/$mode"
 python3 "$maker" "$case_dir" "$mode"
 blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || exit 23
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || exit 24
 grep -q 'Mesh OK.' "$case_dir/mesh-check.log" || exit 24
 postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$case_dir/volumes.log" 2>&1 || {
 tail -80 "$case_dir/volumes.log"; exit 25;
 }
 timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || {
 tail -90 "$case_dir/solve.log"; exit 31;
 }
 python3 "$checker" "$case_dir" "$case_dir/solve.log" "$mode"
done
# Genuine two-rank source-volume reduction and time integration check.
for exe in decomposePar mpirun reconstructPar; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing MPI $exe"; exit 43; }
done
parallel="$work/pulse-parallel"
python3 "$maker" "$parallel" pulse
blockMesh -case "$parallel" > "$parallel/blockMesh.log" 2>&1 || exit 44
checkMesh -case "$parallel" -allGeometry -allTopology > "$parallel/mesh-check.log" 2>&1 || exit 45
grep -q 'Mesh OK.' "$parallel/mesh-check.log" || exit 45
postProcess -case "$parallel" -func writeCellVolumes -time 0 > "$parallel/volumes.log" 2>&1 || exit 46
cat > "$parallel/system/decomposeParDict" <<'DICT'
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
decomposePar -case "$parallel" > "$parallel/decompose.log" 2>&1 || {
 tail -90 "$parallel/decompose.log"; exit 47;
}
for rank in 0 1; do
 test -f "$parallel/processor$rank/constant/polyMesh/owner" || exit 48
done
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 180s mpirun --oversubscribe -np 2 \
 "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$parallel" > "$parallel/solve.log" 2>&1 || {
 tail -90 "$parallel/solve.log"; exit 49;
}
reconstructPar -case "$parallel" -latestTime > "$parallel/reconstruct.log" 2>&1 || {
 tail -90 "$parallel/reconstruct.log"; exit 50;
}
for name in HT HTO; do
 test -f "$parallel/1/$name" || exit 51
 for rank in 0 1; do
  test -f "$parallel/processor$rank/1/$name" || exit 52
 done
done
grep -Eq 'nProcs[[:space:]]*:[[:space:]]*2' "$parallel/solve.log" || {
 echo "FAIL: not actual two-rank MPI source case"; exit 53;
}
python3 "$checker" "$parallel" "$parallel/solve.log" pulse
python3 - "$work/pulse" "$parallel" <<'PYMPI'
from pathlib import Path
import re,math,sys
a,b=[Path(p) for p in sys.argv[1:]]
num=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
def field(path):
 s=path.read_text()
 u=re.search(r"internalField\s+uniform\s+("+num+r")\s*;",s)
 if u: return [float(u.group(1))]*24
 v=re.search(r"internalField\s+nonuniform\s+List<scalar>\s+24\s*\((.*?)\)\s*;",s,re.S)
 if not v: raise SystemExit("FAIL: MPI native final scalar field malformed")
 result=[float(t) for t in v.group(1).split()]
 if len(result)!=24: raise SystemExit("FAIL: MPI native wrong field count")
 return result
def inventory(path):
 s=path.read_text()
 pat=r"CTEAS_M1C activity_HT_Bq=("+num+r") activity_HTO_Bq=("+num+r")"
 vals=[tuple(map(float,v)) for v in re.findall(pat,s)]
 if len(vals)!=20: raise SystemExit("FAIL: MPI real 20-step history absent")
 return vals
maxfield=max(abs(x-y) for sp in ("HT","HTO") for x,y in zip(field(a/"1"/sp),field(b/"1"/sp)))
maxinv=max(abs(x-y) for xrow,yrow in zip(inventory(a/"solve.log"),inventory(b/"solve.log")) for x,y in zip(xrow,yrow))
if maxfield>1e-6 or maxinv>1e-6:
 raise SystemExit(f"FAIL: MPI source parity field={maxfield} inventory={maxinv}")
print(f"PASS: M1C-5 genuine two-rank source MPI final field maxDelta={maxfield:.12g} Bq/m3 20-step inventory maxDelta={maxinv:.12g} Bq")
PYMPI
cp "$work/pulse/solve.log" "$work/negative-emission.log"
python3 - "$work/negative-emission.log" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]); s=p.read_text()
s,n=re.subn(r"(emitted_HT_Bq=)([-+0-9.eE]+)",lambda m:m.group(1)+str(float(m.group(2))+1000),s,count=1)
if n!=1: raise SystemExit("FAIL: missing native emission log negative control")
p.write_text(s)
PY
if python3 "$checker" "$work/pulse" "$work/negative-emission.log" pulse > "$work/negative-log-result.log" 2>&1; then
 echo "FAIL: M1C-5 tampered emission log accepted"; exit 35
fi
grep -q 'FAIL: M1C-5' "$work/negative-log-result.log" || exit 36
echo "PASS: M1C-5 changed emitted activity log rejected"
field="$work/pulse/1/HT"
cp "$field" "$work/pristine-HT"
python3 - "$field" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]); s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+24\s*\(\s*|internalField\s+uniform\s+)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: no real native final HT field")
p.write_text(s[:m.start(2)]+str(float(m.group(2))+1000)+s[m.end(2):])
PY
if python3 "$checker" "$work/pulse" "$work/pulse/solve.log" pulse > "$work/negative-field-result.log" 2>&1; then
 echo "FAIL: M1C-5 tampered HT native field accepted"; exit 37
fi
grep -q 'FAIL: M1C-5' "$work/negative-field-result.log" || exit 38
echo "PASS: M1C-5 changed native final HT field rejected"
for mode in negative_rate reverse_window missing_window wrong_dimensions; do
 case_dir="$work/$mode"
 python3 "$maker" "$case_dir" "$mode"
 blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || exit 23
 if timeout 30s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1; then
  echo "FAIL: invalid emission case $mode accepted"; exit 39
 fi
 case "$mode" in
  missing_window) grep -q 'requires sourceStart and sourceEnd' "$case_dir/solve.log" || { tail -40 "$case_dir/solve.log"; exit 40; };;
  wrong_dimensions) grep -Eq 'FOAM FATAL|incompatible dimensions' "$case_dir/solve.log" || { tail -40 "$case_dir/solve.log"; exit 41; };;
  *) grep -q 'requires nonnegative finite Bq/s' "$case_dir/solve.log" || { tail -40 "$case_dir/solve.log"; exit 42; };;
 esac
 echo "PASS: M1C-5 invalid $mode explicitly rejected"
done
echo "PASS: M1C-5 native OpenFOAM finite-duration Bq/s HT/HTO emission"
