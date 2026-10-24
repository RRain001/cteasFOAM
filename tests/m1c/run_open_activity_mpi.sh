#!/usr/bin/env bash
# M1C-4: actual warped 3D PISO inlet/outlet HT/HTO under x-axis 2-rank MPI.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1C-4 setup_rc=$setup_rc version=${WM_PROJECT_VERSION:-unset}"
[[ "${WM_PROJECT_VERSION:-}" == v2606 || "${WM_PROJECT_VERSION:-}" == 2606 ]] || exit 18
for exe in python3 wmake wclean blockMesh checkMesh postProcess decomposePar mpirun reconstructPar timeout; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
if [[ -n "${CTEAS_M1C_EVIDENCE_DIR:-}" ]]; then
 mkdir -p "$CTEAS_M1C_EVIDENCE_DIR"
 work=$(mktemp -d "$CTEAS_M1C_EVIDENCE_DIR/cteasfoam-m1c4-XXXXXX")
else
 work=$(mktemp -d /tmp/cteasfoam-m1c4-XXXXXX)
fi
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1C-4 rc=$rc work=$work"
 for f in "$work"/*.log; do [[ -f "$f" ]] && { echo "== $f =="; tail -70 "$f"; }; done
fi
if [[ -n "${CTEAS_M1C_EVIDENCE_DIR:-}" ]]; then
 echo "M1C-4 retained diagnostic directory=$work"
else
 rm -rf "$work"
fi' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { tail -100 "$work/build.log"; exit 21; }
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 21

for mode in warped reverse zero; do
 serial="$work/serial-$mode"
 parallel="$work/parallel-$mode"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$serial" "$mode"
 if [[ "$mode" == warped ]]; then
  blockMesh -case "$serial" > "$work/blockMesh.log" 2>&1 || exit 23
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$serial" "$mode"
  checkMesh -case "$serial" -allGeometry -allTopology > "$work/check.log" 2>&1 || exit 24
  grep -q "Mesh OK." "$work/check.log" || exit 24
  postProcess -case "$serial" -func writeCellCentres -time 0 > "$work/centres.log" 2>&1 || exit 25
  postProcess -case "$serial" -func writeCellVolumes -time 0 > "$work/volumes.log" 2>&1 || exit 26
 else
  cp -a "$work/serial-warped/constant/polyMesh" "$serial/constant/"
  cp -a "$work/serial-warped/0/C" "$serial/0/"
  cp -a "$work/serial-warped/0/V" "$serial/0/"
 fi
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$serial" "$mode"
 python3 /workspace/tests/m1c/make_open_activity.py "$serial" "$mode"
 cp -a "$serial" "$parallel"
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
 timeout 240s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$serial" > "$work/serial-$mode.log" 2>&1 || { tail -90 "$work/serial-$mode.log"; exit 31; }
 decomposePar -case "$parallel" > "$work/decompose-$mode.log" 2>&1 || { tail -90 "$work/decompose-$mode.log"; exit 32; }
 for rank in 0 1; do test -f "$parallel/processor$rank/constant/polyMesh/owner" || exit 33; done
 OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 240s mpirun --oversubscribe -np 2 "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$parallel" > "$work/parallel-$mode.log" 2>&1 || { tail -90 "$work/parallel-$mode.log"; exit 34; }
 reconstructPar -case "$parallel" -latestTime > "$work/reconstruct-$mode.log" 2>&1 || { tail -90 "$work/reconstruct-$mode.log"; exit 35; }
 echo "M1C-4 completed mode=$mode genuine_serial_MPI=2"
done
checker=/workspace/tests/m1c/check_open_activity.py
python3 "$checker" "$work/check.log" "$work/serial-warped" "$work/serial-warped.log" "$work/serial-reverse" "$work/serial-reverse.log" "$work/serial-zero" "$work/serial-zero.log"
python3 "$checker" "$work/check.log" "$work/parallel-warped" "$work/parallel-warped.log" "$work/parallel-reverse" "$work/parallel-reverse.log" "$work/parallel-zero" "$work/parallel-zero.log"
mpi_checker=/workspace/tests/m1c/check_open_activity_mpi.py
python3 "$mpi_checker" "$work"
field="$work/parallel-warped/0.4/HT"
cp "$field" "$work/pristine-HT"
python3 - "$field" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]); s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+960\s*\(\s*)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: missing reconstructed native HT for negative control")
s=s[:m.start(2)]+str(float(m.group(2))+100000)+s[m.end(2):]
p.write_text(s)
PY
if python3 "$mpi_checker" "$work" > "$work/tamper.log" 2>&1; then
 echo "FAIL: M1C-4 corrupted parallel native HT field accepted"; exit 36
fi
grep -q "FAIL: M1C-4" "$work/tamper.log" || { cat "$work/tamper.log"; exit 37; }
echo "PASS: M1C-4 deliberately corrupted MPI native HT rejected"
echo "PASS: M1C-4 official v2606 warped PISO HT/HTO open-boundary two-rank MPI"
