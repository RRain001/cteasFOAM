#!/usr/bin/env bash
# M1C-7 end-to-end native v2606 PISO 3D localized Bq/s source.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1C-7 setup_rc=$setup_rc version=${WM_PROJECT_VERSION:-unset}"
[[ "${WM_PROJECT_VERSION:-}" == v2606 || "${WM_PROJECT_VERSION:-}" == 2606 ]] || exit 18
for exe in python3 wmake wclean blockMesh checkMesh postProcess decomposePar mpirun reconstructPar timeout; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c7-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1C-7 rc=$rc work=$work"
 for f in "$work"/build.log "$work"/*.log; do
  if [[ -f "$f" ]]; then echo "== $f =="; tail -90 "$f"; fi
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
serial="$work/serial"
parallel="$work/parallel"
zero="$work/zero"
python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$serial" warped
blockMesh -case "$serial" > "$work/blockMesh.log" 2>&1 || {
 tail -80 "$work/blockMesh.log"; exit 23;
}
python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$serial" warped
checkMesh -case "$serial" -allGeometry -allTopology > "$serial.mesh.log" 2>&1 || {
 tail -80 "$serial.mesh.log"; exit 24;
}
grep -q "Mesh OK." "$serial.mesh.log" || exit 24
postProcess -case "$serial" -func writeCellCentres -time 0 > "$work/centres.log" 2>&1 || exit 25
postProcess -case "$serial" -func writeCellVolumes -time 0 > "$work/volumes.log" 2>&1 || exit 26
python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$serial" warped
python3 /workspace/tests/m1c/make_open_activity.py "$serial" warped
python3 /workspace/tests/m1c/make_open_gaussian_release.py "$serial" release
cp "$serial.mesh.log" "$serial/mesh-check.log"
cp -a "$serial" "$parallel"
cp -a "$serial" "$zero"
# M1C-7 negative control uses precisely the same pressure and mesh with
# zero source Bq/s; keep all other material and FVM inputs identical.
python3 - "$zero/constant/transportProperties" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]);s=p.read_text()
for name,r in (("sourceHT_Bq_s",60000),("sourceHTO_Bq_s",12000)):
 before=f"{name} [0 0 -2 0 0 0 0] {r};"
 if s.count(before)!=1:raise SystemExit(f"FAIL: missing {name} before zero-rate control")
 s=s.replace(before,f"{name} [0 0 -2 0 0 0 0] 0;")
p.write_text(s)
PY
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
timeout 240s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$serial" > "$work/serial.log" 2>&1 || {
 tail -100 "$work/serial.log"; exit 31;
}
timeout 240s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$zero" > "$work/zero.log" 2>&1 || {
 tail -100 "$work/zero.log"; exit 32;
}
decomposePar -case "$parallel" > "$work/decompose.log" 2>&1 || {
 tail -90 "$work/decompose.log"; exit 33;
}
for rank in 0 1; do test -f "$parallel/processor$rank/constant/polyMesh/owner" || exit 34; done
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 240s mpirun --oversubscribe -np 2 \
 "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$parallel" > "$work/parallel.log" 2>&1 || {
 tail -100 "$work/parallel.log"; exit 35;
}
reconstructPar -case "$parallel" -latestTime > "$work/reconstruct.log" 2>&1 || {
 tail -80 "$work/reconstruct.log"; exit 36;
}
for rank in 0 1; do
 for sp in U HT HTO; do test -f "$parallel/processor$rank/0.4/$sp" || exit 37; done
done
for sp in U HT HTO; do test -f "$parallel/0.4/$sp" || exit 38; done
checker=/workspace/tests/m1c/check_open_gaussian_release.py
python3 "$checker" "$serial" "$work/serial.log" release serial "$parallel" "$work/parallel.log" release mpi
python3 "$checker" "$zero" "$work/zero.log" zero_rate serial
python3 - "$serial" "$zero" "$work/serial.log" "$work/zero.log" <<'PYCONTROL'
from pathlib import Path
import sys
sys.path.insert(0,"/workspace/tests/m1c")
from check_open_gaussian_release import check,fail
s,z=map(Path,sys.argv[1:3])
h=check(s,sys.argv[3],"release","serial")
g=check(z,sys.argv[4],"zero_rate","serial")
err=max(abs(x-y) for rowa,rowb in zip(h["U"],g["U"]) for x,y in zip(rowa,rowb))
if err>.003:fail(f"source back-reacted on native pressure-only momentum field {err}")
print(f"PASS: M1C-7 zero-release pressure PISO wind unaffected max_U_delta={err:.12g} m/s")
PYCONTROL
cp "$work/serial.log" "$work/corrupt-source.log"
python3 - "$work/corrupt-source.log" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
s,n=re.subn(r"(emitted_HT_Bq=)([-+0-9.eE]+)",lambda m:m.group(1)+str(float(m.group(2))+1000),s,count=1)
if n!=1:raise SystemExit("FAIL: source log did not contain physical emission")
p.write_text(s)
PY
if python3 "$checker" "$serial" "$work/corrupt-source.log" release serial > "$work/tamper-source.log" 2>&1; then
 echo "FAIL: tampered true emitted HT passed"; exit 40
fi
grep -q 'FAIL: M1C-7' "$work/tamper-source.log" || exit 41
echo "PASS: M1C-7 native emitted Bq log tampering rejected"
python3 - "$serial/0.4/HT" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+960\s*\(\s*)([-+0-9.eE]+)",s)
if not m:raise SystemExit("FAIL: native PISO final HT 960-cell vector missing")
p.write_text(s[:m.start(2)]+str(float(m.group(2))+1000)+s[m.end(2):])
PY
if python3 "$checker" "$serial" "$work/serial.log" release serial > "$work/tamper-field.log" 2>&1; then
 echo "FAIL: corrupted native 3D PISO final HT accepted"; exit 42
fi
grep -q 'FAIL: M1C-7' "$work/tamper-field.log" || exit 43
echo "PASS: M1C-7 960-cell native final HT tampering rejected"
echo "PASS: M1C-7 genuine warped 3-D PISO localized emission with signed HT/HTO flux and MPI"
