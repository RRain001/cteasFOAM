#!/usr/bin/env bash
# M1C-2: genuine 2-rank OpenFOAM v2606 HT/HTO transport and budget parity.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1C-2 setup_rc=$setup_rc version=${WM_PROJECT_VERSION:-unset}"
[[ "${WM_PROJECT_VERSION:-}" == v2606 || "${WM_PROJECT_VERSION:-}" == 2606 ]] || exit 18
for app in python3 blockMesh checkMesh postProcess wmake wclean decomposePar reconstructPar mpirun timeout; do
 command -v "$app" >/dev/null || { echo "FAIL: missing $app"; exit 18; }
done
set -u
# The frozen M1C-2 admission case uses x decomposition by default.
# Optional y/z are diagnostic controls: they do NOT replace x admission.
axis=${CTEAS_M1C_DECOMPOSITION:-x}
case "$axis" in
 x) partitions='2 1 1' ;;
 y) partitions='1 2 1' ;;
 z) partitions='1 1 2' ;;
 *) echo "FAIL: unsupported decomposition $axis (x/y/z)"; exit 18 ;;
esac
echo "M1C-2 decomposition=$axis partitions=($partitions)"
work=$(mktemp -d /tmp/cteasfoam-m1c-mpi-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1C-2 rc=$rc work=$work"
 for f in "$work"/build.log "$work"/*.decompose.log "$work"/*.mpi.log "$work"/*.reconstruct.log "$work"/*.negative.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -90 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { tail -100 "$work/build.log"; exit 21; }
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then exit 22; fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for entry in coarse:8x20x6:960 fine:12x30x9:3240; do
 IFS=: read -r label resolution count <<< "$entry"
 serial="$work/serial-$label"
 parallel="$work/parallel-$label"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$serial" zero "$resolution"
 python3 /workspace/tests/m1b1/make_warped_periodic_advection_diffusion.py prepare "$serial"
 python3 /workspace/tests/m1c/make_spatial_activity.py prepare "$serial"
 # Tighten the iterative species solve (both serial and MPI) without
 # changing the predeclared parity, physics or discretization thresholds.
 python3 - "$serial/system/fvSolution" <<'PYTOL'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]); s=p.read_text()
for name in ('HT','HTO'):
    pattern=rf'({name} \{{[^}}]*tolerance )1e-12;'
    s,n=re.subn(pattern,r'\g<1>1e-15;',s)
    if n!=1: raise SystemExit(f'FAIL: missing {name} solver tolerance')
p.write_text(s)
PYTOL
 blockMesh -case "$serial" > "$work/$label.mesh.log" 2>&1 || { tail -90 "$work/$label.mesh.log"; exit 23; }
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$serial" zero "$resolution"
 checkMesh -case "$serial" -allGeometry -allTopology > "$work/$label.check.log" 2>&1 || { tail -90 "$work/$label.check.log"; exit 24; }
 grep -q 'Mesh OK.' "$work/$label.check.log" || exit 24
 postProcess -case "$serial" -func writeCellCentres -time 0 > "$work/$label.centres.log" 2>&1 || exit 25
 postProcess -case "$serial" -func writeCellVolumes -time 0 > "$work/$label.volumes.log" 2>&1 || exit 26
 python3 /workspace/tests/m1c/make_spatial_activity.py seed "$serial" "$count"
 cp -a "$serial" "$parallel"
 cat > "$parallel/system/decomposeParDict" <<DICT
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
 n ($partitions);
 delta 0.001;
}
distributed no;
roots ();
DICT
 timeout 240s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$serial" > "$work/$label.serial.log" 2>&1 || { tail -90 "$work/$label.serial.log"; exit 31; }
 decomposePar -case "$parallel" > "$work/$label.decompose.log" 2>&1 || { tail -90 "$work/$label.decompose.log"; exit 32; }
 for rank in 0 1; do
  test -f "$parallel/processor$rank/constant/polyMesh/owner" || { echo "FAIL: missing actual processor mesh $rank"; exit 32; }
 done
 OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 240s mpirun --oversubscribe -np 2 \
  "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$parallel" > "$work/$label.mpi.log" 2>&1 || { tail -100 "$work/$label.mpi.log"; exit 33; }
 reconstructPar -case "$parallel" -latestTime > "$work/$label.reconstruct.log" 2>&1 || { tail -90 "$work/$label.reconstruct.log"; exit 34; }
 for species in HT HTO; do
  test -f "$parallel/0.4/$species" || { echo "FAIL: missing reconstructed $species"; exit 34; }
 done
done
# The original M1C-1 checker must independently accept both real serial and
# reconstructed MPI cases against continuum waves and FVM species budgets.
base_checker=/workspace/tests/m1c/check_spatial_activity.py
python3 "$base_checker" "$work/serial-coarse" "$work/coarse.check.log" "$work/coarse.serial.log" \
 "$work/serial-fine" "$work/fine.check.log" "$work/fine.serial.log"
python3 "$base_checker" "$work/parallel-coarse" "$work/coarse.check.log" "$work/coarse.mpi.log" \
 "$work/parallel-fine" "$work/fine.check.log" "$work/fine.mpi.log"
checker=/workspace/tests/m1c/check_mpi_activity.py
args=("$work/serial-coarse" "$work/parallel-coarse" "$work/coarse.serial.log" "$work/coarse.mpi.log" \
 "$work/serial-fine" "$work/parallel-fine" "$work/fine.serial.log" "$work/fine.mpi.log")
python3 "$checker" "${args[@]}"
python3 - "$work/parallel-fine/0.4/HT" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+3240\s*\(\s*)([-+0-9.eE]+)",s)
if not m: raise SystemExit("FAIL: reconstructed native HT missing")
s=s[:m.start(2)]+f"{float(m.group(2))+10:.17g}"+s[m.end(2):]
p.write_text(s)
PY
if python3 "$checker" "${args[@]}" > "$work/tamper.negative.log" 2>&1; then
 echo "FAIL: corrupted reconstructed HT passed"; exit 35
fi
grep -q 'FAIL: M1C-2' "$work/tamper.negative.log" || { cat "$work/tamper.negative.log"; exit 36; }
echo 'PASS: deliberately corrupted reconstructed MPI HT field rejected'
if [[ "$axis" == x ]]; then
 echo 'PASS: M1C-2 frozen x-decomposition official v2606 two-rank species parity'
else
 echo "PASS: M1C-2 $axis-decomposition diagnostic-only parity (frozen x gate NOT_RUN in this invocation)"
fi
