#!/usr/bin/env bash
# Actual two-rank PISO pressure-velocity solve vs serial TGV.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B1 MPI setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for exe in wmake blockMesh decomposePar reconstructPar mpirun python3 timeout; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-piso-mpi-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "M1B1 MPI FAILED rc=$rc";
 for f in "$work"/build.log "$work"/mesh*.log "$work"/decompose.log "$work"/solve*.log "$work"/reconstruct.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -120 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(
 cd "$work/atmosphere"
 wclean
 wmake > "$work/build.log" 2>&1 || { cat "$work/build.log"; exit 21; }
)
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then
 echo "FAIL: hidden compile error"; exit 22
fi
for name in serial parallel; do
 python3 /workspace/tests/m1b1/make_tgv_case.py "$work/$name" 12
 # Temporal-refinement stress gate: same 0.1s physical interval, twice
 # the number of implicit Euler steps compared to analytic baseline.
 sed -i 's/deltaT 0.005;/deltaT 0.0025;/' "$work/$name/system/controlDict"
 sed -i 's/writeInterval 20;/writeInterval 40;/' "$work/$name/system/controlDict"
 grep -q 'deltaT 0.0025;' "$work/$name/system/controlDict" || exit 28
 grep -q 'writeInterval 40;' "$work/$name/system/controlDict" || exit 28
 blockMesh -case "$work/$name" > "$work/mesh$name.log" 2>&1 || {
  cat "$work/mesh$name.log"; exit 23;
 }
done
echo "[M1B-1b] reference 12x12 pressure-corrected serial solve"
timeout 120s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$work/serial" \
 > "$work/solveserial.log" 2>&1 || { cat "$work/solveserial.log"; exit 24; }
cat > "$work/parallel/system/decomposeParDict" <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
distributed no;
roots ();
EOF
echo "[M1B-1b] split both U and p across two real subdomains"
decomposePar -case "$work/parallel" > "$work/decompose.log" 2>&1 || {
 cat "$work/decompose.log"; exit 25;
}
for rank in 0 1; do
 test -d "$work/parallel/processor$rank" || exit 25
done
echo "[M1B-1b] two-rank coupled U/p PDE"
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
timeout 120s mpirun --oversubscribe -np 2 \
 "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$work/parallel" \
 > "$work/solveparallel.log" 2>&1 || {
 cat "$work/solveparallel.log"; exit 26;
}
echo "[M1B-1b] reconstruct U and p"
reconstructPar -case "$work/parallel" -latestTime \
 > "$work/reconstruct.log" 2>&1 || {
 cat "$work/reconstruct.log"; exit 27;
}
python3 /workspace/tests/m1b1/check_piso_mpi.py \
 "$work/serial" "$work/parallel" \
 "$work/solveserial.log" "$work/solveparallel.log"
echo "PASS: M1B1 serial / actual MPI pressure-correction equivalence"
