#!/usr/bin/env bash
# Run actual OpenFOAM PDE in serial and two-rank MPI.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1A MPI setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
if [[ "$WM_PROJECT_VERSION" != v2606 && "$WM_PROJECT_VERSION" != 2606 ]]; then
  echo "FAIL: expected v2606"; exit 18
fi
for exe in wmake wclean blockMesh decomposePar reconstructPar mpirun python3 timeout; do
  command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-mpi-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
  echo "M1A two-rank regression failed with rc=$rc";
  for f in "$work"/build.log "$work"/mesh*.log "$work"/decompose.log "$work"/solve*.log "$work"/reconstruct.log; do
    if [[ -f "$f" ]]; then echo "=== $f ==="; tail -100 "$f"; fi
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
  echo "FAIL: hidden compiler diagnostics"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for type in serial parallel; do
  cp -a /workspace/cases/m1a/diffusion1d "$work/$type"
  blockMesh -case "$work/$type" > "$work/mesh$type.log" 2>&1 || {
    cat "$work/mesh$type.log"; exit 23;
  }
done
echo "[M1A MPI] solving serial 40-cell, 200-step diffusion reference"
timeout 120s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$work/serial" \
  > "$work/solveserial.log" 2>&1 || { cat "$work/solveserial.log"; exit 24; }
echo "[M1A MPI] decompose 40-cell mesh along x into 2 domains"
decomposePar -case "$work/parallel" > "$work/decompose.log" 2>&1 || {
  cat "$work/decompose.log"; exit 25;
}
for rank in 0 1; do
  test -d "$work/parallel/processor$rank" || {
    echo "FAIL: actual decomposed processor mesh $rank absent"; exit 25;
  }
done
echo "[M1A MPI] actual 2-rank OpenFOAM PDE execution"
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
timeout 120s mpirun --oversubscribe -np 2 \
"$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$work/parallel" \
> "$work/solveparallel.log" 2>&1 || {
  cat "$work/solveparallel.log"; exit 26;
}
echo "[M1A MPI] reconstruct final field"
reconstructPar -case "$work/parallel" -latestTime \
  > "$work/reconstruct.log" 2>&1 || {
  cat "$work/reconstruct.log"; exit 27;
}
python3 /workspace/tests/m1a/check_parallel_diffusion.py \
  "$work/serial" "$work/parallel" \
  "$work/solveserial.log" "$work/solveparallel.log"
echo "PASS: 2-rank MPI reconstructed PDE equivalence"
