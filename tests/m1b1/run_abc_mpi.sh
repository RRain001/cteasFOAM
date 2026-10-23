#!/usr/bin/env bash
# Run exact same nonlinear ABC initial state serial and with 2 MPI ranks
# decomposed along two orthogonal true-3D mesh directions.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B-1d-2 3D MPI setup_rc=$setup_rc OpenFOAM=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for cmd in wmake wclean blockMesh decomposePar reconstructPar mpirun python3 timeout; do
 command -v "$cmd" >/dev/null || { echo "FAIL: missing $cmd"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-abc3dmpi-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "M1B-1d-2 actual 3D MPI FAILED rc=$rc";
 for f in "$work"/compile.log "$work"/mesh*.log "$work"/decompose*.log "$work"/solve*.log "$work"/reconstruct*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -130 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(
 cd "$work/atmosphere"
 wclean
 wmake > "$work/compile.log" 2>&1 || { cat "$work/compile.log"; exit 21; }
)
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/compile.log"; then
 echo "FAIL: hidden solver compiler error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for name in serial x z; do
 python3 /workspace/tests/m1b1/make_abc_case.py "$work/$name" 8
 blockMesh -case "$work/$name" > "$work/mesh$name.log" 2>&1 || {
   cat "$work/mesh$name.log"; exit 23;
 }
done
echo "[M1B-1d-2] real serial ABC 8^3 exact reference run"
timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$work/serial" \
 > "$work/solveserial.log" 2>&1 || { cat "$work/solveserial.log"; exit 24; }

for axis in x z; do
  if [[ "$axis" == x ]]; then parts="2 1 1"; else parts="1 1 2"; fi
  cat > "$work/$axis/system/decomposeParDict" <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 2;
method simple;
simpleCoeffs { n ($parts); delta 0.001; }
distributed no;
roots ();
EOF
  echo "[M1B-1d-2] real 3D 2-rank decomposition along $axis"
  decomposePar -case "$work/$axis" > "$work/decompose$axis.log" 2>&1 || {
    cat "$work/decompose$axis.log"; exit 25;
  }
  for rank in 0 1; do
    test -d "$work/$axis/processor$rank" || {
      echo "FAIL: actual 3D MPI processor$rank missing"; exit 25;
    }
  done
  echo "[M1B-1d-2] actual 3D nonlinear U/p two-rank MPI axis=$axis"
  OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    timeout 180s mpirun --oversubscribe -np 2 \
    "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$work/$axis" \
    > "$work/solve$axis.log" 2>&1 || {
       cat "$work/solve$axis.log"; exit 26;
    }
  reconstructPar -case "$work/$axis" -latestTime \
    > "$work/reconstruct$axis.log" 2>&1 || {
       cat "$work/reconstruct$axis.log"; exit 27;
    }
done
python3 /workspace/tests/m1b1/check_abc_mpi.py \
  "$work/serial" "$work/solveserial.log" \
  "$work/x" "$work/solvex.log" \
  "$work/z" "$work/solvez.log"
echo "PASS: true 3D ABC x/z 2-rank MPI serial equivalence"
