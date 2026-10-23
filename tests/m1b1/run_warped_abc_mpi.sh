#!/usr/bin/env bash
# Real warped 8^3 ABC with two extra nonorth pressure corrections,
# compare the same exact polyMesh, centroids and initial U/p across 3 jobs.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
rc=$?
set -e
echo "M1B-1f3 setup_rc=$rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for cmd in wmake wclean blockMesh checkMesh postProcess decomposePar reconstructPar mpirun python3 timeout cmp; do
  command -v "$cmd" >/dev/null || { echo "FAIL: missing $cmd"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-warpedMPI-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: warped 3D MPI scientific regression rc=$rc";
 for f in "$work"/compile.log "$work"/mesh*.log "$work"/check*.log "$work"/centres*.log "$work"/decompose*.log "$work"/solve*.log "$work"/reconstruct*.log; do
   if [[ -f "$f" ]]; then echo "=== $f ==="; tail -140 "$f"; fi
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
  echo "FAIL: hidden wmake error"; exit 22;
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22

for name in serial x z; do
 case_dir="$work/$name"
 python3 /workspace/tests/m1b1/make_abc_case.py "$case_dir" 8
 blockMesh -case "$case_dir" > "$work/mesh$name.log" 2>&1 || {
   cat "$work/mesh$name.log"; exit 23;
 }
 python3 /workspace/tests/m1b1/warp_abc_mesh.py warp "$case_dir" 8 0.06 || exit 24
 checkMesh -case "$case_dir" -allGeometry -allTopology \
   > "$work/check$name.log" 2>&1 || { tail -120 "$work/check$name.log"; exit 25; }
 grep -q "Mesh OK." "$work/check$name.log" || exit 25
 postProcess -case "$case_dir" -func writeCellCentres -time 0 \
   > "$work/centres$name.log" 2>&1 || { tail -120 "$work/centres$name.log"; exit 26; }
 test -s "$case_dir/0/C" || exit 26
 python3 /workspace/tests/m1b1/warp_abc_mesh.py seed "$case_dir" 8 || exit 27
 sed -i '/^nPressureCorrectors 3;$/a nNonOrthogonalCorrectors 2;' \
   "$case_dir/constant/transportProperties"
 grep -q '^nNonOrthogonalCorrectors 2;$' \
   "$case_dir/constant/transportProperties" || exit 27
done
for f in constant/polyMesh/points 0/C 0/U 0/p; do
  cmp -s "$work/serial/$f" "$work/x/$f" || {
    echo "FAIL: nonidentical x-split physical input $f"; exit 28;
  }
  cmp -s "$work/serial/$f" "$work/z/$f" || {
    echo "FAIL: nonidentical z-split physical input $f"; exit 28;
  }
done
echo "PASS: identical warped polyMesh, OpenFOAM centroids and initial ABC fields"

echo "[M1B-1f3] serial warped 8^3 nonorth pressure solve"
timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$work/serial" \
  > "$work/solveserial.log" 2>&1 || { tail -150 "$work/solveserial.log"; exit 31; }

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
  decomposePar -case "$work/$axis" > "$work/decompose$axis.log" 2>&1 || {
    tail -150 "$work/decompose$axis.log"; exit 32;
  }
  for rank in 0 1; do
    test -d "$work/$axis/processor$rank" || exit 32
  done
  echo "[M1B-1f3] genuine deformed 3D MPI solver axis=$axis"
  OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    timeout 180s mpirun --oversubscribe -np 2 \
    "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$work/$axis" \
    > "$work/solve$axis.log" 2>&1 || { tail -150 "$work/solve$axis.log"; exit 33; }
  reconstructPar -case "$work/$axis" -latestTime \
    > "$work/reconstruct$axis.log" 2>&1 || {
      tail -150 "$work/reconstruct$axis.log"; exit 34;
    }
done

for name in serial x z; do
  steps=$(grep -c '^CTEAS_M1B1 nonorthogonal_correctors=2$' "$work/solve$name.log" || true)
  psolves=$(grep -c 'Solving for p,' "$work/solve$name.log" || true)
  echo "M1B-1f3 nonorthMPI name=$name steps=$steps pressure_solves=$psolves"
  [[ "$steps" -eq 20 && "$psolves" -eq 180 ]] || {
    echo "FAIL: real inner nonorth pressure solves count mismatch"; exit 35;
  }
done
python3 /workspace/tests/m1b1/check_warped_abc_mpi.py \
  "$work/serial" "$work/solveserial.log" \
  "$work/x" "$work/solvex.log" \
  "$work/z" "$work/solvez.log"
echo "PASS: warped 3D serial/x-MPI/z-MPI true pressure-velocity PDE agreement"
