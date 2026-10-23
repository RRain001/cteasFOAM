#!/usr/bin/env bash
# M1B-1g-1: real three-dimensional non-periodic, interior-warped
# pressure-duct MPI comparison; serial / x split / z split / zero dp.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B-1g-1 official OpenFOAM setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for exe in wmake wclean blockMesh checkMesh postProcess decomposePar reconstructPar mpirun python3 timeout cmp; do
 command -v "$exe" >/dev/null || { echo "FAIL: missing $exe"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-duct-mpi-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1B-1g-1 warped-pressure MPI rc=$rc";
 for f in "$work"/build.log "$work"/mesh*.log "$work"/check*.log "$work"/centres*.log "$work"/decompose*.log "$work"/solve*.log "$work"/reconstruct*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -100 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { cat "$work/build.log"; exit 21; }
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then
 echo "FAIL: hidden compiler error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for tag in serial zero-z; do
 if [[ "$tag" == serial ]]; then mode=warped; else mode=zero; fi
 case_dir="$work/$tag"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" "$mode"
 blockMesh -case "$case_dir" > "$work/mesh$tag.log" 2>&1 || { tail -120 "$work/mesh$tag.log"; exit 23; }
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" "$mode" || exit 24
 checkMesh -case "$case_dir" -allGeometry -allTopology > "$work/check$tag.log" 2>&1 || {
  tail -120 "$work/check$tag.log"; exit 25;
 }
 grep -q 'Mesh OK.' "$work/check$tag.log" || exit 25
 postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$work/centres$tag.log" 2>&1 || {
  tail -120 "$work/centres$tag.log"; exit 26;
 }
 test -s "$case_dir/0/C" || exit 26
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$case_dir" "$mode" || exit 27
done
# Copy the exact OpenFOAM physical mesh and actual C/initial values instead
# of regenerating mathematically similar-but-different MPI test cases.
for tag in x z; do
 cp -a "$work/serial" "$work/$tag"
 for f in constant/polyMesh/points 0/C 0/U 0/p 0/c; do
  cmp -s "$work/serial/$f" "$work/$tag/$f" || {
   echo "FAIL: nonidentical x/z MPI physical input $tag/$f"; exit 28;
  }
 done
done
echo 'PASS: physically identical warped OpenFOAM points, C, U, p, c before MPI'
timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$work/serial" \
 > "$work/solveserial.log" 2>&1 || { tail -120 "$work/solveserial.log"; exit 31; }
for tag in x z zero-z; do
 if [[ "$tag" == x ]]; then parts='2 1 1'; else parts='1 1 2'; fi
 cat > "$work/$tag/system/decomposeParDict" <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 2;
method simple;
simpleCoeffs { n ($parts); delta 0.001; }
distributed no;
roots ();
EOF
 decomposePar -case "$work/$tag" > "$work/decompose$tag.log" 2>&1 || {
  tail -120 "$work/decompose$tag.log"; exit 32;
 }
 for rank in 0 1; do test -d "$work/$tag/processor$rank" || exit 32; done
 OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  timeout 180s mpirun --oversubscribe -np 2 \
  "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -parallel -case "$work/$tag" \
  > "$work/solve$tag.log" 2>&1 || { tail -150 "$work/solve$tag.log"; exit 33; }
 reconstructPar -case "$work/$tag" -latestTime > "$work/reconstruct$tag.log" 2>&1 || {
  tail -120 "$work/reconstruct$tag.log"; exit 34;
 }
done
for tag in serial x z zero-z; do
 steps=$(grep -c '^CTEAS_M1B1 nonorthogonal_correctors=2$' "$work/solve$tag.log" || true)
 solves=$(grep -c 'Solving for p,' "$work/solve$tag.log" || true)
 echo "M1B-1g-1 tag=$tag steps=$steps pressure_solves=$solves"
 [[ "$steps" -eq 20 && "$solves" -eq 180 ]] || exit 35
done
python3 /workspace/tests/m1b1/check_warped_open_pressure_mpi.py \
 "$work/serial" "$work/checkserial.log" "$work/solveserial.log" \
 "$work/x" "$work/solvex.log" "$work/z" "$work/solvez.log" \
 "$work/zero-z" "$work/checkzero-z.log" "$work/solvezero-z.log"
echo 'PASS: M1B-1g-1 official-v2606 full three-dimensional open-pressure MPI regression'
