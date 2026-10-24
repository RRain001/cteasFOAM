#!/usr/bin/env bash
# Called by run_bounded_diffusion.sh with work/exe/normal native 8³.
for setting in "x:2:2_1_1" "y:2:1_2_1" "xy:4:2_2_1"; do
  IFS=: read -r direction ranks coords <<< "$setting"
  IFS=_ read -r nx ny nz <<< "$coords"
  c="$work/mpi_$direction"
  cp -a "$work/normal" "$c"
  rm -rf "$c/0.4"
  cat > "$c/system/decomposeParDict" <<DICT
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $ranks;
method simple;
simpleCoeffs { n ($nx $ny $nz); delta 0.001; }
distributed no;
roots ();
DICT
  decomposePar -case "$c" -time 0 > "$work/decompose_$direction.log" 2>&1 || { tail -80 "$work/decompose_$direction.log";exit 81; }
  for ((rank=0;rank<ranks;rank++));do test -f "$c/processor$rank/constant/polyMesh/owner" || exit 82;done
  OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 240s \
    mpirun --oversubscribe -np "$ranks" "$exe" -parallel -case "$c" > "$c/solve.log" 2>&1 || { tail -90 "$c/solve.log";exit 83; }
  reconstructPar -case "$c" -latestTime > "$work/reconstruct_$direction.log" 2>&1 || { tail -75 "$work/reconstruct_$direction.log";exit 84; }
  python3 /workspace/tests/m1c/check_robustness.py mpi "$work/normal" "$c" "$c/solve.log" "$ranks:$direction"
  python3 /workspace/tests/m1c/check_bounded_diffusion.py "$c" 8 80 0.005
done
echo "PASS: M1C21 true native x/y 2-rank, xy 4-rank reconstructed HT/HTO and diffusive face budget tests"
