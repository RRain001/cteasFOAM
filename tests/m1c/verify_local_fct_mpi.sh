#!/usr/bin/env bash
# Sourced by run_local_fct.sh; requires work,exe and tested serial 8_face.
parallel="$work/8_face_mpi"
cp -a "$work/8_face" "$parallel"
cat > "$parallel/system/decomposeParDict" <<'DICT'
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
distributed no;
roots ();
DICT
decomposePar -case "$parallel" -time 0 > "$work/decompose.log" 2>&1 || { tail -65 "$work/decompose.log"; exit 51; }
for rank in 0 1; do test -f "$parallel/processor$rank/constant/polyMesh/owner" || exit 52; done
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 timeout 240s mpirun --oversubscribe -np 2 "$exe" -parallel -case "$parallel" > "$work/mpi.log" 2>&1 || { tail -90 "$work/mpi.log"; exit 53; }
reconstructPar -case "$parallel" -latestTime > "$work/reconstruct.log" 2>&1 || { tail -60 "$work/reconstruct.log"; exit 54; }
python3 - "$work/8_face" "$parallel" "$work/mpi.log" <<'PYMPI'
from pathlib import Path
import sys,re,math
sys.path.insert(0,'/workspace/tests/m1c')
from check_gaussian_green import field,DIM,INV,HTF,ZTF
serial,parallel=map(Path,sys.argv[1:3]); log=Path(sys.argv[3]).read_text()
if len(re.findall(r'(?m)^Time = ',log))!=80 or len(INV.findall(log))!=80 or len(HTF.findall(log))!=80 or len(ZTF.findall(log))!=80:
 raise SystemExit('FAIL: M1C19 native MPI 80-step activity and signed boundary Bq log missing')
for sp in ('HT','HTO'):
 a=field(serial/'0.4'/sp,512,'scalar',DIM)
 b=field(parallel/'0.4'/sp,512,'scalar',DIM)
 delta=max(abs(x-y) for x,y in zip(a,b))
 print(f'M1C19_MPI_PARITY {sp} max_abs_Bq_m3={delta:.12g}')
 if delta>3e-4:raise SystemExit('FAIL: M1C19 native MPI reconstructed species field exceeds fixed 3e-4 Bq/m3 threshold')
serial_inv=INV.findall((serial/'solve.log').read_text())
parallel_inv=INV.findall(log)
for k,(a,b) in enumerate(zip(serial_inv,parallel_inv)):
 if max(abs(float(x)-float(y)) for x,y in zip(a,b))>1e-5:
  raise SystemExit(f'FAIL: M1C19 native MPI signed HT/HTO integrated Bq mismatch at step {k}')
print('PASS: M1C-19 genuine 2-rank OpenFOAM coupled-face/periodic field and 80-step Bq parity')
PYMPI
