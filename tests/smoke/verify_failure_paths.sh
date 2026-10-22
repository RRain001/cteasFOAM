#!/usr/bin/env bash
# Run inside the development image; reject the two observed false-success paths.
set -euo pipefail
expect_rejected() {
  local name=$1 expected=$2 injection=$3 rc
  echo "[cteasFOAM] negative regression: $name"
  set +e
  bash -c "$injection; bash /workspace/tests/smoke/verify_image.sh"
  rc=$?
  set -e
  echo "CASE=$name EXPECTED=$expected OBSERVED=$rc"
  if [[ $rc != "$expected" ]]; then
    echo "FAIL: injected toolchain failure was not rejected correctly"
    exit 1
  fi
}
expect_rejected broken-mpi-version 42 \
  'mpirun() { echo "INJECTED: MPI failure" >&2; return 42; }; export -f mpirun'
expect_rejected hidden-wmake-error 14 \
  'wmake() { echo "Make/options:1: error: injected preprocessing failure" >&2; return 0; }; export -f wmake'
echo "PASS: both false-success regression checks rejected the injected failure"
