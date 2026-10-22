#!/usr/bin/env bash
# Compile and run a REAL OpenFOAM transient FVM analytic diffusion case.
# This is distinct from the M0 compiler/MPI-only smoke.
set -eo pipefail
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -u
work=$(mktemp -d /tmp/cteasfoam-m1a-XXXXXX)
trap 'rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
cp -a /workspace/cases/m1a/diffusion1d "$work/case"
(
  cd "$work/atmosphere"
  wclean
  wmake > "$work/compile.log" 2>&1 || { cat "$work/compile.log"; exit 21; }
)
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/compile.log"; then
  echo "FAIL: wmake emitted errors despite exit zero"
  exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam"
(
  cd "$work/case"
  blockMesh > "$work/blockMesh.log" 2>&1 || { cat "$work/blockMesh.log"; exit 23; }
  "$FOAM_USER_APPBIN/cteasAtmosphereFoam" > "$work/solver.log" 2>&1 || {
    cat "$work/solver.log"; exit 24;
  }
)
python3 /workspace/tests/m1a/check_diffusion.py "$work/case" "$work/solver.log"
echo "PASS: M1A closed-domain finite-volume diffusion verification"
