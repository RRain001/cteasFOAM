#!/usr/bin/env bash
# Compile and run a REAL OpenFOAM transient FVM analytic diffusion case.
# This is distinct from the M0 compiler/MPI-only smoke.
set -eo pipefail

# OpenFOAM v2606's etc/bashrc is known to return a nonzero status on some
# official APT installations when optional utilities are absent. It still
# successfully exports the required OpenFOAM environment. Do not let errexit
# silently abort before verification; check the actual tools and version.
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1A OpenFOAM setup_rc=$setup_rc WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
for command_name in wmake wclean blockMesh python3; do
    command -v "$command_name" || { echo "FAIL: missing $command_name"; exit 18; }
done
if [[ "${WM_PROJECT_VERSION:-}" != "v2606" && "${WM_PROJECT_VERSION:-}" != "2606" ]]; then
    echo "FAIL: expected OpenFOAM v2606"; exit 19
fi
set -u
work=$(mktemp -d /tmp/cteasfoam-m1a-XXXXXX)
trap 'rc=$?; if (( rc != 0 )); then
        echo "M1A failed with exit $rc; retained diagnostics follow:";
        for f in "$work/compile.log" "$work/blockMesh.log" "$work/solver.log"; do
            if [[ -f "$f" ]]; then echo "=== $f ==="; tail -100 "$f"; fi
        done
    fi
    rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
cp -a /workspace/cases/m1a/diffusion1d "$work/case"
echo "[M1A] compiling real finite-volume scalar solver"
(
  cd "$work/atmosphere"
  wclean
  wmake > "$work/compile.log" 2>&1 || { cat "$work/compile.log"; exit 21; }
)
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/compile.log"; then
  echo "FAIL: wmake emitted errors despite exit zero"
  exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || { echo "FAIL: solver missing"; exit 22; }
echo "[M1A] generating benchmark mesh"
(
  cd "$work/case"
  blockMesh > "$work/blockMesh.log" 2>&1 || { cat "$work/blockMesh.log"; exit 23; }
  echo "[M1A] solving transient finite-volume diffusion"
  "$FOAM_USER_APPBIN/cteasAtmosphereFoam" > "$work/solver.log" 2>&1 || {
    cat "$work/solver.log"; exit 24;
  }
)
echo "[M1A] comparing to analytic Neumann eigenmode and conservation"
python3 /workspace/tests/m1a/check_diffusion.py "$work/case" "$work/solver.log"
echo "PASS: M1A closed-domain finite-volume diffusion verification"
