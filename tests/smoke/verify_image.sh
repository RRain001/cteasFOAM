#!/usr/bin/env bash
# Runs inside the official-package OpenFOAM v2606 development container.
set -eo pipefail
echo "[cteasFOAM] official v2606 development-image check"

setup=/usr/lib/openfoam/openfoam2606/etc/bashrc
if [[ ! -f "$setup" ]]; then
  echo "FAIL: official OpenFOAM v2606 environment script absent: $setup"
  exit 9
fi

# Explicitly tolerate an optional environment-shell return status while
# separately verifying every build requirement below.
set +e
source "$setup"
setup_rc=$?
set -e
echo "OpenFOAM setup status: $setup_rc"

if ! command -v wmake >/dev/null 2>&1; then
  echo "FAIL: wmake absent"
  exit 10
fi
for tool in g++ mpirun; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "FAIL: $tool absent"
    exit 11
  fi
done
echo "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-<unset>}"
echo "WM_PROJECT_DIR=${WM_PROJECT_DIR:-<unset>}"
echo "wmake=$(command -v wmake)"
echo "g++=$(g++ --version | head -n1)"
echo "mpirun=$(mpirun --version | head -n1)"
if [[ "${WM_PROJECT_VERSION:-}" != "v2606" && "${WM_PROJECT_VERSION:-}" != "2606" ]]; then
  echo "FAIL: not OpenFOAM v2606"
  exit 12
fi
if [[ ! -f "${FOAM_SRC:-/absent}/finiteVolume/lnInclude/fvCFD.H" && ! -f "${LIB_SRC:-/absent}/finiteVolume/lnInclude/fvCFD.H" ]]; then
  echo "FAIL: fvCFD.H is absent"
  exit 13
fi

export FOAM_USER_APPBIN=/tmp/cteasfoam-m0-bin
mkdir -p "$FOAM_USER_APPBIN" /tmp/cteasfoam-m0-build
for component in atmosphere soil; do
  echo "[cteasFOAM] compile $component"
  cp -a "/workspace/src/$component" "/tmp/cteasfoam-m0-build/$component"
  (
    cd "/tmp/cteasfoam-m0-build/$component"
    wmake
  )
done
"$FOAM_USER_APPBIN/cteasAtmosphereFoam"
"$FOAM_USER_APPBIN/cteasSoilFoam"
echo "PASS: official-package v2606 tools, headers, compiler, MPI and both placeholder builds"
