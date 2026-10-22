#!/usr/bin/env bash
# Runs *inside* a candidate OpenFOAM development container.
set -eo pipefail
echo "[cteasFOAM] image runtime check"
if ! command -v wmake >/dev/null 2>&1; then
  for setup in \
    /usr/lib/openfoam/openfoam2606/etc/bashrc \
    /opt/openfoam2606/etc/bashrc \
    /opt/OpenFOAM/OpenFOAM-v2606/etc/bashrc \
    /home/openfoam/OpenFOAM/OpenFOAM-v2606/etc/bashrc
  do
    if [[ -f "$setup" ]]; then
      echo "Source $setup"
      # OpenFOAM env setup may inspect unset variables.
      set +u
      source "$setup" >/dev/null
      break
    fi
  done
fi
if ! command -v wmake >/dev/null 2>&1; then
  echo "FAIL: wmake absent after standard OpenFOAM startup"
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
if [[ -z "${LIB_SRC:-}" && -z "${FOAM_SRC:-}" ]]; then
  echo "WARN: LIB_SRC/FOAM_SRC not defined; check wmake setup"
fi
if [[ ! -f "${FOAM_SRC:-/does-not-exist}/finiteVolume/lnInclude/fvCFD.H" && ! -f "${LIB_SRC:-/does-not-exist}/finiteVolume/lnInclude/fvCFD.H" ]]; then
  echo "FAIL: fvCFD.H not found in installed finiteVolume headers"
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
echo "PASS: OpenFOAM v2606 tools+headers+compiler+MPI+two placeholder builds"
