#!/usr/bin/env bash
# Runs inside the official-package OpenFOAM v2606 development container.
set -eo pipefail
echo "[cteasFOAM] official v2606 development-image check"
source /etc/os-release
[[ "$ID" == ubuntu && "$VERSION_ID" == 22.04 && $(uname -m) == x86_64 ]]
dpkg-query -W -f='${Package}=${Version}\n' \
  openfoam2606 openfoam2606-dev openfoam2606-default g++ openmpi-bin libopenmpi-dev

setup=/usr/lib/openfoam/openfoam2606/etc/bashrc
if [[ ! -f "$setup" ]]; then
  echo "FAIL: official OpenFOAM v2606 environment script absent: $setup"
  exit 9
fi
# Upstream setup is not nounset-safe and can return an optional-tool status.
set +e
source "$setup"
setup_rc=$?
set -e
echo "OpenFOAM setup status: $setup_rc"
set -u
for tool in wmake wclean g++ mpirun mpicxx timeout; do
  command -v "$tool"
done
echo "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-<unset>}"
echo "WM_PROJECT_DIR=${WM_PROJECT_DIR:-<unset>}"
echo "WM_OPTIONS=${WM_OPTIONS:-<unset>} WM_MPLIB=${WM_MPLIB:-<unset>}"
# Invoke directly: command substitutions inside echo can hide tool failure.
g++ --version
mpirun --version
mpicxx --showme
if [[ "${WM_PROJECT_VERSION:-}" != "v2606" && "${WM_PROJECT_VERSION:-}" != "2606" ]]; then
  echo "FAIL: not OpenFOAM v2606"
  exit 12
fi
for header in finiteVolume/lnInclude/fvCFD.H meshTools/lnInclude/cyclicAMIPolyPatch.H; do
  if [[ ! -f "$WM_PROJECT_DIR/src/$header" ]]; then
    echo "FAIL: header absent: $header"
    exit 13
  fi
done

build_dir=$(mktemp -d /tmp/cteasfoam-m0-XXXXXX)
export FOAM_USER_APPBIN="$build_dir/bin"
mkdir -p "$FOAM_USER_APPBIN"
for component in atmosphere soil; do
  echo "[cteasFOAM] compile $component"
  cp -a "/workspace/src/$component" "$build_dir/$component"
  (
    cd "$build_dir/$component"
    # Ignore any copied build artefacts: this probe must compile from source.
    wclean
    compile_log="$build_dir/$component-wmake.txt"
    wmake 2>&1 | tee "$compile_log"
    # wmake can return zero even when preprocessing Make/options reports an
    # error. Reject diagnostics as well as the pipeline's nonzero exit code.
    if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$compile_log"; then
      echo "FAIL: $component emitted build errors despite wmake exit status"
      exit 14
    fi
  )
done
for app in cteasAtmosphereFoam cteasSoilFoam; do
  ldd "$FOAM_USER_APPBIN/$app"
  "$FOAM_USER_APPBIN/$app"
done

echo "[cteasFOAM] independent MPI compiler/runtime smoke (two ranks)"
mpicxx /workspace/tests/smoke/mpi_probe.C -o "$build_dir/mpi_probe"
# The disposable test container runs as root, so scope OpenMPI consent here.
OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  timeout 60s mpirun --oversubscribe -np 2 "$build_dir/mpi_probe"
echo "PASS: official-package v2606 tools, headers, two placeholder builds and MPI toolchain"
echo "NOT_RUN: solver parallel numerics, PDE physics and coupling conservation (NOT IMPLEMENTED)"
