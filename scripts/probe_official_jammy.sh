#!/usr/bin/env bash
# Actual build-and-compile acceptance for official OpenCFD v2606 packages.
# Run on a native x86-64 Linux Docker host; Mac is for editing/static checks.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

mkdir -p reports
# Keep every attempt, including failed builds, rather than truncating evidence.
run_dir=$(mktemp -d "$PWD/reports/official-jammy-$(date -u +%Y%m%dT%H%M%SZ)-XXXXXX")
report="$run_dir/probe.txt"
exec > >(tee "$report") 2>&1
trap 'rc=$?; echo "PROBE_EXIT_STATUS=$rc"' EXIT
image="cteasfoam/official-v2606-dev:jammy-m0"

echo "UTC=$(date -u +%FT%TZ)"
echo "SOURCE_HEAD=$(git rev-parse HEAD)"
git status --short
git diff --binary HEAD > "$run_dir/source.patch"
git ls-files --cached --others --exclude-standard -z | xargs -0 shasum -a 256 > "$run_dir/source-sha256.txt"
echo "BASE: ubuntu:22.04 / official dl.openfoam.com v2606 packages"
if [[ $(uname -s) != Linux || $(uname -m) != x86_64 ]]; then
  echo "NOT_RUN: native x86-64 Linux host required; no Mac/emulated acceptance"
  exit 2
fi
cat /etc/os-release
docker info --format 'DOCKER={{.ServerVersion}} OS={{.OSType}} ARCH={{.Architecture}}'
case "$(docker info --format '{{.OSType}}/{{.Architecture}}')" in
  linux/x86_64|linux/amd64) ;;
  *) echo "NOT_RUN: Docker daemon must be native Linux x86-64"; exit 2 ;;
esac

build_args=()
# Optional host-local proxy for package downloads. Docker's daemon pull proxy
# is configured separately by the operator, never embedded in this image.
if [[ -n ${CTEASFOAM_BUILD_PROXY:-} ]]; then
  # Buildx obtains registry tokens from its client process as well as dockerd.
  export HTTP_PROXY="$CTEASFOAM_BUILD_PROXY" HTTPS_PROXY="$CTEASFOAM_BUILD_PROXY"
  export http_proxy="$CTEASFOAM_BUILD_PROXY" https_proxy="$CTEASFOAM_BUILD_PROXY"
  build_args+=(--network host --build-arg "HTTP_PROXY=$CTEASFOAM_BUILD_PROXY"
    --build-arg "HTTPS_PROXY=$CTEASFOAM_BUILD_PROXY")
fi
# Explicit pull also makes base-image metadata available on fresh BuildKit
# hosts, where the build cache alone need not create a local image tag.
base_image=$(awk '$1 == "FROM" { print $2; exit }' docker/Dockerfile.official-jammy)
docker pull --platform linux/amd64 "$base_image"
docker build "${build_args[@]}" --pull --platform linux/amd64 --progress=plain \
  -f docker/Dockerfile.official-jammy -t "$image" .

echo "IMAGE METADATA"
docker image inspect "$image" --format 'ID={{.Id}} OS={{.Os}} ARCH={{.Architecture}}'
docker image inspect "$base_image" --format 'BASE_ID={{.Id}} RepoDigests={{json .RepoDigests}}'
echo "COMPILER/MPI/HEADERS/ACTUAL WMAKE COMPILE"
docker run --rm --platform linux/amd64 \
  --mount "type=bind,src=$PWD,dst=/workspace,readonly" \
  --entrypoint /bin/bash "$image" \
  -c 'bash /workspace/tests/smoke/verify_image.sh'
echo "EXPECTED FAILURE REGRESSIONS (injected errors below must be rejected)"
docker run --rm --platform linux/amd64 \
  --mount "type=bind,src=$PWD,dst=/workspace,readonly" \
  --entrypoint /bin/bash "$image" \
  -c 'bash /workspace/tests/smoke/verify_failure_paths.sh'
echo "SUCCESS: two M0 stub builds, MPI toolchain and failure regressions passed; physical solvers remain NOT IMPLEMENTED"
