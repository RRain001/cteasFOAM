#!/usr/bin/env bash
# Actual build-and-compile acceptance for Ubuntu 22.04 + official OpenCFD v2606 packages.
# Usage from repository root: bash scripts/probe_official_jammy.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

mkdir -p reports
report="reports/official-jammy-dev-probe.txt"
: > "$report"
image="cteasfoam/official-v2606-dev:jammy-m0"

echo "BASE: ubuntu:22.04 / official dl.openfoam.com v2606 packages" | tee -a "$report"
docker build --pull --progress=plain \
  -f docker/Dockerfile.official-jammy \
  -t "$image" . 2>&1 | tee -a "$report"

echo "COMPILER/MPI/HEADERS/ACTUAL WMAKE COMPILE" | tee -a "$report"
docker run --rm --platform linux/amd64 \
  --mount "type=bind,src=$PWD,dst=/workspace,readonly" \
  --entrypoint /bin/bash "$image" \
  -c "bash /workspace/tests/smoke/verify_image.sh" 2>&1 | tee -a "$report"

echo "IMAGE METADATA" | tee -a "$report"
docker image inspect "$image" --format 'ID={{.Id}} RepoDigests={{json .RepoDigests}}' | tee -a "$report"
echo "SUCCESS: official package-based dev image compiled two M0 stubs; physical solvers remain NOT IMPLEMENTED" | tee -a "$report"
