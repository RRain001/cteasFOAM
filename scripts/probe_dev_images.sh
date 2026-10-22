#!/usr/bin/env bash
# Usage: bash scripts/probe_dev_images.sh
# Probe upstream official development image first, community development image second.
set -eo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
images=(opencfd/openfoam-dev:2606 microfluidica/openfoam:2606)
mkdir -p reports
report="reports/docker-image-probe.txt"
: > "$report"
selected=""
for candidate in "${images[@]}"; do
  echo "===== CANDIDATE $candidate =====" | tee -a "$report"
  if ! docker pull "$candidate" 2>&1 | tee -a "$report"; then
    echo "PULL FAILED: $candidate" | tee -a "$report"
    continue
  fi
  if docker run --rm --platform linux/amd64 \
      --mount "type=bind,src=$PWD,dst=/workspace,readonly" \
      --entrypoint /bin/bash "$candidate" \
      -c "bash /workspace/tests/smoke/verify_image.sh" 2>&1 | tee -a "$report"; then
    selected="$candidate"
    echo "SELECTED: $selected" | tee -a "$report"
    docker image inspect "$candidate" --format 'DIGEST: {{json .RepoDigests}}' | tee -a "$report"
    break
  else
    echo "TOOLCHAIN/COMPILE FAILED: $candidate" | tee -a "$report"
  fi
done
if [[ -z "$selected" ]]; then
  echo "FAIL: neither candidate provided a working v2606 C++ build environment" | tee -a "$report"
  exit 1
fi
echo "SUCCESS: development-image compilation smoke completed; no physical validation performed" | tee -a "$report"
