# OpenFOAM v2606 Docker image audit

## Documentation-level findings (2026-10-22)

- `opencfd/openfoam-run:2606`: **runtime only**. Official Docker Hub documentation explicitly says **NO source code or build tools**, and no tutorials. It should **not** be selected for building custom C++ solvers.
  - https://hub.docker.com/r/opencfd/openfoam-run
- `opencfd/openfoam-dev`: official development **image family** has source code and build tools, but no tutorials. **CI run 39806054597 (2026-10-22) verified that `docker pull opencfd/openfoam-dev:2606` fails: `manifest unknown`.** Only the image family's general development-tool support is established; do not use a nonexistent tag. The CI probes the community image as alternative. Actual build tests must pass before M0 acceptance.
  - https://hub.docker.com/r/opencfd/openfoam-dev
- `microfluidica/openfoam:2606`: community-provided non-slim image explicitly supports extension builds; alternative if official `-dev:2606` tag is absent or fails smoke. Review provenance before using in production.
  - https://hub.docker.com/r/microfluidica/openfoam
  - https://github.com/gerlero/docker-openfoam

**No Docker image has been runtime-tested merely by authoring this file.** GitHub Actions `docker-dev-image-probe.yml` records runtime evidence on PR/manual runs; don't label images verified before reading CI logs.

## Acceptance gate M0

1. Docker pull succeeds for exact tag and linux/amd64.
2. OpenFOAM setup resolves v2606 and `wmake`.
3. `g++`, `mpirun`, `fvCFD.H` and FV library are available.
4. Both OpenFOAM-linked placeholder executables compile and run.
5. Record immutable image digest. Compile/run tests are not physical validation.
6. User Ubuntu server later repeats build and MPI smoke, as GitHub-hosted CI does not validate user's production CPU/NUMA/I/O.

## Initial live probe receipt (2026-10-22)

- [GitHub Actions run 39806054597](https://github.com/RRain001/cteasFOAM/actions/runs/39806054597): **FAILED**.
- `opencfd/openfoam-dev:2606`: pull reported `manifest unknown` (tag absent).
- `microfluidica/openfoam:2606`: successfully pulled as `sha256:cd572ab70e3f5af55151dfa6bbc6bfdb0f963a91dfe653716568a8537bf24df7` (manifest digest reported by Docker). The initial smoke script stopped during OpenFOAM bashrc sourcing before compiler inspection; not proof of toolchain failure.
- Follow-up M0 probe fixes the startup-script failure handling. **Do not mark developer image compiled/accepted before follow-up passes.**
