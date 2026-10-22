# OpenFOAM v2606 official Docker development-image audit

## Verified distribution sources

- [OpenFOAM v2606 release announcement](https://www.openfoam.com/news/main-news/openfoam-v2606): precompiled Ubuntu 22.04 LTS packages are available.
- [OpenCFD jammy v2606 package directory](https://dl.openfoam.com/repos/deb/dists/jammy/main/pool/2606_0/binary-amd64/): the official `openfoam2606`, `openfoam2606-tools`, `openfoam2606-dev` and `openfoam2606-default` package names are published.
- [Official APT repository bootstrap](https://dl.openfoam.com/add-debian-repo.sh) documents the OpenCFD/Keysight APT feed. The project Dockerfile imports the key with `signed-by`; it does **not** execute the remote bootstrap shell script.
- `opencfd/openfoam-run:2606` is a **runtime-only** Docker image and does not supply the full build environment.
- The exact official Docker development-image tag `opencfd/openfoam-dev:2606` was tested on [Actions run 39806054597](https://github.com/RRain001/cteasFOAM/actions/runs/39806054597) and returned `manifest unknown`. That tag is **not** the chosen build dependency.

## Selected project build

`docker/Dockerfile.official-jammy` begins with `ubuntu:22.04`, configures the **official signed OpenCFD APT repository**, and installs precompiled `openfoam2606-default` + `openfoam2606-dev` plus C++/OpenMPI build dependencies. It does **not** rebuild upstream OpenFOAM from source.

`scripts/probe_official_jammy.sh` is the **only** Docker image probe and build acceptance route. CI attempts the Docker build and verifies OpenFOAM version/environment, `fvCFD.H`, `wmake`, compiler, MPI and actual compilation/execution of both cteasFOAM placeholder executables.

Acceptance evidence is the successful workflow log, not Dockerfile existence. After CI passes, pin base-image digest and package versions and repeat essential parallel/MPI checks on the Ubuntu 22.04 production server.

## Operating system

- Production host: **Ubuntu 22.04 LTS**.
- Docker container: **Ubuntu 22.04 LTS**, using official `jammy` packages.
- Docker containers share their host Linux kernel; container and host Ubuntu releases are not inherently required to match.

## Validation boundary

Neither a successful image build nor successful `wmake` compilation means the 3-D air PDE, 3-D soil PDE, or conservative transfer interface exists or has passed scientific verification.
