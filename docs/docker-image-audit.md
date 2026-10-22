# OpenFOAM v2606 Docker image audit

## Official evidence

- [OpenFOAM v2606 release announcement](https://www.openfoam.com/news/main-news/openfoam-v2606): **official precompiled Ubuntu 22.04 LTS packages** are available.
- [OpenCFD jammy v2606 package directory](https://dl.openfoam.com/repos/deb/dists/jammy/main/pool/2606_0/binary-amd64/): `openfoam2606`, `openfoam2606-tools`, `openfoam2606-dev`, `openfoam2606-default` are published.
- [Official signed APT repo bootstrap](https://dl.openfoam.com/add-debian-repo.sh) documents the OpenCFD/Keysight APT feed. Our Dockerfile uses a scoped `signed-by` key rather than remotely executing a bootstrap script.
- `opencfd/openfoam-run:2606` is **runtime**, not a full custom C++ development environment.
- `opencfd/openfoam-dev:2606` was **actually tested absent**: [Actions 39806054597](https://github.com/RRain001/cteasFOAM/actions/runs/39806054597), `manifest unknown`. No Docker Hub tag guessing.

## Primary chosen development route

`docker/Dockerfile.official-jammy` starts from `ubuntu:22.04`, imports the official scoped APT key, installs **prebuilt** `openfoam2606-default` and `openfoam2606-dev`, and tests with `scripts/probe_official_jammy.sh` via GitHub Actions. Does **not** compile upstream OpenFOAM from source.

M0 must verify actual APT resolution, release/version, `fvCFD.H`, `wmake`, `g++`, `mpirun`, and compiling/running both cteasFOAM placeholder executables. A successful CI is not evidence of physics correctness or the user's production-host MPI scaling. After acceptance, pin image digest and Debian package versions.

## Third-party alternative (not production primary)

`microfluidica/openfoam:2606` is based on the OpenCFD official package repository; its [public Dockerfile](https://github.com/gerlero/docker-openfoam/blob/main/Dockerfile) and [CI matrix](https://github.com/gerlero/docker-openfoam/blob/main/.github/workflows/ci.yml) provide provenance. It was downloaded by [Actions 39806054597](https://github.com/RRain001/cteasFOAM/actions/runs/39806054597) as:
`sha256:cd572ab70e3f5af55151dfa6bbc6bfdb0f963a91dfe653716568a8537bf24df7`.
Initial test script failed in shell startup before compiler verification. Community image testing remains optional and is **not** a validated production build.

## OS distinction

- **Production host**: Ubuntu 22.04 LTS (per user).
- **Container**: Ubuntu 22.04 (selected for predictable deployment and official `jammy` packages).
- Containers and hosts **do not need the same Ubuntu release**: Linux containers share the host kernel. The user does **not** need to upgrade Ubuntu to run OpenFOAM v2606.
