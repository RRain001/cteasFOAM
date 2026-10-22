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

Acceptance evidence is an actual successful server probe or workflow log, not Dockerfile existence. Each probe preserves its source SHA-256 manifest, working-tree patch, complete build/compiler/MPI output and exit status under a unique `reports/official-jammy-*` directory. The Dockerfile pins the Ubuntu manifest digest and OpenFOAM package version; the remaining Ubuntu dependencies are not a full distribution snapshot. An independent two-rank MPI compiler/runtime probe does not establish OpenFOAM solver parallel correctness.

## Operating system

- Production host: **Ubuntu 22.04 LTS**.
- Docker container: **Ubuntu 22.04 LTS**, using official `jammy` packages.
- Docker containers share their host Linux kernel; container and host Ubuntu releases are not inherently required to match.

## Validation boundary

Neither a successful image build nor successful `wmake` compilation means the 3-D air PDE, 3-D soil PDE, or conservative transfer interface exists or has passed scientific verification.

## M0 follow-up checks

The official public key is ASCII-armored despite its `.gpg` URL, so the Dockerfile dearmors it before installing the binary `.gpg` keyring and retains `signed-by` scoping. Both solver `Make/options` include and link `meshTools`, required by the cyclic AMI headers transitively included by `fvCFD.H`. Tool version commands fail the probe directly; a failing `mpirun --version` must not be hidden inside a successful `echo`.

## Locked baseline

- Ubuntu 22.04 manifest: `sha256:0a9887931b0e38af2e3487ed39e29436b0535c4f3f52e57745baaf79d1898956`.
- Official jammy OpenFOAM package version: `2606.0~rc2-1`, preserved literally from the signed repository metadata. The name `v2606` must not conceal this package version string.
- `openfoam2606-default` and `openfoam2606-dev` are installed with exact version constraints. Their OpenFOAM dependencies use matching exact-version constraints in the official package metadata.

The native build also exposed invalid `#` comments in both `Make/options` files.
Those files are preprocessed before Make, so they now use C-style comments.
Because wmake returned zero despite preprocessing errors, the smoke probe also
rejects compiler error and missing-dependency diagnostics. Failure injection
checks cover this behavior and the previously hidden MPI-version failure.
