# cteasFOAM OpenFOAM v2606 development image

## Selected route: Ubuntu 22.04 + official OpenCFD v2606 APT development packages

```bash
bash scripts/probe_official_jammy.sh
```

The script builds `docker/Dockerfile.official-jammy` using:
- `ubuntu:22.04`;
- the **OpenCFD / Keysight** signed APT repository `https://dl.openfoam.com/repos/deb` (distribution `jammy`);
- precompiled `openfoam2606-default` and `openfoam2606-dev` packages;
- C++ compiler, OpenMPI, `wmake`, OpenFOAM FVM headers and libraries.

It then compiles and runs the M0 atmospheric and soil **placeholder** applications. Compilation **does not** validate any PDE or coupling physics.

Ubuntu 22.04 host and container are the selected baseline. GitHub-hosted CI may run its Docker build on a newer Ubuntu host while still testing the **22.04 container**. The target Ubuntu 22.04 server requires separate final MPI, NUMA and I/O acceptance tests.

## Image acceptance

Only this official-package build path is maintained. We do not compile the whole OpenFOAM source tree or fork upstream during M0.

The Dockerfile pins the Ubuntu base manifest and the official OpenFOAM package
version. The remaining Ubuntu dependency versions are recorded by the probe;
they are not a complete distribution snapshot. Bind-mount code and results
outside the container writable layer. Compilation establishes no numerical
correctness.

## Server execution and evidence

Use a native x86-64 Linux server with Docker Engine and its Buildx plugin.
On Ubuntu 22.04 these are supplied by `docker.io` and `docker-buildx`.
The script refuses Mac/ARM hosts and non-x86-64 Docker daemons. It does not
start or install Docker. Run from a normal repository checkout with Docker
access, using a narrowly scoped sudo wrapper if group access is unavailable.

Each attempt creates `reports/official-jammy-<UTC>-<unique>/`, containing the
complete build/test log, exit status, source commit, working-tree patch and
source SHA-256 manifest. Keep these ignored local reports; CI uploads them as
artifacts. The image build context excludes project files; sources are mounted
read-only and compilation happens in the disposable container's `/tmp`.

Where package downloads require a host-local HTTP proxy, explicitly set
`CTEASFOAM_BUILD_PROXY` when invoking the script. This enables host networking
for the build, sets the Buildx client HTTP(S) proxy environment, and supplies
Docker's predefined HTTP(S) proxy build arguments. A sudo wrapper must explicitly
preserve those proxy variables for the Docker client.
Configure the Docker daemon's pull proxy separately on the server. Do not put
proxy credentials in Git or in captured command output. No proxy is required
or configured by default.

The MPI smoke compiles a small independent MPI program and checks a two-rank
`MPI_Allreduce` result. This verifies the installed MPI compiler/runtime only;
parallel OpenFOAM solver numerics and all physical models remain NOT_RUN.

The acceptance probe also runs `tests/smoke/verify_failure_paths.sh` in a separate
container: an MPI version command returning 42 and a fake wmake command returning
zero while emitting an error must both be rejected. These are expected failures
in the log, explicitly marked as injected regression cases.
