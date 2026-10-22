# cteasFOAM OpenFOAM v2606 development image

## Primary (selected): Ubuntu 22.04 + OpenCFD official v2606 Debian packages

```bash
bash scripts/probe_official_jammy.sh
```

The script builds `docker/Dockerfile.official-jammy` using:
- `ubuntu:22.04`;
- the **OpenCFD / Keysight** signed APT repository `https://dl.openfoam.com/repos/deb` (distribution `jammy`);
- precompiled `openfoam2606-default`, `openfoam2606-dev`;
- C++ compiler, OpenMPI, `wmake`, OpenFOAM FVM headers and library linkage.

It then **actually compiles and runs** the M0 atmospheric and soil **placeholders**. Successful M0 compilation is **not** an implementation/validation of the PDEs or coupling.

Ubuntu 22.04 host is directly compatible with this baseline; a host OS matching the container is convenient but not required by Docker. GitHub-hosted CI on Ubuntu 24.04 still builds and tests a **22.04 container**; final production testing on user's Ubuntu 22.04 server is a separate gate.

## Optional comparative images

```bash
bash scripts/probe_dev_images.sh
```

Official `opencfd/openfoam-dev:2606` has been tested **absent** (`manifest unknown`, Actions #39806054597). The third-party `microfluidica/openfoam:2606` image did pull successfully but is not selected as the long-lived production base. This optional comparison does not block the official APT-based route.

No OpenFOAM source recompilation and no OpenFOAM upstream fork are necessary at M0.

## After M0

Pin actual digest(s), record package versions/compiler/MPI ABI, use Docker bind mounts for source/results and avoid container writable-layer output. On the production Ubuntu server, benchmark MPI affinity, NUMA and I/O separately.
