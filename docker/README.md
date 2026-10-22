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

When the image is validated, pin its base image digest and installed OpenFOAM package versions, record compiler/MPI ABI, bind-mount code and results outside the container writable layer, and avoid assuming a successful compilation establishes numerical correctness.
