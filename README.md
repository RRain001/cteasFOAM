# cteasFOAM

OpenFOAM v2606-based research software for coupled atmospheric and soil tritium transport.

**Status: M0 bootstrap. No atmospheric, soil, or coupling physics has passed validation.**

## Intended architecture

- `src/atmosphere/`: independent OpenFOAM 3-D transient URANS + HT/HTO solver.
- `src/soil/`: independent OpenFOAM 3-D porous-media flow/transport and HT/HTO solver.
- `src/coupling/`: 2-D air–soil interface, conservative activity transfer (Bq), nonmatching mesh mapping and multirate scheduling.
- `src/models/`: common physics, units and contracts.
- `cases/`, `tests/`, `docker/`, `docs/`: reproducible demonstrations, verification, build and validation evidence.

Both solvers are based on **OpenFOAM v2606 FVM libraries**, but have independent meshes, state, and time steps. A soil solver is **not** a from-scratch FVM library.

## Selected Docker development baseline

**Ubuntu 22.04 base image + official OpenCFD/Keysight v2606 APT packages.** Your Ubuntu 22.04 server needs no OS upgrade. `docker/Dockerfile.official-jammy` installs prebuilt OpenFOAM development packages, C++ and MPI; it does not compile upstream OpenFOAM source.

From a Linux machine with Docker:

```bash
bash scripts/probe_official_jammy.sh
```

The probe builds the image and attempts real `wmake` compilation of both **M0 placeholders**. See CI and `docs/docker-image-audit.md` for actual acceptance status. **Only the official-package Docker build is maintained.**

## Validation principles

A compiling executable is **not** a validated solver. Future milestones require separate atmospheric, soil and interface conservation tests, followed by observational validation. Only publish datasets/software where rights permit it.
