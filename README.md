# cteasFOAM

OpenFOAM v2606-based research software for coupled atmospheric and soil tritium transport.

**Status: M0 bootstrap. No atmospheric, soil, or coupling physics has passed validation.**

## Intended architecture

- `src/atmosphere/`: independent OpenFOAM 3-D transient URANS + HT/HTO solver.
- `src/soil/`: independent OpenFOAM 3-D porous-media flow/transport and HT/HTO solver.
- `src/coupling/`: 2-D air–soil interface, conservative activity transfer (Bq), nonmatching mesh mapping and multirate scheduling.
- `src/models/`: common physics, units and contracts.
- `cases/`, `tests/`, `docker/`, `docs/`: reproducible demonstrations, verification, build and validation evidence.

Both solvers will be based on **OpenFOAM v2606 libraries**, but have independent meshes, state and timesteps. A soil solver **is not** a separate from-scratch FVM library.

## Docker image selection

The official `opencfd/openfoam-run:2606` contains OpenFOAM runtime binaries/libraries **but explicitly omits source and build tools**. The official `opencfd/openfoam-dev` series is described as including sources and development tools. **Availability and contents of the exact `2606` tag require a live pull and compile test**. Do not infer that a tag exists from an image family overview.

The initial M0 PR includes a GitHub Actions probe that tests `opencfd/openfoam-dev:2606` first and independently probes `microfluidica/openfoam:2606` as a fallback. No homemade OpenFOAM source build is needed unless both development-image choices fail.

See `docs/architecture.md` and `docs/docker-image-audit.md` after M0 is merged.

## Validation principles

A compiling executable is **not** a validated solver. Future milestones require separate atmospheric, soil, and interface conservation tests, followed by observational validation. Only publish datasets/software where rights permit it.
