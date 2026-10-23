# Agent development instructions

## Owner's current CI and merge decision (2026-10-23)
- Retain automatic pull_request GitHub Actions OpenFOAM v2606 checks; no owner-initiated manual click is required. Prior manual-only Draft PR #21 was closed.
- Continue local Ubuntu Docker compile/smoke/full and targeted analytic conservation checks first; each candidate PR then triggers automatic relevant GitHub Actions on its exact HEAD.
- After every applicable software, analytic, independent validation and GitHub check passes, automatically mark PR Ready and squash merge or enable native auto-merge without a separate owner request. Do not infer HT/HTO, soil physics, coupling or URANS scientific PASS from existing M0/M1 tests. Draft PR #20 remains unmerged until its actual scope is accepted.


Project: RRain001/cteasFOAM. All code is public. Never upload proprietary experimental data, restricted PDFs, personal credentials, or copyrighted upstream OpenFOAM source snapshots without rights review.

## Fixed baseline and boundaries
- OpenFOAM **OpenCFD/Keysight v2606**, Ubuntu 22.04 x86-64 Linux server, Docker-preferred.
- Two **separate executables**, both linked against OpenFOAM v2606 C++ FVM libraries:
  - 3-D atmospheric URANS + HT and HTO transport.
  - 3-D porous soil PDE: HT oxidation, water/HTO transport, gas–liquid exchange.
- Actual 2-D air–soil physical interface with nonmatching meshes, conservative HT/HTO activity exchange and multirate synchronization.
- Plant, food chain and dose are separately validated kinetic models.
- Do not rewrite OpenFOAM FVM/MPI infrastructure or fork upstream until a specific upstream limitation is demonstrated.
- Source of truth for thesis scope and timing remains private `career-plan`; do not copy private documents here.

## Verification discipline
- Mark skeletons as **NOT IMPLEMENTED**. Never imply a solver exists because `wmake` succeeds.
- Progress must be evidenced: compilation -> analytic/conservation tests -> independent field validation -> HPC scaling.
- HT and HTO have independent activity balances; net conversion preserves total tritium activity before radioactive decay and boundary fluxes.
- Positive interface sign is **air to soil**. Transfer signed integrated activity in Bq, not just concentration, on each coupling window.
- On public GitHub, do not import restricted MUST/JAEA/other nonredistributable experimental documents.

## M0: only official OpenFOAM packages
- Build `docker/Dockerfile.official-jammy` from `ubuntu:22.04` and install the precompiled OpenCFD v2606 development packages from the official signed APT repository.
- Verify `wmake`, `g++`, `mpirun`, OpenFOAM headers, and actual compilation of two clearly labeled placeholder applications.
- Do not rebuild the complete OpenFOAM source tree unless a demonstrated issue requires it.
- Do not add external OpenFOAM container distributions as an alternate deployment path without a new explicit decision.
