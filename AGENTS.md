# Agent development instructions

Project: RRain001/cteasFOAM. All code is public. Never upload proprietary experimental data, restricted PDFs, personal credentials, or copyrighted upstream OpenFOAM source snapshots without rights review.

## Fixed baseline and boundaries
- OpenFOAM **OpenCFD/Keysight v2606**, Ubuntu x86-64 Linux server, Docker-preferred.
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

## M0
Confirm exact Docker tag pull, OpenFOAM version, `wmake`, `g++`, `mpirun`, headers, and compilation of two clearly labeled placeholder applications. Prefer official `opencfd/openfoam-dev:2606`; if unavailable/incomplete, explicitly record this and test third-party `microfluidica/openfoam:2606`. Avoid a custom full source build until necessary.
