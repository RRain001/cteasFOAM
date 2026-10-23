# M1A-3: Actual OpenFOAM MPI PDE equivalence

This is a genuine two-rank finite-volume diffusion benchmark, not the M0 MPI communication smoke test. The same 40-cell Neumann-wall diffusion case is solved once serially and once through OpenFOAM decomposePar along x and mpirun -np 2.

The regression calls reconstructPar -latestTime and an independent checker. The checker requires:

- Maximum reconstructed vs serial cell discrepancy <= 1e-9.
- Serial vs MPI final integrated concentration difference <= 1e-9.
- Global inventory drift across all logged time steps <= 2e-8.
- Maximum reconstructed error against the analytic diffusion solution <= 1.5e-3.
- Exactly 40 finite-volume cell outputs and a complete 200-step history.

Build/decomposition/solver/reconstruction errors are displayed in CI logs. The existing M1A diffusion and nonzero periodic-advection regression gates remain mandatory.

Explicitly NOT_RUN: nonlinear URANS flow, parallel advection, heterogeneous 3-D winds, restart and larger decomposition sizes, HT/HTO activities, soil PDE, and air-soil interface conservation. Track subsequent gates in Issue #3.
