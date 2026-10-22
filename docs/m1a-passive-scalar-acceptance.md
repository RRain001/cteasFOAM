# M1A — prescribed-flow passive scalar: scoped implementation and gates

This is the **first physical finite-volume component**, not an atmospheric
URANS model, and **not** an HT/HTO radiological transport model.

## Implemented contract
- A single OpenFOAM `volScalarField c` (dimensionless) on a three-dimensional
  `fvMesh`, with a prescribed, incompressible velocity field `U` in m/s.
- Time-dependent FVM: `ddt(c) + div(phi,c) - laplacian(D,c) = 0`;
  `phi = flux(U)`, diffusivity `D` in m²/s.
- `D >= 0`, explicit mandatory field/dictionary reads and field output.
- Integrated `c dV` log at start and after each time step. This integral is
  **only invariant for zero net boundary flux with divergence-free velocity**.

## Initial numerical benchmark (CI gate)
Case: `cases/m1a/diffusion1d`, 40 x 1 x 1 FV cells embedded in OpenFOAM,
no-flow, no-flux domain, `D=0.05 m²/s`, `dt=0.01 s`, `t_end=2 s`.
Initial condition `c(x,0)=1+0.2 cos(pi x)`.
Analytic solution: `c(x,t)=1+0.2 exp(-D pi² t) cos(pi x)`.
Independent Python checker requires max error <= 1.5e-3, closed-domain mean
drift <= 2e-8, and nontrivial amplitude decay. This verifies transient
**diffusion only**; the advective term is present but not physically validated.
All compilation and numerical PASS results require actual CI logs.

Run: `bash tests/m1a/run_diffusion_regression.sh` **inside the official
v2606 Docker image**, with the repository read-only mounted at `/workspace`.
The program and case are copied to disposable writable paths.

## Explicitly NOT_RUN / out of scope
- Advective translation and MMS convergence; nonzero inflow/outflow flux.
- 3-D obstacle wake, unsteady momentum and URANS turbulence closures.
- Two independent HT/HTO activity fields, decay, conversion, oxidation,
  wet/dry deposition, HTO boundary exchange or any radiological sources.
- Atmospheric solver MPI, decomposed mesh / restart verification.
- Soil PDE and air-soil coupling, Bq bookkeeping and multirate synchronization.

No M1A acceptance without independent GitHub Actions job SUCCESS and inspection
of actual analytic-check numerical output. Do not merge a red or unverified PR.
