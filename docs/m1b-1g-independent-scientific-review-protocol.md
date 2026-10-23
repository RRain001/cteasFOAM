# M1B-1g independent scientific review protocol

Status: **REVIEW REQUEST PREPARATION — NOT APPROVED / NOT VALIDATED**.

This document describes what a distinct, technically qualified reviewer
would need to verify before M1B-1g can be independently accepted.
It is **not** a review performed by the author, ChatGPT, CI, or a bot.
The GitHub `scientific-admission` control requires a current-HEAD
non-bot independent APPROVED review and a separately justified
`science-validated` determination. Do not apply that label on the
strength of this document or of software PASS alone.

## Exact source and physical scope

- Obtain the **live** PR #20 commit SHA from GitHub. Never substitute an
  old hash; compare `git rev-parse HEAD`, downloaded artifact checksums
  and manifest `source_sha` to the live PR HEAD.
- Environment: x86-64 Ubuntu, official signed OpenCFD/Keysight OpenFOAM
  **v2606** packages (`openfoam2606=2606.0~rc2-1`), not an unofficial
  distribution. Official locally pinned image ID is recorded in the PR
  evidence. Do not publish private host logs or restricted materials.
- Actual physical test is **incompressible laminar** unit-cube 3-D duct,
  nonperiodic pressure patches xMin=0.1 and xMax=0 **m²/s²** (kinematic
  pressure), four no-slip solid walls, `nu=0.05 m²/s`, D=0,
  initial U=(0,0,0), passive scalar c=1, true interior-only warp
  `x'=x+0.04 sin(pi*x)sin(pi*y)sin(pi*z)`. Its exterior walls, inlet
  and outlet are unwarped.
- Numerical method in candidate: OpenFOAM finite volume, Euler time,
  Gauss upwind convection, Gauss linear corrected Laplacian and
  snGrad, PISO 3 outer pressure correctors and 2 additional inner
  nonorthogonal correctors, actual flux-corrected `phi`.
  Source: `src/atmosphere/cteasAtmosphereFoam.C`.
- **Conspicuously excluded:** RAS/URANS, turbulence closures,
  buoyancy, dynamic meteorological inflow, building-wake validation,
  HTO/HT activity conversion, dry/wet deposition, soil PDE,
  air–soil nonmatching coupling, external observational validation.

## Evidence reproduction matrix

| Independent reviewer action | Minimum observable evidence | Interpretation |
|---|---|---|
| Inspect all C++ source, test cases and initial fields | No source-code modification between cases, correct units and physical boundary values | Verification of intended equations; not field validation |
| Execute `run_warped_open_pressure.sh` | Native `blockMesh/checkMesh/writeCellCentres`, true warped and flat 3-D cases plus Δp=0 negative control; physical signed inlet/outlet flux, 20 steps and 180 p solves | Conservation/analytic smoke |
| Execute `run_warped_open_pressure_mpi.sh` | Genuine two-rank x/z MPI and zero-drive parallel control; real reconstructions of U, absolute p, flux, scalar | Decomposition consistency |
| Execute `run_warped_open_pressure_timestep.sh` | Same physical mesh and BCs, Δt=.02/.01/.005, 20/40/80 full steps, 180/360/720 p solves; predeclared U and Q difference ratios ≤0.85 | Fixed-grid temporal self-convergence only |
| Execute `run_warped_open_pressure_spatial.sh` | Three **actual** warped 3-D grids 960/3240/7680 cells, common Δt=.01, 40 steps; native-centroid axial RMS vs independent series; adjacent error ratios <0.9, finest/coarsest ≤0.7 | Controlled spatial error trend, not formal GCI |
| Execute `run_icofoam_cross_solver.sh` | Run real unmodified official upstream `icoFoam` binary with byte-identical physical inputs, 20×9 PISO p solves; compare final U, **absolute** p and zero-drive response; injected corruption must fail checker | Upstream-implementation cross-check; not external/independent scientific approval |
| Execute `run_warped_abc_mpi.sh` | True 3-D periodic warped ABC with orthogonal MPI decompositions and existing analytic/continuity gates | Prior-model regression |
| Inspect managed Ubuntu Docker `local-ci full` manifest | `PASS_SCOPED/0` for exact source SHA, image ID, all selected steps, evidence hash | Selected software regression only |

Reviewer should run from a clean clone at **their own controlled
environment** when feasible rather than relying solely on author's
logs or CI status. Cases are standalone in the public repository;
external third-party numerical or field data are not embedded.

## Questions requiring explicit scientific judgement

1. Are governing momentum and pressure-flux equations, pressure
   boundary conditions and nonorthogonal correction placement
   consistent with the intended OpenFOAM v2606 formulation?
2. Is the Stokes square-duct analytic solution appropriate for a
   finite-length **nonlinear** Navier–Stokes computation? Where does
   it omit end effects, convection and actual pressure corrections?
3. Are each of the analytic-velocity comparisons at true native
   OpenFOAM cell centres; does comparing differently sampled RMS
   errors sufficiently establish a refinement trend, rather than
   true continuum convergence or a formal order of accuracy?
4. Why does fixed-grid temporal refinement slightly **increase**
   total analytic U error even as successive numerical differences
   decrease? Distinguish spatial and temporal error contributions.
5. Does the official `icoFoam` agreement constitute implementation
   consistency only, because both use the same OpenFOAM finite-volume
   algorithms? What genuinely *external* benchmark would be required
   for the proposed eventual atmospheric URANS application?
6. Are there unintended changes to `momentumMode prescribed` or
   `laminarShear` and earlier validation paths on current HEAD?
7. Is there adequate basis for generalizing a mild, interior-warped
   hexahedral case (max nonorthogonality about 7°) to higher
   nonorthogonal terrain/body-conforming meshes? **No such
   generalization is presently authorized.**

## Outcomes and GitHub process

- **Accept scoped numerical evidence** only after the reviewer
  reproduces/inspects all applicable tests, documents assumptions,
  unresolved limitations and cites the precise PR HEAD. Acceptance
  must explicitly restrict claims to these canonical laminar tests;
  it is not a certification of planned future physics.
- **Request changes** if the numerical physics, tolerances, budgets,
  input provenance, boundary conditions or generalization claims are
  unsupported. The author must correct and resubmit a new exact-HEAD
  verification; old review approval and statuses must not carry over.
- If no independent human/non-bot reviewer is available, record the
  blocker honestly as **NOT_REVIEWED**. Never self-label, set
  `merge-verified`, disable required checks or merge the Draft PR.
- The PR author must decide separately whether the project-specific
  `science-validated` scope is limited to canonical M1B-1g numerical
  verification or requires an external independent benchmark. That
  decision cannot be inferred from an automated passing test.

This review checklist contains no reviewer signature or endorsement.

## Reversed-gradient reproducibility addendum (M1B-1g-5)

The same exact-head reviewer execution of
`tests/m1b1/run_icofoam_cross_solver.sh` additionally requires
a negative kinematic inlet pressure of -0.1 m²/s² with zero at xMax
on the unchanged 960-cell warped duct. Confirm actual reverse-signed
inlet/outlet flow, uniform passive-scalar budget, 180 real native
pressure solves, signed Stokes reference, full-field U/absolute p
cross-solver parity, and the forward/reverse bounded low-Re residuals.
Observe that both the custom and official implementations still share
the OpenFOAM finite-volume infrastructure; this is not independent
external validation. The predeclared numerical bounds and caveats are
in `docs/m1b-1g-warped-open-pressure-acceptance.md`.
The need for an actual independent non-bot review, a separate
`science-validated` judgement, and Draft status remains unchanged.

## M1B-1g-6: nonuniform tracer and discrete mass-accounting evidence

A qualified independent reviewer should also execute
`bash tests/m1b1/run_open_pressure_tracer.sh` on the exact PR HEAD in
official v2606. Inspect the **actual** unmodified C++ scalar equation
`ddt(c)+div(phi,c)-laplacian(D,c)=0`, the optional
`auditPassiveScalarFlux` implementation, `fvSchemes`
`fluxRequired(c)`, the pressure-patch inflow switches in the scalar
`inletOutlet` conditions, and dimensional consistency of all terms.

Required reproducibility observations: three real warped 3-D
forward/reverse/zero pressure cases, zero initial c, switched c=1
*only on the active inflow*, D=0, nonuniform final c, same native mesh
and initial U, 20 steps/180 p solves each, correct signed inward and
outward volume/solute fluxes, whole-domain M(t) and discrete
`M_new-M_old+dt ΣF_c=0` to the a priori tolerance. Confirm that a
deliberate +0.01 m³ modification of the **printed** budget residual is
detected by the independent Python parser. Evaluate whether the
matrix-flux diagnostic is sufficiently independent from the solved
equation for the claims being made; consider externally computed face
fluxes as a distinct later validation. Neither this diagnostic nor
zero-control positivity is independent external field evidence.

This test uses a dimensionless passive tracer, not Bq/m³ activity
or any verified HT/HTO physics. The independent review, external
experimental validation and `science-validated` gate remain
**NOT_RUN / NOT_APPROVED**; do not merge based on a Python checker.

## M1B-1g-7: positive-diffusivity reproducibility addendum (NOT REVIEWED)

An independent reviewer must also run the new D=0 versus D=0.01
six-case warped-duct comparison on the exact PR SHA, inspect native
`cEqn.flux()` for both advective and diffusive contributions, and
independently verify the signed integrated inventory and physical
boundary conditions. The declared sensitivity thresholds and passive
flow non-feedback criteria are in the M1B-1g acceptance document.
The checker intentionally rejects a falsified D=0.01 input dictionary;
this is a test of the checker, not independent scientific validation.
No reviewer approval or external field validation is implied.
