# cteasFOAM

OpenFOAM v2606-based research software for coupled atmospheric and soil tritium transport.

**M0 infrastructure and staged M1A prescribed-flow scalar regressions have
passed their scoped gates. M1B-1f laminar pressure–velocity canonical tests
through three-dimensional warped ABC meshes have passed exact-HEAD v2606 CI.
M1B-1g distorted three-dimensional nonperiodic fixed-pressure duct is an
experimental acceptance candidate. M1C-0 adds an **experimental** HT/HTO
activity reaction/transport branch with a closed-box numerical gate. M1C-1
adds a separately gated warped-3D periodic HT/HTO advection–diffusion–reaction
benchmark and finite-volume face-flux activity budget. M1C-2/3/4 add
scope-limited MPI and signed pressure-open activity checks, including real
2-rank x-decomposition and zero-flow controls; none is an independent
scientific or environmental validation. Atmospheric URANS, complete HT/HTO
physics, soil physics and air–soil coupling remain NOT IMPLEMENTED / NOT VALIDATED.**

## Intended architecture

- `src/atmosphere/`: OpenFOAM executable with prescribed-U passive scalar,
  special Couette and experimentally pressure-corrected **laminar** momentum
  modes. Optional M1C HT/HTO advection–diffusion–reaction prototype has
  **restricted numerical** verification in closed boxes, periodic warped
  3-D ducts and open-pressure 3-D ducts, but no independent field validation.
  URANS closure, realistic HT/HTO parameters and deposition are future milestones.
  Optional `auditActivityFlux true` reports the signed HT/HTO matrix
  face flux (Bq/s) and per-species/combined discrete activity residuals;
  M1C-1 verifies a periodic domain. M1C-3 adds a **separately gated**
  nonperiodic pressure-open HT/HTO inflow/outflow and actual finite-volume
  signed Bq budget benchmark, including flow reversal and zero-flow
  negative control. M1C-2 verifies periodic HT/HTO MPI and M1C-4 additionally
  checks serial/MPI pressure-open species fields, signed fluxes and budgets
  for all three boundary modes. M1C-5 introduces optional piecewise-constant
  total HT/HTO release rates (`sourceHT_Bq_s`, `sourceHTO_Bq_s`) integrated
  over the exact overlap of each Euler time step and the dimensioned
  [`sourceStart`, `sourceEnd`) window. M1C-5 defaults to a uniform volumetric
  rate. M1C-6 optionally selects `sourceShape gaussian` with dimensioned
  `sourceCenter` and `sourceSigma`, distributing Bq/s across actual native
  cell volumes by a globally MPI-normalized 3-D Gaussian spatial kernel.
  This is **not** a resolved stack or an accident-release model; source-total
  normalization does not establish plume grid convergence. See
  `docs/m1c-5-timed-emission-acceptance.md`,
  `docs/m1c-5-mpi-emission-addendum.md` and
  `docs/m1c-6-localized-gaussian-acceptance.md`. M1C-7 independently verifies
  the same localized timed source inside a genuinely warped 3-D, pressure-open
  laminar PISO flow, including signed species flux budgets, two-rank MPI and
  no-release controls; see `docs/m1c-7-open-localized-piso-acceptance.md`.
  Source + momentum integration is numerically scoped; it does **not**
  constitute atmospheric URANS, plume grid convergence or independent field
  validation. M1C-8 adds a separate five-case 3-D grid/time refinement
  audit (576/960/2880 cells and 0.04/0.02/0.01 s steps) with fixed
  Gaussian source/receiver geometry and native signed Bq budgets. Its
  independently checked R1 concentration differences shrink on both
  refinement axes, but several R2/R3 results remain non-monotonic;
  see `docs/m1c-8-convergence-acceptance.md`. This is **not** a global
  mesh-independent, meteorological or externally validated result.
  M1C-9 independently repeats the M1C-8 five native cases and three
  directional X/Y/Z refinement cases; it audits receptor averaging widths
  0.05/0.10/0.15 m and ±0.025 m streamwise receptor shifts. The original
  four non-monotone receiver metrics are preserved against independent
  prior-SHA evidence. The R2 anomaly is observation-operator-sensitive,
  while R3 spatial non-monotonicity persists across all tested averaging
  widths; see `docs/m1c-9-receiver-diagnostics-acceptance.md`. No physical
  cause is uniquely identified, and scientific validation remains NOT_RUN.
  M1C-10 independently refines x-only cells 8→12→16→24→32 with y/z fixed,
  checks the truncated-cube Gaussian source centroid/second moments against
  analytic integrals, and retains exact native HT/HTO activity conservation.
  R3 HT/HTO show a limited decreasing-difference trend, but the final
  x refinement (24→32) still changes R3 concentrations 1.89%/1.35%; source shape and
  transport discretization remain confounded. No global convergence or
  independent physics claim; see `docs/m1c-10-x-refinement-findings.md`.
  M1C-11 separately verifies source-free fully three-dimensional periodic
  HT/HTO advective-diffusive-reactive trigonometric eigenmodes with analytic
  continuum reference on native warped 8³/12³/16³ cells, with nonzero
  Ux/Uy/Uz and zero-flow controls, Bq budgets, and negative tests; see
  `docs/m1c-11-source-free-3d-analytic-findings.md`. This isolates
  transport-discretization error from M1C-10 Gaussian source quadrature,
  without establishing site-scale plume validation or URANS/soil physics.
  M1C-12 pairs native source-free/finite-duration uniform-source three-dimensional
  analytic HT/HTO transport fields on 8³/12³/16³ grids, verifying partial-step
  emissions, exact 3840/768 Bq totals, finite-step two-species reactions,
  time-order 0.01/0.005/0.0025 s and negative controls with no C++ changes.
  It does not validate spatially localized Gaussian or real-world stack release;
  see `docs/m1c-12-pulsed-analytic-findings.md`.
  M1C-13 independently compares genuine native localized 3D Gaussian HT/HTO
  plume fields to a distinct periodic Green-function time convolution
  (wrapped Gaussian; outside-cube probability 8.57e-8), auditing
  source shape, finite release Bq, stepwise budgets and real nearest-cell
  receptor diagnostics. It does not establish whole-field or fixed-point
  mesh convergence, or real stack/weather/URANS/soil physics; see
  `docs/m1c-13-local-gaussian-green-findings.md`.
  M1C-14 compares fixed physical R1/R2/R3 points via native warped-centre
  hexahedral trilinear interpolation, 3D mesh 8³→12³→16³→24³, affine
  reconstruction, split observation error and independent Gaussian Green
  response. R2 HT still differs 33.07% at 24³, and R3 HT is nonmonotone;
  no full grid independence or realistic source/soil validation is claimed.
  See `docs/m1c-14-fixed-receptor-findings.md` and geometry erratum.
  M1C-15 optionally enables `sourceSpatialQuadrature cellTet4` for native
  real-polyhedron 4-node positive tetra Gaussian volume quadrature, with
  `writeSourceProfile true` for source diagnostics. The original cell-centre
  kernel remains the default. Six native 3-D serial and one two-rank MPI
  paired comparisons conserve prescribed source Bq, but R2 HT/HTO fixed
  point errors become slightly worse on 8³/16³/24³. The alternate method
  is experimental, NOT an accuracy improvement or a validated stack source;
  see `docs/m1c-15-gaussian-cell-quadrature-findings.md`.
  M1C-16 keeps `Gauss upwind` as default but independently compares it
  with `Gauss limitedLinear 1` for HT/HTO in six true 3D native
  Gaussian-source warped-grid cases using fixed R1/R2/R3 Green references.
  LimitedLinear substantially reduces R2 concentration error but produces
  intermediate **negative native HT/HTO activity concentrations** at
  8³/16³/24³ and is SCIENTIFICALLY REJECTED. Optional
  `auditActivityExtrema true` logs exact all-step global primitive-cell
  extrema; it never clips or corrects the fields. The diagnostic software
  test only PASSes when the physical candidate is flagged REJECTED; see
  `docs/m1c-16-convection-scheme-findings.md`. No real plume/URANS/soil
  science is accepted.
  M1C-17 adds an explicitly opt-in `enforceNonnegativeActivity true`
  fail-closed safety gate for initial and post-solve native HT/HTO primitive-cell
  extrema, with no clipping, flux adjustment or change to the default solver.
  An actual v2606 8³ warped-grid Gaussian test preserves exactly identical
  upwind HT/HTO fields and budgets with the gate on/off, but aborts genuine
  `limitedLinear 1` negative HT at step 15 and negative initial HT at time
  zero. This diagnostic rejects unsafe results; it is **not** an FCT/MULES
  positivity-preserving high-resolution method or independent science
  validation. See `docs/m1c-17-positivity-gate-findings.md`.
  M1C-18 experimentally adds optional `activityConvectionBlend globalPositive`:
  two conservative native FV Euler solutions (upwind and a named high-order
  limitedLinear candidate) combined by one all-cell nonnegativity-limited
  convex coefficient per HT/HTO species, blending their native signed face
  fluxes with the same coefficient. Only serial prescribed-U cyclic warped
  meshes are in scope; unsupported MPI/open modes reject. It is **not**
  local FCT/MULES, proven general positive transport, or external science.
  See `docs/m1c-18-positivity-flux-blend-findings.md`.
  M1C-19 experimentally adds `activityConvectionBlend faceMules`:
  native OpenFOAM MULES limiter synchronizes **per-face correction fluxes**
  against the implicit positive upwind HT and HTO fields, applies the same
  limited signed flux to activity budget accounting, and keeps physical
  noncoupled inlet/outlet corrections zero. Synthetic warped 3-D Gaussian,
  2-rank periodic MPI and open nonperiodic duct benchmarks are documented
  in `docs/m1c-19-local-fct-findings.md`. Independent environmental/URANS,
  soil, real air-soil interface and non-Bot science remain NOT_RUN.
  M1C-20 adds genuine y/z 2-rank and xy 4-rank decomposition validation,
  reversed xMax inletOutlet native Bq checks, and time-step robustness probes.
  **Important negative finding:** `faceMules` with the existing corrected
  nonorthogonal diffusion becomes low-state HT negative at dt=0.0025 s;
  the solver correctly aborts. Separate uncorrected two-point diffusion
  tests at dt=0.0025/0.01 s pass, but are NOT equivalent diffusion schemes
  or proof of corrected-scheme time robustness. The original frozen
  time-step-only gate remains PARTIAL. See
  `docs/m1c-20-numerical-robustness-findings.md`.
- `src/soil/`: independent OpenFOAM 3-D porous-media and HT/HTO solver
  (M0 placeholder only).
- `src/coupling/`: nonmatching 2-D air-soil conservative exchange
  (not implemented).
- `src/models/`, `cases/`, `tests/`, `docker/`, `docs/`: units, case
  contracts, benchmarks and provenance.


## Managed Ubuntu development and merge gates (2026-10-23)

For the shared development host, first read `/home/lhy/dev/AGENTS.md`, `/home/lhy/dev/README.md` and this repository's `AGENTS.md`. Work from the actual GitHub PR HEAD in a clean isolated branch/worktree. Ubuntu Docker runs the primary software tests and exact-SHA numerical checks; GitHub provides the mandatory lightweight `repository-admission` / `scientific-admission` PR checks and selected additional CI. This host operates **Mihomo at `127.0.0.1:7897`**, with no MacBook proxy forwarding. These `/home/lhy/dev` commands are host-specific rather than prerequisites for independent clones.

The cteasFOAM verification environment reuses the pinned local OpenFOAM v2606 image; `local-ci-build` currently supports only CTEAS-next and cufotri, not cteasFOAM.

```bash
/home/lhy/dev/bin/local-ci smoke cteasFOAM --wait
/home/lhy/dev/bin/local-ci full cteasFOAM --wait
/home/lhy/dev/bin/local-ci scientific cteasFOAM --wait
```

The `local-ci` entry points currently require the repository's fixed `/home/lhy/dev/cteasFOAM/worktrees/local-docker-dev` checkout to contain the intended clean HEAD, and image/inputs/results must refer to that SHA. `full` is only the selected **software** suite; `scientific` can legitimately report `PARTIAL`/`NOT_RUN`. A source/test PR needs a verified `ubuntu-docker-full` commit status for its live SHA; the project also retains any applicable original-model, analytical, scientific or independent-review gates.

The `main` Ruleset requires PR review through `repository-admission` and `scientific-admission`. When every applicable check passes and any required independent scientific acceptance is documented, the agent adds `merge-verified` and GitHub's event-driven native auto-merge processes the PR. Draft, failed, pending, missing-data, stale-SHA or scientifically unaccepted changes remain blocked. Heavy scientific workflows can be agent-dispatched without manual GitHub UI interaction; they are **not** automatically PASS when skipped.

## Official OpenFOAM environment

Ubuntu 22.04 + official signed OpenCFD/Keysight OpenFOAM v2606 APT packages,
no third-party OpenFOAM images and no source rebuild. Although v2606 was
officially released 2026-06-26, its official Jammy Debian pool packages are
still named `2606.0~rc2-1`; this exact tested version is intentionally pinned.

On a **native x86-64 Ubuntu Docker host**, execute M0 toolchain checks:

```bash
bash scripts/probe_official_jammy.sh
```

The M1A closed-domain analytic diffusion regression belongs to the Ubuntu Docker smoke/full verification profile; the original GitHub OpenFOAM probe is retained for agent-initiated `workflow_dispatch` reproduction rather than automatic PR builds. See `docs/m1a-passive-scalar-acceptance.md` for equation, tolerance and limitations.

## Validation

Compilation is not a physical test. Gate every increment separately with
analytic or conservation benchmarks, field checks, then MPI/scaling and
independent observations. No proprietary data in this public repository.
