# M1C-17 — Opt-in fail-closed native HT/HTO positivity gate

The acceptance contract was frozen first at `8c2354f` from M1C-16 Draft PR #44 (`864e0510765239849e0974a8d95065434afdd47a`). This branch adds `enforceNonnegativeActivity true` only to the existing `htHto` prototype. The option is default-off; it checks real primitive-cell HT/HTO at initial time and after every coupled species solve. It globally reduces counts of nonfinite values and minima/maxima over MPI ranks, and terminates with `CTEAS_M1C_POSITIVITY_FAILURE` on nonfinite or below `-1e-8 Bq/m3`. The threshold is an **unsafe-event detection criterion**, not a positivity theorem. No values are clipped, no conservative flux correction is implemented, and the default `Gauss upwind` discretization is untouched.

## First actual v2606 focused native result (pre-freeze source)

The official `cteasfoam/official-v2606-dev:local-python3` image ran `bash tests/m1c/run_positivity_gate.sh` on a genuine 8³ warped cyclic OpenFOAM polyMesh. Finite-window Gaussian source: HT 60000 Bq/s, HTO 12000 Bq/s, t=[0.073,0.137] s, prescribed nonzero 3-D U, 80 steps at dt=0.005 s, native Bq balances and fixed physical Green-function receptor checks. Executed tests yielded:

- Unmodified upwind with switch off and on: both real 80-step solves completed; HT and HTO final native field maximum differences **exactly 0 Bq/m3**; signed per-species/combined Bq maximum residual **4.919172625e-9 Bq**; genuine injected source totals 3840/768 Bq; R1/R2/R3 independent reference checks PASS in both cases.
- `Gauss limitedLinear 1` with fail-closed gate: actual native `HT` concentration reached **-2.52317072302958 Bq/m3** on step **15**, phase `postSolve`. Native solver halted with nonzero exit and explicit reason, before reporting an accepted 80-step solution. This is the M1C-16 physically rejected candidate, NOT a bounded/acceptable high-order result.
- Real native time-zero HT `-1 Bq/m3`: rejected in phase `initial`, with **zero** time steps attempted.
- Tampering a genuine native failure log to remove the marker is rejected by the independent Python checker; a successful exit for either real negative case would make the shell test fail.

The direct development test returned exit 0 after explicitly detecting the two *expected physical failure exits*. It is **PASS for fail-closed guard behavior**, **not PASS for limitedLinear science**. These numbers characterize the pre-final C++ version; exact-SHA retest and managed full results must be recorded separately in PR and Ubuntu run manifests before any software status is published.

## Scientific boundary

No high-resolution positivity-preserving correction (FCT/MULES), no general flux sign proof, no positive boundary interpolation theorem, no independent outdoor plume validation, no atmospheric URANS, no soil PDE and no air–soil exchange are established. M1C-16 limitedLinear remains SCIENTIFICALLY REJECTED. The appropriate next physics step is a conservative bounded high-resolution flux correction with complete signed HT/HTO activity accounting and actual MPI/open-boundary tests, not value clamping.
