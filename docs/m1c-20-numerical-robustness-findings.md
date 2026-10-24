# M1C-20 — MULES multiple-MPI / time-step / reverse-open robustness investigation

Parent M1C-19 Draft PR #47 exact SHA `7ef9eb4ba814a38d24f6edfc335e1c9bb2727066`. Acceptance frozen first in `7199bc484fce0fcaaa819528b73dab8af9b25b3d`. **No transport solver C++ code or M1C-19 acceptance thresholds were changed.** This is a numerical investigation and extension of real OpenFOAM v2606 research tests, NOT atmospheric/soil science acceptance.

## 1. MPI partition directions (PASS_SCOPED)

Using actual warped periodic 8³ native Gaussian, 80 × 0.005 s, independent decompositions, real native OpenFOAM MPI solves and `reconstructPar`, HT and HTO serial-vs-MPI cell fields, all time-step species/summed Bq residuals, emitted Bq and native face limiter logs:

| Ranks | Decomposition | Max HT field difference (Bq/m³) | Max HTO field difference (Bq/m³) |
|---|---|---:|---:|
| 2 | y (1,2,1) | 5.90451963944e-6 | 6.1523239836e-7 |
| 2 | z (1,1,2) | 2.72433999271e-6 | 5.64828397387e-7 |
| 4 | xy (2,2,1) | 2.85300484393e-6 | 6.12675599143e-7 |

All passed unchanged 3e-4 Bq/m³ reconstructed field and 1e-5 Bq inventory/budget tolerances. Original 2-rank x decomposition and all 20 M1C-0..19 real native suites remain required.

## 2. Genuine corrected-diffusion time-step failure (IMPORTANT, NEGATIVE)

The frozen standard M1C-19 `faceMules` case uses OpenFOAM `laplacianSchemes { default Gauss linear corrected; }`. Running the *same native 8³ Gaussian mesh, identical sources, source-window, physics, positivity gate and coefficients*, but shrinking Δt from 0.005 to **0.0025 s**, triggers `CTEAS_M1C19_INVALID_LOW species=HT minimum=-2.75241736514867e-06 Bq/m3` near time t=0.08 s after the implicit low-order HT step. The model correctly fails closed rather than clipping. This is a counterexample to the assertion that *upwind convection plus corrected nonorthogonal implicit diffusion* is always a nonnegative low-order MULES baseline. The first two actual failing development logs are private Ubuntu `reports/m1c20/dev1.log` and `dev2.log`.

This counterexample **must not be omitted or reclassified as PASS**. M1C-20 preserves it as an expected-failure native negative test. The standard corrected diffusion scheme did NOT pass the frozen time-step-only matrix. Thus M1C-20 is **PARTIAL against the initially frozen scientific robustness acceptance** even when all its extended software tests pass.

## 3. Remedial numerical comparison, not equivalent diffusion physics (PASS_SCOPED_DIFFERENT_SCHEME)

We tested the same exact 8³ Gaussian physical parameters at two additional *constant* time steps with **`laplacianSchemes { default Gauss linear uncorrected; }` for these exploratory cases only**. This exchanges the nonorthogonal corrected Laplacian for a two-point uncorrected positive-diffusion stencil; it can sacrifice spatial accuracy on skewed grids and is NOT interchangeable with the default corrected diffusion physics. The tests separately verify all real time records, finite window emissions HT 3840 Bq and HTO 768 Bq, per-step signed Bq budgets ≤1e-5 Bq, genuine native MULES per-face lambdas in [0,1], 3D cell extrema and independent physical Green R1/R2/R3 receptor errors under the existing 1.5 cap:

| Scheme | Δt (s) | Genuine steps | R2 HT rel. error | R2 HTO rel. error |
|---|---:|---:|---:|---:|
| Gauss linear uncorrected | 0.0025 | 160 | 52.8145% | 43.3303% |
| Gauss linear uncorrected | 0.01 | 40 | 54.1404% | 44.7143% |

For comparison, standard corrected-diffusion M1C-19 at 8³/0.005 passed its original test, while the 0.0025 corrected-diffusion run fails as above. The different schemes must never be collapsed into a claimed three-point time-convergence order. A scientifically better next step is a **positivity-bounded nonorthogonal diffusion correction**, with a genuinely M-matrix low-order diffusion operator and facewise limited high-order diffusive correction, rather than relying exclusively on a convective MULES flux limiter.

## 4. Physically reversed nonperiodic inlet/outlet (PASS_SCOPED)

One actual 960-cell nonperiodic OpenFOAM duct uses prescribed U=(-0.2,0,0) m/s, real boundary xMax incoming HT=1e6 Bq/m³ and HTO=2e5 Bq/m³, xMin outgoing. Both species use real inletOutlet boundaries. At 20 time steps it verifies xMax **negative signed entering flux** (max incoming magnitude approx 358713 Bq/s), genuine HT, HTO and summed time-step signed boundary activity residuals ≤1e-5 Bq and nonnegative final real scalar cells. The previous M1C-19 xMin inlet direction remains an independent required test. This demonstrates both steady boundary-flow directions, not arbitrary time-varying meteorological wind or dynamic patch switching.

## 5. Formal gate and next direction

This only tests synthetic prescribed velocity, analytic/reaction/diffusion numerical transport, not independent field science, full URANS, soil PDE, deposition or nonmatching air–soil conservative interface. The same-head Ubuntu Docker Full and 21 historical+new targeted suites must PASS before publishing software evidence. An external non-Bot science review is still NOT_RUN; PR remains Draft, scientific-admission blocking, no merge or scientific labels. Next M1C-21 must explicitly fix the corrected nonorthogonal diffusion low-state negativity with a conservative positive low diffusion stencil + bounded face-local diffusive correction, *then* repeat the original `dt=0.0025` case without changing its mathematical high-order target.

Retain the source scheme and records; no broad mesh/time-step robustness has been claimed.
