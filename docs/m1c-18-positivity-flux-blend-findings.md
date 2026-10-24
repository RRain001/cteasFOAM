# M1C-18 native global convex positivity baseline — findings

Frozen acceptance: `docs/m1c-18-positivity-flux-blend-acceptance.md` (first committed at `c14ba8d`). Parent M1C-17 Draft #45 exact SHA `13efcd2569260e3336bcbae934aa6fc548e0ab7a`. M1C-18 is a research baseline and **NOT** a local flux-corrected transport (FCT), MULES implementation or proven universally bounded advection method. Default solver remains unchanged.

## Physical algorithm and limited proof

For HT, solve two actual native Euler finite-volume systems with identical prior timestep, spatial source profile/finite source-window overlap, `kOx+lambdaT` implicit losses and diffusivity; one uses `Gauss upwind`, the other the opt-in extra scheme `div(phi,HTHigh) Gauss limitedLinear 1`. Let their internal-cell solutions and native matrix face fluxes be `H_L,F_L` and `H_H,F_H`. Determine a single global coefficient

`alpha = min(1, min_{cells with H_H < 0} H_L/(H_L-H_H))` (margin `1-1e-12` only when limited),

then form `H = (1-alpha) H_L + alpha H_H` and `F = (1-alpha) F_L + alpha F_H`. Low-order negative/nonfinite inputs or candidate nonfinite values fail closed. Every corrected primitive-cell value must be >= -1e-12 Bq/m3; the existing stricter opt-in field gate also runs. No cell value clipping/rescaling.

HTO is then solved with `kOx*H` referencing the **final blended HT**; its own high/low implicit transport and final flux use the same procedure with a *separate global alpha*. Because each candidate obeys a linear native FVM budget with identical source and reaction coefficients, and its flux is blended with the **same alpha as its species inventory**, the combined method preserves each species' actual native signed activity budget up to solver roundoff. Uniform HT/HTO decay and linear conversion keep the summed tritium balance consistent. The flux is a convex combination of the two *candidate flux solutions*, **not** a recomputed single constitutive stencil on the blended state; there is no local FCT limiting or general monotonicity theorem on arbitrary cases. This mode is deliberately serial-only, prescribed-velocity-only and cyclic-only, and aborts on unsupported MPI/open boundary/momentum combinations.

## First genuine official OpenFOAM v2606 findings (pre-candidate)

The 8³/16³ native warped periodic true-cell Gaussian experiment with 80 genuine 0.005-s steps to 0.4 s, nonzero 3-D U, HT source 3840 Bq and HTO source 768 Bq, independent fixed Green-function receiver R1/R2/R3 and signed FVM budgets has passed the targeted checker (initial development run log: `/home/lhy/dev/cteasFOAM/reports/m1c18/dev-run2.log`). These are frozen previous-test numerical tolerances, not newly relaxed values.

| 3-D cells | Mode | Worst signed HT/HTO/total residual (Bq) | HT global alpha minimum | HTO global alpha minimum | R2 HT error | R2 HTO error |
|---|---|---:|---:|---:|---:|---:|
| 8³ | Baseline upwind | 4.9191726248e-9 | n/a | n/a | 62.2353% | 54.3747% |
| 8³ | Global convex blend | 8.91384388524e-9 | 0.1350374231 (7/80 constrained steps) | 0.5548088171 (4/80) | 51.4038% | 41.8192% |
| 16³ | Baseline upwind | 4.65767957536e-9 | n/a | n/a | 41.4269% | 35.0288% |
| 16³ | Global convex blend | 1.16131397832e-8 | 1 (0/80 constrained) | 1 (0/80 constrained) | 19.6055% | 13.9968% |

A later larger diagnostic (`dev-run3.log`) additionally verified 24³: upwind worst signed budget `1.3573128399e-8 Bq`, blend worst signed budget `9.86499770761e-9 Bq`; alpha=1 for both species for all 80 steps. R2 HT error decreased **33.07094% → 11.04437%** and R2 HTO **27.50314% → 7.70968%**. This is comparable but NOT identical to M1C-16's unsafe direct high-order experiment; here the limiter is assembled from a provisional upwind state. All 24³ true-cell extrema and source/Bq budgets met the unchanged thresholds, and genuine 2-rank decomposed MPI was explicitly rejected by the native serial-only guard (not claimed as MPI parity).

The 8³ global limiter *actually activates* on the native solver's high-order candidate; on 16³ it does not activate because the high-candidate assembled using the low-order field is nonnegative on that test. Intermediate extrema logs were compared against physical final fields and all 80 steps. R1/R2 relative errors improved in these restricted tests; R3 16³ HT worsened from 2.13% to 2.98%, so no broad accuracy superiority is claimed. Unlike M1C-16's direct limitedLinear solution, the high-candidate here uses the low provisional field for the non-linear limiter coefficient: their numeric solutions are not identical.

## Required remaining gates

The exact-SHA full validation, complete M1C-0..17 regression suite, open-boundary rejection evidence, and external model science/independent non-Bot review must be reported separately from the development result. The 24³ scoped and genuine two-rank rejection development evidence was collected; that evidence does not establish MPI coupled flux parity. No general-site plume, actual URANS, deposition, soil PDE or nonmatching conservative air–soil exchange is implemented or accepted. Keep Draft, no `science-validated`, no `merge-verified`, no merge.

Next numerical stage, if this baseline is accepted for comparison, is truly facewise *local*, conservative FCT/MULES with explicit coupled-boundary flux proof and MPI parity rather than a global alpha.

## Final focused development matrix before candidate freeze

`/home/lhy/dev/cteasFOAM/reports/m1c18/dev-run4.log` is a real official-v2606 container invocation, exit **0**, after genuine 8³, 16³, 24³ baseline and global-blend solves (six independent meshed cases, 80 scalar steps each). Final state and all 80 per-species Bq budget records were checked against real native `0/C`, `0/V`, `0/HT`, `0/HTO`, `0/U`, `0.4/HT`, `0.4/HTO`, native sources, physical Green-function R1/R2/R3 and the frozen 1e-5 Bq tolerance. All six cases passed. The actual limit activated 7/80 (HT) and 4/80 (HTO) steps on 8³. None on 16³/24³.

Five fail-closed controls also passed: invalid additional native FV high-order discretization; invalid `activityConvectionBlend` mode; **genuinely decomposed two-rank MPI** rejected for serial-only mode by explicit model guard; modified genuine native log with a negative alpha rejected by independent parser; **genuinely nonperiodic wall/patch mesh** rejected by native cyclic-only guard. A rejected unsupported experiment is not evidence of MPI/open-boundary solver capability.

The 24³ result improves R2 with global blending, but the scheme requires two transport solves per isotope on each timestep (four solves instead of two). It uses a **single global alpha** across all cells, and this synthetic result is not an efficiency/physics validation. It is only a stepping stone to a future conservative *facewise* FCT/MULES algorithm. Old defaults and the entire preceding scientific-rejection status of unguarded `Gauss limitedLinear 1` are unchanged.
