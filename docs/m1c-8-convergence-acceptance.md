# M1C-8: one-factor-at-a-time spatial and temporal convergence audit

**This contract was frozen before writing M1C-8 case generation/analysis and running native cases.** Parent Draft PR #35 exact HEAD `55c30aa2abec7f6062d6df304ea5a42553433649`. Scope: realistic *native* OpenCFD v2606 laminar PISO numerical sensitivity of the existing idealized Gaussian activity source; does **not** assert atmospheric CFD/grid independence or scientific field validation.

## Fixed physical case

Warped unit 3-D pressure-open duct; smooth interior warp is x += 0.04 sin(pi*x) sin(pi*y) sin(pi*z) at interior mesh vertices; all outer physical faces unchanged; no-slip lateral walls; pressure xMin = +0.1 m²/s², xMax = 0; viscosity = 0.05 m²/s. Initial velocity Ux = 0.05 sin(pi*y) sin(pi*z) m/s, Uy=Uz=0, initialized from actual native mesh centers; always zero HT/HTO initial and zero external HT/HTO inlet at xMin AND xMax. kOx=.4 s⁻¹, *artificial numerical* lambdaT=.2 s⁻¹, DHT=.01 and DHTO=.02 m²/s. Gaussian source centered (0.38,0.62,0.55)m, sigma=.13m, Q_HT=60000 Bq/s, Q_HTO=12000 Bq/s, **exact time window [0.075,0.225)s**, 9000/1800 Bq in all runs. Total clock [0,0.4]s, PISO 3 pressure correctors and 2 nonorthogonal correctors. Solver, source, FVM schemes, wall/open BC, geometry function all unchanged. Same official image ID `sha256:0b40c3d9d1fac4471730a6081b85d43c8b5974e1efaa953c60ba6ea82e7b2d2a`.

## Fixed 5-run experiment matrix (independent axes)

| name | native mesh NX x NY x NZ | cells | dt s | steps | purpose |
| :-- | :-- | --: | --: | --: | :-- |
| space_lo | 6x16x6 | 576 | .02 | 20 | spatial coarse |
| baseline | 8x20x6 | 960 | .02 | 20 | shared reference |
| space_hi | 12x30x8 | 2880 | .02 | 20 | spatial fine |
| time_lo | 8x20x6 | 960 | .04 | 10 | temporal coarse |
| time_hi | 8x20x6 | 960 | .01 | 40 | temporal fine |

All 5 independent case directories must have real `blockMesh`, physical interior mesh warp and `checkMesh Mesh OK.`; `writeCellCentres` and `writeCellVolumes` from OpenFOAM. Do not reuse wrong-resolution `0/C` or `0/V`. Native HT, HTO and U final fields exist and each native physical V and mesh cell number match. Real `wmake` must compile official C++ solver; no Python-emulated CFD.

## Independent numerical measurements and gates

- Independent Python checker parses actual 3-D polyMesh vertices/cell centers/volumes, all solved time records, emitter logs, native FVM boundary HT/HTO fluxes, true native final cell fields, pressure continuity and processor-equation counts. Cross-check finite-window **per-step** Bq against analytic overlap ≤1e-9 Bq and total integrated Bq ≤1e-9 Bq, all 5 cases; native final inventory vs volume sum ≤1e-4 Bq. Each signed per-species/combined Bq residual from native inventory, release, rates and flux and native-reported budget ≤1e-4 Bq on every step; pressure continuity max ≤1e-7 m³/s. Native fields finite and not significantly negative.
- Receptors are defined **independently of mesh** by three fixed continuous physical positions `R1=(.60,.62,.55)`, `R2=(.75,.50,.50)`, `R3=(.88,.62,.55)` m and the *fixed physical Gaussian receptor support* s=.10 m. For each species and receptor, compute 0.4s concentration as `sum_i V_i exp(-|C_i-R|²/(2s²)) field_i / sum_i V_i exp(-...)`; do NOT select different arbitrary nearest cells across grids. Also capture per-cell nearest center distances for interpretation. The receptor-weighted result is **not** a validated 1.5m receptor.
- Independently compute cumulative signed **net boundary export** in Bq from native per-time flux records using actual dt, total HT/HTO final inventory Bq via native field-volume integrals and analytic source counts. Fixed metrics (3 receptors × 2 species, inventories 2 species, cumulative net boundary export 2 species) for spatial and temporal sequences. Record wall-clock per case from timed official native solver; timings are informational and not benchmark claims under shared load.
- For each sequence (low/base/high) and each recorded scalar observable, report low→base relative change `|base-low| / max(|high|,1e-9)`, base→high relative change `|high-base| / max(|high|,1e-9)`, and ratio of absolute deltas `|high-base|/max(|base-low|,1e-9)`. Report raw values and **per-observable** outcome. **Strict evidence criterion:** in both axes and for BOTH species at primary receptor R1, final concentration changes must decrease strictly: `|C_high-C_base| < |C_base-C_low|`, with low-base difference >1e-9 Bq/m³; otherwise explicitly label `CONVERGENCE_NOT_DEMONSTRATED` (not PASS), do NOT silently adjust source, mesh, thresholds, or claim overall mesh independence. Even if this passes, three-level trend is only limited *observed asymptotic improvement* and NOT broad field verification.
- Fail-closed negative controls: (a) tampering one actual final HT cell by +1000 Bq/m³ must cause independent native field-volume audit rejection, and (b) tampering one native logged emitted_HT_Bq value by +1000 Bq must fail independent source-integral check. Keep prior M1C-0..7 regression thresholds untouched. Full exact-commit official Docker 9 targeted scripts and managed Ubuntu 7-step software `full` are separate software gates; external science and independent human review remain NOT_RUN, all related PRs Draft.

**Explicit limitations:** 0.4 s 1 m synthetic duct, artificial decay, prescribed Gaussian volume kernel, laminar PISO and fixed physical receptor averaging. Even perfect convergence here would not establish actual chimney, URANS, field-scale transport or regulated dose predictions.
