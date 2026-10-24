# M1C-20: frozen face-MULES robustness acceptance

Parent M1C-19 Draft PR #47 exact SHA `7ef9eb4ba814a38d24f6edfc335e1c9bb2727066`. Freeze this acceptance before modifying code or numerical validation. M1C-19's actual synthetic field and budget tests remain unchanged and required, not a stand-in for external atmospheric science.

## Required new native cases

- On genuine 8x8x8 warped, cyclic Gaussian HT/HTO initial/source data, repeat native 80-step reference at `dt=.005` and verify full 3-D serial vs reconstructed two-rank MPI for three **independently generated** decompositions: `(2,1,1)`, `(1,2,1)`, `(1,1,2)`; four ranks `(2,2,1)` if supported. Reuse the same exact solver binary, source and custom mesh. Check all 80 consecutive time records, independent signed budget/HT/HTO total checks; final reconstructed HT/HTO primitive fields, all time-step inventory compared against serial with frozen tolerances `3e-4 Bq/m3` and `1e-5 Bq`, respectively. No changing a failing tolerance.
- True 3-D warped 8³ and (if possible) 16³ periodic cases with `dt=.0025` (160 steps to 0.4 s) and `dt=.01` (40 steps to 0.4 s), change only native controlDict deltaT; confirm same total HT+HTO emissions (HT 3840 Bq, HTO 768 Bq) with actual finite-window source overlap, real field extrema nonnegative, actual signed conservation per step within `1e-5 Bq` and final R1/R2/R3 independent analytic Green checks at unchanged tolerances. Document accuracy trend: compare time refinement at identical final t, do not claim order of accuracy from two points.
- Open boundary: native 960-cell inletOutlet wall/patch geometry must retain **nonzero** signed boundary transport and per species Bq budgets, detect any boundary flux recomputation regression. Add a **sign-reversed wind** open case (`U_x=-0.2 m/s`) with physically consistent inlet/outlet BCs and source-free initial species. Positive/negative time-dependent inflows require explicit user BC; do not silently impose a wrong sign.
- Audit facewise limiter with native 80/160/40 per-step logs, ranges `[0,1]`, independent tests that tampered negative lambda or missing species/step fails, and ensure actual localized high-vs-low field differences are nonzero. No clipping or activity renormalization.
- M1C-0..19 (20 target suites) still PASS on **same HEAD**, plus new M1C-20 native targeted suite, Ubuntu `local-ci full cteasFOAM --wait` 7/7 and `local-ci-report` exact SHA. Formally scientific remains PARTIAL/NOT_RUN until independent site / field / coupled soil and URANS validation. Keep Draft unmerged and `scientific-admission` blocking.

## Limits

Variable-dt tests use *constant dt per run*, not an adaptive time-step controller. MPI scaling/performance, multi-mesh 3-D wind, heterogeneous terrain and true reactive/dry/wet soil interface remain out of scope. Any unsupported case is NOT_RUN, never PASS. This stage may add numerical software tests but not claim validated atmospheric/soil physics.
