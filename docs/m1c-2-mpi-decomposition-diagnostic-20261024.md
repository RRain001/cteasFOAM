# M1C-2 (2026-10-24): corrected-Laplacian MPI decomposition diagnostic

Status: **FAIL / NOT SCIENTIFICALLY ACCEPTED** for the frozen x-axis two-rank requirement. This report documents experimental, reproducible developer runs, not an external independent review. Parent: Draft PR #29, commit b25cdc709fa709b104ff3cb32fbb08f3f802768c; M1C-2 validation was predeclared in 392197e3a881dcce9027a362f9c08bcf78d4ca25.

Official OpenCFD v2606 Jammy container: `cteasfoam/official-v2606-dev:local-python3` (image ID prefix `0b40c3d9d1fa`). Actual 960 and 3240 cell warped 3-D meshes, U=(0.25,0,0) m/s, 40 dt=0.01 s steps, HT/HTO diffusion, oxidation and decay. Decompositions use real `decomposePar`, `mpirun -np 2`, `reconstructPar`; all compared to their original serial mesh. Solver residual was tightened to 1e-15 *in both serial and MPI* solely to exclude linear-solver stopping as the explanation; frozen physical field thresholds were not changed.

## Reproduced native field comparisons

| Decomposition / scheme | Cells | HT max absolute (Bq/m3) | HT RMS (Bq/m3) | Interpretation |
| --- | ---: | ---: | ---: | --- |
| x=(2,1,1), Gauss linear corrected | 960 | 0.151387732942 | 0.029491326311 | **FAIL** frozen L-inf 0.001/RMS 0.0005 |
| x=(2,1,1), Gauss linear corrected | 3240 | 0.061155370087 | 0.010373603773 | **FAIL** frozen gate |
| y=(1,2,1), same corrected scheme | 960 | 3.0035153e-08 | 6.153086e-09 | PASS of full diagnostic checker |
| y=(1,2,1), same corrected scheme | 3240 | 1.8044375e-08 | 2.780912e-09 | PASS of full diagnostic checker |
| z=(1,1,2), same corrected scheme | 960 | 8.9639798e-09 | 3.624753e-09 | PASS only for coarse field |
| z=(1,1,2), same corrected scheme | 3240 | 0.025818246999 | 0.002003709008 | **FAIL** full diagnostic checker |
| x=(2,1,1), Gauss linear uncorrected | 960 | 8.0326572e-09 | 2.957146e-09 | Ablation ONLY, changed discretization; NOT frozen gate |
| x=(2,1,1), corrected + Gauss linear gradient (instead of leastSquares) | 960 | 0.150616778061 | 0.029468725362 | FAIL, gradient switch alone insufficient |

For y decomposition, the **actual full targeted runner** additionally verified both native analytic M1C-1 baselines, HT/HTO full 40-step inventories, signed xMin/xMax/net flux records, residuals, two processor meshes, reconstructed fields, and the deliberately corrupted HT field rejection. Full checker values: 960 max inventory delta 3.027e-9 Bq, max flux delta 7.858e-10 Bq/s, max MPI activity residual 2.821e-9 Bq; 3240 max inventory delta 4.424e-9 Bq, max flux delta 5.530e-10 Bq/s, max MPI residual 5.159e-9 Bq.

For the *x* case, independent species budgets and native M1C-1 analytic checks pass, yet the **cell-wise** native concentration comparison fails and signed cyclic-face flux discrepancies reach 0.00271547 Bq/s, exceeding the 0.001 Bq/s gate. This rules out claiming scientific acceptance merely from conservation.

## Mechanism and boundaries

The tested evidence localizes the defect to a decomposition-sensitive term of the **explicit non-orthogonal diffusion correction** in the corrected Laplacian. X splitting also creates `processorCyclic` patches where the original periodic xMin/xMax cross ranks. The y and coarse-z controls suggest coupled-patch handling, but the fine-z failure means the exact implementation fault has **not** been proven to be exclusive to `processorCyclic`. Neither upstream OpenFOAM nor the local solver has yet been proved defective. An additional faceCorrected alternative aborted on an OpenFOAM oriented/unoriented operator mismatch; that variant is also **FAIL** and is not a validated replacement.

Do **not** silently change the frozen default `Gauss linear corrected` to `uncorrected`, alter thresholds, or substitute y-split for the required x-split. `CTEAS_M1C_DECOMPOSITION=y` in the targeted runner is strictly a diagnostic control, not M1C-2 scientific admission. Default x retains its original expected-to-fail numerical challenge and must become PASS before M1C-2 accepts.

Evidence is preserved on the authorized Ubuntu host under `/home/lhy/dev/cteasFOAM/reports/m1c2-diagnostic-20261024/`, including original native meshes, serial/MPI fields, Docker solver logs, y full checker, z failure and source/variant runner logs. These private scratch contents are not pushed to GitHub.

## Remaining acceptance

- Reconcile the physical `corrected` non-orthogonal surface contribution at decomposed cyclic and processor faces without changing the continuum equation or frozen input.
- Rerun default **x** for 960 and 3240 cells; both fields, all 40 budgets/flux records and the corruption negative control must pass.
- Re-run official v2606 exact-commit Docker smoke/full and broader regressions; obtain a genuinely independent scientific review tied to that SHA. These are **NOT_RUN** for the next candidate commit.
- Draft PR only; never apply science-validated/merge-verified or auto-merge from the diagnostic y PASS.
