# M1C-2 (2026-10-24): corrected-Laplacian MPI decomposition diagnostic

Status: **FIX CANDIDATE / NOT INDEPENDENTLY SCIENTIFICALLY ACCEPTED**. The frozen x-axis two-rank case failed before an explicit species coupled-boundary refresh and passed the targeted official-v2606 developer run after that minimal fix. Exact-commit regression and independent scientific review remain separate requirements. This report documents reproducible developer runs, not an external independent review. Parent: Draft PR #29, commit b25cdc709fa709b104ff3cb32fbb08f3f802768c; M1C-2 validation was predeclared in 392197e3a881dcce9027a362f9c08bcf78d4ca25.

Official OpenCFD v2606 Jammy container: `cteasfoam/official-v2606-dev:local-python3` (image ID prefix `0b40c3d9d1fa`). Actual 960 and 3240 cell warped 3-D meshes, U=(0.25,0,0) m/s, 40 dt=0.01 s steps, HT/HTO diffusion, oxidation and decay. Decompositions use real `decomposePar`, `mpirun -np 2`, `reconstructPar`; all compared to their original serial mesh. Solver residual was tightened to 1e-15 *in both serial and MPI* solely to exclude linear-solver stopping as the explanation; frozen physical field thresholds were not changed.

## Reproduced native field comparisons (before coupled-boundary refresh)

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

For the *original failing x* case (before coupled-patch synchronization), independent species budgets and native M1C-1 analytic checks passed, yet the **cell-wise** native concentration comparison failed and signed cyclic-face flux discrepancies reached 0.00271547 Bq/s, exceeding the 0.001 Bq/s gate. This rules out claiming scientific acceptance merely from conservation.

## Mechanism and boundaries

The tested evidence localizes the defect to a decomposition-sensitive term of the **explicit non-orthogonal diffusion correction** in the corrected Laplacian. X splitting also creates `processorCyclic` patches where the original periodic xMin/xMax cross ranks. The y and coarse-z controls suggest coupled-patch handling, but the fine-z failure means the original defect was **not** exclusive to `processorCyclic`. Directly refreshing both species' boundary values before equation assembly restores x-split parity and identifies a missing coupled-field synchronization step in the local solver; no upstream OpenFOAM fault has been established. An additional faceCorrected alternative aborted on an OpenFOAM oriented/unoriented operator mismatch; that variant is also **FAIL** and is not a validated replacement.

Do **not** silently change the frozen default `Gauss linear corrected` to `uncorrected`, alter thresholds, or substitute y-split for the required x-split. `CTEAS_M1C_DECOMPOSITION=y` in the targeted runner is strictly a diagnostic control, not M1C-2 scientific admission. Default x remains the *unchanged original frozen scientific gate*. It now passes the corrected candidate's targeted local run, but formal acceptance additionally requires exact-SHA regression and independent review.

Evidence is preserved on the authorized Ubuntu host under `/home/lhy/dev/cteasFOAM/reports/m1c2-diagnostic-20261024/`, including original native meshes, serial/MPI fields, Docker solver logs, y full checker, z failure and source/variant runner logs. These private scratch contents are not pushed to GitHub.

## Corrected boundary-state candidate: frozen x restored

A minimal solver-only change adds `HT.correctBoundaryConditions();` and `HTO.correctBoundaryConditions();` immediately **before** HT matrix assembly. It does not change the equations, discretization, case inputs, or frozen tolerances. This synchronizes coupled patch/processorCyclic values used by the explicit corrected surface-gradient term. The independent original x=(2,1,1) runner was then re-executed with the same official v2606 image:

| Mesh cells | Corrected candidate max HT/HTO L-inf Bq/m3 | RMS Bq/m3 | Max inventory delta Bq | Max flux delta Bq/s | Max serial/MPI activity residual Bq | Gate |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 960 | 8.0326572e-09 | 2.9192801e-09 | 3.0267984e-09 | 4.0745363e-10 | 2.8530849e-09 | PASS |
| 3240 | 7.1013346e-09 | 1.4980089e-09 | 4.8894435e-09 | 2.9103830e-10 | 6.0863385e-09 | PASS |

Both M1C-1 analytic checks, all 40 time-step records for two species, and the +10 Bq/m3 reconstructed HT tamper rejection passed. The original broken test (before the fix) remains in the archived local negative evidence. This is evidence of **a solver boundary-state consistency omission**, not evidence of a defect in upstream OpenFOAM. The earlier uncorrected-diffusion diagnostic did not become the shipped solver.

## Remaining acceptance

- **Candidate fixed / targeted PASS:** explicit HT/HTO coupled-patch refresh restores frozen x-processor parity on both meshes, keeps Gauss linear corrected, and passes the 40-step budgets/flux and negative control.
- **Still required:** commit this code; run official v2606 exact-commit Docker targeted tests, smoke/full and all applicable broader regressions; obtain genuinely independent scientific review tied to that SHA. No complete science admission has been granted.
- Retain Draft PR until all GitHub required checks, scientific evidence, and independent review truly pass; never apply science-validated/merge-verified from a targeted developer PASS alone.
