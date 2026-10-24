# M1C-5 addendum: uniform Bq/s source under two-rank MPI

This additional scoped acceptance is frozen **before** running the MPI release case. Same official OpenCFD v2606, same 4×3×2 actual 3-D unit cube, HT/HTO pulse, kOx, lambdaT, 20 time steps and 0.075–0.225 s window as the separately frozen M1C-5 serial acceptance. Decompose x into exactly 2 ranks (simple (2 1 1)); native decomposePar, mpirun -np 2 and reconstructPar (not a synthetic log). Both deconstructed meshes and final species fields must exist.

1. Independently run existing M1C-5 20-step source/inventory/closed-budget/initial/native-volume checker on the actual parallel reconstructed source case; all of its frozen BE and Bq tolerances apply unchanged (1e-6 Bq for species).
2. Independently compare serial and MPI final 24 native HT/HTO cell values, maximum absolute difference ≤1e-6 Bq/m³; 20 logged HT/HTO inventories differ by ≤1e-6 Bq. MPI source-overlap and emission histories independently obey the same 0.075–0.225 s integrated input.
3. A missing processor partition, nonzero MPI solver return or missing reconstructed field is FAIL. No change to source physical rate, mesh, concentration or frozen numerical tolerances to secure PASS. Keep this as scoped numerical research, never external scientific validation.
