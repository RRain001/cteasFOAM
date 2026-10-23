# M1C-2: two-rank MPI HT/HTO transport parity (predeclared)

Status: CANDIDATE / NOT SCIENTIFICALLY ACCEPTED. Stacked on Draft PR #29 exact SHA `b25cdc709fa709b104ff3cb32fbb08f3f802768c`; prior PR #28 and #20 remain untouched. Commit this gate before numerical implementation/runs; no tolerance changes to obtain PASS.

## Physical and numerical scope

Use the already predeclared M1C-1 real OpenFOAM v2606 warped periodic 3-D HT/HTO advection–diffusion–reaction case at both 960 and 3240 cells: U=(0.25,0,0) m/s, DHT=0.01 and DHTO=0.02 m²/s, kOx=0.4 s⁻¹, lambdaT=0.2 s⁻¹, Euler dt=0.01 s, 40 steps to 0.4 s. Run a serial and a genuine two-rank MPI solve with the same mesh, fields, schemes and decomposition along x. Preserve both independent species budgets including FVM boundary fluxes. Use `decomposePar`, `mpirun -np 2` and `reconstructPar`; verify processor0/1 meshes and real rank-specific field files exist. Do not replace MPI with two serial runs or synthetic results.

## Frozen acceptance

For each mesh, the original independent M1C-1 checker must pass the serial and reconstructed parallel cases (including analytical errors, cyclic flux cancellation and activity budgets). Each rank must advance 40 real steps; master log must have 40 HT/HTO solves and 40 species flux/budget records, without MPI fatal errors. Independently compare *all* reconstructed native final HT/HTO cell values to the serial values in the same cell order: max absolute difference <= 1e-3 Bq/m³, RMS difference <= 5e-4 Bq/m³. Compare all 40 per-step global HT/HTO inventories: max difference <= 1e-4 Bq; and all 40 xMin/xMax and net HT/HTO flux records: max difference <= 1e-3 Bq/s. Both serial and MPI reconstructed per-step activity residuals must remain <= 1e-5 Bq. Negative control: deliberately change one reconstructed final HT cell by +10 Bq/m³; independent comparison must fail. All comparisons require actual official native output and physical Bq/m³ field dimensions.

Run exact-SHA Docker targeted test, smoke/full, M1C-0 closed reaction regression and existing M1B original icoFoam same-input flow comparison as applicable. This verifies only scoped MPI decomposition/reconstruction and numerics; no real environmental chemistry, open boundaries, deposition, soil PDE, air–soil coupling or independent science review. These remain NOT_RUN. Keep Draft; no science-validated/merge-verified/auto-merge without genuinely independent non-Bot current-SHA science review and all applicable checks.
