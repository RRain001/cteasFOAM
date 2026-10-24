# M1C-4: two-rank warped 3-D pressure-open HT/HTO parallel conservation (predeclared)

Status: PROPOSED / NOT SCIENTIFICALLY ACCEPTED. Branch based on Draft PR #30 open-boundary solver tests, including the correction from Draft PR #31. This gate must be committed *before* any M1C-4 numerical runs. This is an additional scoped test, not evidence for atmospheric URANS, deposition, soil PDE or nonmatching interfaces.

## Identical numerical setup

Use official OpenCFD v2606 Jammy Docker and the existing M1C-3 warped duct, 960 physical cells, 665 moved vertices, true 3-D no-slip y/z walls, fixed p at xMin/xMax. Forward/reverse/zero driving drops +0.1/-0.1/0; 3 outer PISO and 2 nonorthogonal inner correction passes (9 pressure solves per step), 20 Euler steps dt=0.02 s to 0.4 s. Nonzero inlet conditions for HT and HTO only at actual inflow side: 1e6 and 2e5 Bq/m3 respectively, zero initial activity. DHT 0.01 m2/s, DHTO 0.02 m2/s, oxidation kOx 0.4 s^-1, artificial decay 0.2 s^-1. Keep original corrected Laplacian, gradient, convection and solution tolerance.

For each mode perform exactly one serial solve and a real two-rank x=(2,1,1) decomposition via decomposePar and mpirun -np 2, followed by reconstructPar -latestTime. Both must start from identical real mesh/fields and use identical solver code. Run the pre-existing M1C-3 native field/volume/mesh/checkMesh and signed-balance checker on **both** serial and MPI cases.

## Frozen acceptance

1. All 6 solver runs exit 0, 20 time steps, 180 genuine p solves, 20 HT and 20 HTO solves and logs. Both processor meshes and final rank-specific fields exist. Native reconstructed HT/HTO/U fields include exactly 960 cells and retain Bq/m3 units for species.
2. For both serial and MPI, existing M1C-3 independent checker passes all 20 per-step HT, HTO and total activity budgets: absolute residual <= 1e-4 Bq; physical signed inlet fluxes in forward and reversed cases must have expected direction and magnitude, and zero-pressure nonzero inletValue must not inject any activity.
3. Additional direct serial vs MPI parity: final HT pointwise |Δ|_max <= 100 Bq/m3 and RMS <= 30 Bq/m3; HTO pointwise |Δ|_max <= 20 Bq/m3 and RMS <= 6 Bq/m3 (both ≤0.01% of their respective 1e6/2e5 inlet concentrations); final U vector maximum componentwise difference <= 0.003 m/s. For **every one** of the 20 steps in forward and reverse cases, each species inventory absolute difference <= 20 Bq and each signed xMin/xMax/net species flux absolute difference <= 25 Bq/s. Signed physical volume xMin/xMax flux absolute difference <= 1e-6 m3/s at every step. For zero case preserve existing <=1e-8 near-zero activity/velocity/flux requirement.
4. Inject +100000 Bq/m3 into one reconstructed MPI HT cell after successful checks, then require the new independent parity checker to reject it. A negative mode or malformed result must fail, never silently skip.
5. No changes to numerical physics, case input, or thresholds in response to failure; failed x-decomposition remains FAIL until a physical cause is fixed. Run M1C-0, M1C-1, M1C-2 and M1C-3 scoped regressions in the same candidate, then managed exact-commit Docker full, GitHub repository-admission and independent scientific-admission. Do not merge Draft absent genuinely independent review.

These checks establish only a restricted laminar pressure-open HT/HTO integration case; external HT chemistry and environmental fields, URANS, dry/wet deposition, soil physics and air-soil coupling remain NOT_RUN.
