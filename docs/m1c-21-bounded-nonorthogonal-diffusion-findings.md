# M1C-21 — bounded conservative nonorthogonal diffusion correction findings

Frozen prior to implementation: acceptance commit 6512f93c5e8dcef2ba9836c62629723f3ccde841. Stacked on M1C-20 Draft #48, parent SHA f2c879d8c42d45caf81c2f27ce1315df989b56c3. All numerical outcomes below are scoped to prescribed-velocity synthetic native OpenFOAM v2606 cases and do not establish external radiation science.

## Numerical method

M1C-20 found genuine low-HT negativity (-2.75241736514867e-6 Bq/m3) on a warped grid at dt=0.0025 s with standard nonorthogonally corrected implicit diffusion and the existing advective face MULES. That original counterexample and default program are UNCHANGED and remain required as a negative-control regression.

M1C-21 explicitly opts into activityDiffusionCorrection faceMules, which requires activityConvectionBlend faceMules and HT/HTO mode. The implicit low-order FV matrices use Gauss linear uncorrected two-point diffusive fluxes, each from an explicit named scheme. On the actual solved low state, define high-minus-low *signed diffusive face flux*:

    Fdiff,corr = - D * magSf * (snGrad_corrected(C_low) - snGrad_uncorrected(C_low)).

Add this face correction to the native high-minus-upwind *advective* correction. One synchronized per-face OpenFOAM MULES limiter acts on their SUM; cell activity is updated with -dt*div(lambda_face*Fcorrection). The same actual limited signed flux is added to the frozen low matrix face flux for separate HT, HTO and combined tritium budgets. Oxidation into HTO refers only to the already corrected HT. Noncoupled wall/physical inlet/outlet correction flux remains zero; their native implicit low physical boundary flux is retained. Cyclic and processor-coupled interfaces use native per-face synchronization.

Require explicit named Gauss linear uncorrected positive low Laplacian and corrected/uncorrected native high/low snGrad schemes, verified before ANY timestep. Wrong or missing scheme and incompatible configurations fail closed; no fallback to scheme dictionary defaults. No clipping of cell activity, source retuning or total-activity renormalization. There is no change to the old M1C-18/M1C-19/M1C-20 opt-out branches. The experimental explicit limited correction is not identical to the original unlimited nonlinear high-order corrected implicit solve, nor a proof of universal unconditional positivity.

## Actual native development tests, official v2606

True warped 3D cyclic Gaussian source HT 3840 Bq, HTO 768 Bq, prescribed nonzero U, diffusion DHT=.01 and DHTO=.02, oxidation kOx=.4, decay .2; every-step native nonnegative finite cells and signed species/total Bq <=1e-5, native Green R1/R2/R3 physical receptors, genuine corrected face-amplitude logs:

| Grid | dt (s) | Native steps | Worst signed Bq residual | Correction active, each isotope |
|---|---:|---:|---:|---|
| 8³ | 0.0025 | 160 | 1.042365394e-8 | 131/160 |
| 8³ | 0.005 | 80 | 8.9875697995e-9 | 66/80 |
| 8³ | 0.01 | 40 | 8.14832112894e-9 | 33/40 |
| 16³ | 0.005 | 80 | 7.04608016378e-9 | 66/80 |

Genuine 8³ x and y two-rank and xy four-rank native decompositions, real 80-step signed Bq and reconstructed HT/HTO fields, all passed unchanged serial-vs-MPI 3e-4 Bq/m3 and 1e-5 Bq tolerances. Authentic 960-cell warped nonperiodic xMin inlet and reversed xMax inlet with physical inletOutlet/wall patches: 20-step native HT/HTO/total signed physical boundary budgets and nonnegative primitive fields passed, with nonzero internal nonorthogonal corrected flux.

Seven genuine fail-closed controls also passed: wrong and missing positive low scheme, missing and wrong high snGrad, invalid option, incompatible options and independently forged negative correction-amplitude log. All tested using the *same compiled real OpenFOAM solver* and original genuine meshes; no fake test stubs.

The new scheme resolves the specific fine-clock negative low state without changing original acceptance thresholds. It may alter spatial truncation error, since the nonorthogonal correction is now explicit and per-face bounded rather than original unbounded implicit correction. General convergence rates, larger mesh anisotropy, dynamically updated wind/WRF, genuinely independent field data, URANS, soil HT/HTO PDE and nonmatching air-soil conservation coupling have NOT been demonstrated. Do not mark these PASS.

## Admission discipline

Same-source Ubuntu Docker full 7 steps, previous M1C-0..20 21 genuine native targeted suites, and new M1C-21 runner must all be verified at the exact candidate SHA; local development result alone is insufficient. Formal independent non-Bot same-SHA scientific review and external field validation remain NOT_RUN. Any new PR stays Draft, with scientific-admission blocking; never apply science-validated or merge-verified or merge it without those independent gates.
