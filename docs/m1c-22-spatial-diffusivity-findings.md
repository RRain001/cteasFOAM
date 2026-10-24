# M1C-22 — real 3D spatially varying HT/HTO diffusion coefficients

Frozen prior to implementation at commit 86dc481db0f0d106b9bb03189d1a5ca364e0af8b. Parent M1C-21 Draft PR #49 head 2e9fdd80eeb0f381434b2e332aa0be9b42b495f4. This is a **research-only numerical generalization** under prescribed velocity, not production atmospheric science.

## Numerical change

Opt-in `activityDiffusivityMode field` is permitted only with HT/HTO and the existing bounded nonorthogonal `faceMules` diffusion/convection option. Actual user-provided OpenFOAM volume scalar inputs `0/DHTSpatial`, `0/DHTOSpatial` have dimensions m²/s. The C++ solver reads each independently, checks dimensionality and positivity/finite values of all actual local cells and boundary patches, and checks MPI-global extrema before advancing. There is no hardcoded synthetic spatial shape or fallback to uniform D.

The same per-species spatial coefficient enters the true positive implicit low-order finite-volume Laplacian (`fvm::laplacian(volScalarField, species, namedLow)`) AND the high-minus-low signed nonorthogonal diffusive face correction (`-fvc::interpolate(Dfield)*magSf*(snGrad_high-snGrad_low)`). The combined advective/diffusive face correction is limited by one synchronized native OpenFOAM v2606 MULES limiter and used exactly for the species concentration update and signed Bq budget. HT to HTO oxidation, default scalar-D cases and prior M1C-0..21 experiments stay unchanged. Source is not altered or normalized to force agreement. Unsupported options, absent field files and actual negative coefficients fail closed.

## Real native original warped periodic Gaussian source results

Cases use true 8³ and 16³ nonorthogonal 3-D mesh, 80 Euler dt=.005-s steps, native 3-D prescribed U, HT=3840 Bq/HTO=768 Bq source window, DHT reference .01 m²/s, DHTO .02 m²/s. Independent synthetic spatial input uses the genuine cell-centre coordinates: D(x,z) = D_ref[1+0.4 sin(2πx)cos(2πz)] m²/s, positive and spatially nonconstant. The analysis separately checks all actual step-by-step native HT/HTO/total signed Bq budgets, finite nonnegative real primitive scalar fields, and genuine integral agreement between volume-weighted final physical scalar arrays and model Bq logs. **The analytic Gaussian Green solution with constant D is not used as variable-D accuracy evidence.**

| Mesh | Mode | Max absolute signed Bq budget residual |
|---|---|---:|
| 8³ | Previous scalar baseline | 8.99e-9 Bq |
| 8³ | Separate genuine *uniform* D volume inputs | 8.99e-9 Bq |
| 8³ | Genuine *spatial* D volume inputs | 1.06e-8 Bq |
| 16³ | Previous scalar baseline | 7.05e-9 Bq |
| 16³ | Separate genuine *uniform* D volume inputs | 7.05e-9 Bq |
| 16³ | Genuine *spatial* D volume inputs | 1.26e-8 Bq |

Uniform-field input matches scalar baseline HT/HTO native real-cell fields with maximum deviations ≤1.02e-10 Bq/m³ on both grids and same 80-step Bq inventory within frozen 1e-5 Bq tolerance. A true variable field changes HT by >3968 Bq/m³ and HTO by >1595 Bq/m³ on the 8³ case, independently proving solver is using the user-supplied field rather than the scalar fallback.

A genuine OpenFOAM 2-rank x-decomposed warped 8³ nonuniform-D case passed native 80-step HT/HTO conserved activity and reconstructed true-cell serial-vs-MPI parity using existing 3e-4 Bq/m³ / 1e-5 Bq strict thresholds. Authentic native error controls passed for absent DHTSpatial file, truly negative DHTSpatial cell, invalid mode and incompatible opt-in modes, with no time step advanced.

## Boundaries and admission

This does NOT yet validate variable-D *physical open-boundary* cases, dynamic wind, stronger skew / irregular 3-D meshes, spatially/temporally variable tensor or anisotropic diffusion, turbulent URANS, radiation field observations, or physical air–soil exchange. They are **NOT_RUN**, not PASS. The default constant-D scientific tests (M1C-0..21) remain mandatory and unchanged. Same-HEAD Ubuntu Docker full, all previous 22 targeted tests plus this new focused spatial-D suite, and repository admission are software-only gates. Independent same-SHA non-Bot science review and external source/soil field validation remain NOT_RUN; Draft PR stays unmerged and scientific-admission blocking.
