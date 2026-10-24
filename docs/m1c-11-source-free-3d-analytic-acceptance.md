# M1C-11 pre-implementation acceptance: source-free fully 3-D analytic eigenmode

Frozen *before* implementing or running M1C-11 native cases. Parent stacked Draft PR #38 at SHA c4737ee0fdf4e5d12455ab14c7571087accda488. This is a 3-D prescribed-flow advective-diffusive HT/HTO reaction test with **no source**, to independently measure transport error without spatial Gaussian source quadrature error. No C++ solver change, full-field grid independence or real-world science validation claim.

## PDE and exact continuum reference

True unit-cube three-dimensional periodic mesh: 3 cyclic patch pairs x/y/z, real interior warped mesh (x displacement 0.04*sin(pi*x)*sin(pi*y)*sin(pi*z), zero at all exterior faces). Prescribed constant divergence-free U, two modes: ADVECT U=(0.20,0.10,0.07) m/s; DIFFUSE U=(0,0,0) m/s. No configured Bq/s source of either species and no source shape/window. kOx=.4 s^-1, synthetic lambdaT=.2 s^-1, DHT=.01 m²/s, DHTO=.02 m²/s. Set initial fields from genuine cell centres:
- F(x,y,z)=cos(2πx)cos(2πy)cos(2πz), HT=1e6+2e5 F, HTO=2e5+5e4 F [Bq/m³].
- t>0, q²=3*(2π)², F(t)=cos(2π(x-Ux*t))*cos(2π(y-Uy*t))*cos(2π(z-Uz*t)).
- a=DHT*q²+kOx+lambdaT; b=DHTO*q²+lambdaT.
- HT mean=1e6 exp(-(kOx+lambdaT)t), HT mode=2e5 exp(-a*t).
- HTO mean=exp(-lambdaT*t)*(2e5+1e6*(1-exp(-kOx*t))); HTO mode=5e4 exp(-b*t)+kOx*2e5*(exp(-a*t)-exp(-b*t))/(b-a).

The above continuum field is independent of any source kernel and must be evaluated with independent Python math, not borrowed from native CFD matrix.

## Predeclared numerics

For each mode, meshes 8x8x8 (512 cells), 12x12x12 (1728), 16x16x16 (4096). Native blockMesh, interior point warp, checkMesh -allGeometry -allTopology, postProcess writeCellCentres/writeCellVolumes, actual OpenCFD v2606 wmake C++ solver. dt=.005s, endTime=.2s, 40 real implicit Euler steps and 40 actual HT/HTO and passive c solves; no pressure or velocity equation solve. Upwind Gauss convective, corrected Gauss linear diffusion; native solver tolerance 1e-12. Calculate volume-weighted RMS and max absolute errors for HT and HTO at actual native final cell centres against the continuum closed form. All four adjacent grid RMS ratios (each species/mode, 8→12 and 12→16) strictly <.94. Native per-species max error [Bq/m³] thresholds: 8³ HT<80000 HTO<30000; 12³ HT<60000 HTO<25000; 16³ HT<50000 HTO<20000. Record actual values, never hide failures. Verify final spatially nonconstant genuinely 3-D fields, positive and finite. Actual native 0/U and last time U correspond to the fixed mode, with all three nonzero in ADVECT and exactly zero in DIFFUSE; validate no pressure solve.

Independent native source-free geometry and physics: verify every periodic physical patch and field, actual vertex positions inverted from the exact warp, native volume integral=1m³. For each of 40 steps and each species/total, independently reconstruct activity Bq balance from logged HT/HTO inventories and signed flux: HT Δ+dt*netFlux+dt*(kOx+lambdaT)*newHT; HTO Δ+dt*netFlux+dt*lambdaT*newHTO-dt*kOx*newHT; combined analogue. Max residual <1e-5 Bq, net flux for each periodic species <1e-5 Bq/s, actual final native field inventory matches logged final <1e-5 Bq. Require distinct advective and stationary fields on the same resolution. Zero configured emission mandatory, not guessed from inventory.

## Fail closed negative controls, regression, science boundary

1. Corrupt one actual 4096-cell native final HT field by +1e5 Bq/m³; independent final-field+budget verifier must reject. 2. Tamper a logged HT net face flux by +1e5 Bq/s; stepwise balance check must reject. 3. Alter one genuine prescribed U component in an isolated copy and verify native input/configuration audit rejects. Negatives must fail for intended scientific reasons (not a missing file or exception). All 11 old M1C-0..10 official native Docker suites are unmodified and must still pass plus new M1C-11 (12/12). Candidate clean exact SHA Ubuntu 7/7 software full. Formal external independent full-science validation NOT_RUN; stacked PR Draft unmerged.

Claim at most SOURCE_FREE_3D_ANALYTIC_VERIFICATION_SCOPED. Even a passing 3-D modal benchmark does not validate a localized plume, variable meteorological flow, physical stack, URANS, soil PDE, air–soil coupling or dose.
