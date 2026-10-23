# M1C-0: HT/HTO activity reaction-box verification (predeclared)

Status: **CANDIDATE — NOT SCIENTIFICALLY ACCEPTED**. This gate is frozen before the first numerical run. Parent is the unmerged M1B-1g PR #20 exact head `e5187a5715006205ac8f24696f28101e07d2cf1b`; the M1C work is a separate stacked branch, not an expansion of PR #20.

## Model and scope

An optional `activityMode htHto` reads separate OpenFOAM volume fields `HT`, `HTO` with physical dimensions `[0 -3 -1 0 0 0 0]` (Bq/m³, since Bq=s⁻¹). The original dimensionless passive scalar `c` remains unchanged. Nonnegative rates `kOx` and `lambdaT` have dimensions `[0 0 -1 0 0 0 0]` (s⁻¹); nonnegative `DHT`, `DHTO` are m²/s. The finite-volume equations are

`ddt(HT) + div(phi,HT) - laplacian(DHT,HT) + (kOx+lambdaT)*HT = 0`

`ddt(HTO) + div(phi,HTO) - laplacian(DHTO,HTO) + lambdaT*HTO = kOx*HT`.

The second equation uses the **newly solved HT** so backward-Euler conversion gain equals the HT conversion loss, to within the actual FVM algebraic solver tolerance. `lambdaT` is user-specified; the accelerated numerical test value is not a physical tritium decay constant. No chemical fractionation, HTO/HT exchange with water, deposition, soil or turbulence is modeled here. Activity inventory is integrated in Bq, not m³. Negative rates or diffusivities must be rejected.

## Frozen numerical gates

Use the official signed OpenCFD v2606 Ubuntu Docker image and an actual 3-D 4×3×2 finite-volume unit-cube mesh, six closed `zeroGradient` species patches, U=0, `momentumMode prescribed`, D=0 and `DHT=DHTO=0`. At t=0 initialize uniform `HT=1e6 Bq/m³`, `HTO=2e5 Bq/m³`; take 20 steps with Δt=0.05 s to t=1.0 s. Execute four separately generated cases: `(kOx,lambdaT)` in s⁻¹ = `(0.4,0)`, `(0.4,0.2)`, `(0,0.2)`, `(0,0)`.

The *independent* Python checker derives the exact **backward-Euler recurrence**, not a falsely claimed continuous-time exact solution:

`HT[n]=HT[n-1]/(1+dt*(kOx+lambdaT))`

`HTO[n]=(HTO[n-1]+dt*kOx*HT[n])/(1+dt*lambdaT)`.

Predeclared PASS gates: all four cases run 20 actual steps, 3-D mesh has 24 actual cells, both species' final full-field maximum absolute error vs this recurrence ≤`1e-6 Bq/m³`; every per-step HT and HTO independently reconstructed budget residual ≤`1e-6 Bq` (initial inventory 1.2e6 Bq); combined activity residual ≤`2e-6 Bq`; no negative or nonfinite fields; exact t=1.0 output exists. Both rate=0 control remains unchanged within same tolerance. A deliberate +100 Bq/m³ corruption of final HT must fail the checker. Cases with negative kOx or lambdaT must terminate nonzero with the expected explicit diagnostic, not merely fail to launch. A test case with a wrongly dimensioned HT input must fail OpenFOAM's dimension check.

For this closed reaction-box case, the independently reconstructed budgets are `ΔA_HT + dt*(kOx+lambdaT)*A_HT,new=0`, `ΔA_HTO + dt*lambdaT*A_HTO,new - dt*kOx*A_HT,new=0`, `Δ(A_HT+A_HTO)+dt*lambdaT*(A_HT,new+A_HTO,new)=0` (Bq). No boundary flux exists. **No assertion is made for advective/diffusive open boundaries or realistic physical chemistry until separately tested.**

## Required regression and scientific limits

Preserve M1A/M1B existing `c` behavior; run exact-SHA Ubuntu Docker `smoke` and selected `full`, and the new targeted Docker numerical script. For the new reaction-box, no official stock OpenFOAM application is asserted to implement identical HT→HTO chemistry; upstream same-input comparison is **NOT_RUN**. Independent external radiological chemistry/field validation, MPI species budgets, spatial advection–diffusion, soil PDE and air–soil interface are **NOT_RUN**. Do not apply `science-validated` or `merge-verified` or merge without genuine independent current-SHA scientific review and all applicable checks.
