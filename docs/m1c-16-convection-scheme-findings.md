# M1C-16 findings: upwind numerical diffusion versus limiter overshoot on actual 3D native fields

The scientific acceptance was frozen **before changing C++ or running this experiment** (commit `31e29b5708410e0217dceff626e235d93d63e069`). `docs/m1c-16-positivity-failure-addendum.md` records the genuine first OpenFOAM limitedLinear failure and diagnostic classification protocol (commit `8fd64d2`). This is a *negative scientific result for the bounded candidate* even when targeted diagnostic software exits zero. No HT/HTO physical transport/reaction/source terms were changed. The only C++ addition is opt-in `auditActivityExtrema true`, printing exact global native internal-cell minima/maxima each of 80 steps. It does NOT clip or alter negative values. Previous source sampling remains default `cellCentre`; previous actual transport choice remains `Gauss upwind`.

## Native problem and main result

Six actual OpenCFD v2606 genuine 3D warped cyclic polyMesh CFD scalar solves: paired 8³, 16³ and 24³ real meshes, `Gauss upwind` versus `Gauss limitedLinear 1` for **HT and HTO only**. Passive scalar c, actual source kernel/shape, mesh, time, wind and reaction identical. Prescribed uniform U=(1,.05,.03)m/s; HT 60000 Bq/s and HTO 12000 Bq/s Gaussian source centered (.42,.55,.50)m sigma=.08m, active [.073,.137]s with exact cumulative emissions 3840/768 Bq. kOx=.4/s, synthetic lambda=.2/s, D_HT=.01, D_HTO=.02m²/s, dt=.005s, final .4s (80 actual steps). The independent M1C-14 fixed-physical-point nonorthogonal trilinear reconstruction and pure-Python periodic Gaussian Green reference are unchanged. All six native signed HT/HTO/total Bq budgets satisfy 1e-5 Bq, but **limitedLinear violates pointwise non-negativity at intermediate times**. Integrated activity conservation does not imply positive cell concentrations.

| Mesh | R2 HT upwind relative error | R2 HT limitedLinear1 error (UNSAFE) | R2 HTO upwind error | R2 HTO limitedLinear1 error (UNSAFE) |
|---|---:|---:|---:|---:|
| 8³ | 62.24% | 51.05% | 54.37% | 41.51% |
| 16³ | 41.43% | 19.23% | 35.03% | 13.74% |
| 24³ | 33.07% | 10.95% | 27.50% | 7.67% |

All three meshes improve R2 HT and HTO relative concentration accuracy with the high-resolution flux limiter, as do most other R1/R3 fixed-position comparisons. This provides strong scoped **diagnostic evidence that excess first-order upwind numerical diffusion contributes substantially to the R2 underprediction**. It does NOT prove the remainder of bias is all source quadrature, or that the limiter is suitable for use. The improvement was measured in a non-admissible high-resolution numerical solution.

**Genuine global all-time-step native minimum across 80 solved steps:**

| Mesh | `upwind` HT / HTO minimum Bq/m³ | `limitedLinear 1` HT minimum Bq/m³ | `limitedLinear 1` HTO minimum Bq/m³ |
|---|---:|---:|---:|
| 8³ | 0 / 0 | **−2.523170723** | **−0.4387444011** |
| 16³ | 0 / 0 | **−0.0008501524338** | **−0.00006909870895** |
| 24³ | 0 / 0 | **−0.001002267809** | 0 |

Frozen positivity limit is −1e−8 Bq/m³. Therefore **limitedLinear1 scientific status: REJECTED / FAIL on all three grids**. Some final T=.4s fields are positive while intermediate HT/HTO levels were negative, confirming that checking final data alone would miss the violation. C++ opt-in extrema logs and strict Python verifier capture real intermediate failures.

Native negative controls all reject deliberately modified actual HT final cell field, a negative extrema log, wrong FV HT scheme, and a tampered emitted-HT Bq log. The 8³ *unmodified* limitedLinear strict physical checker explicitly exits nonzero with the expected genuine positivity failure. The separate diagnostic script carries on and returns software exit 0 ONLY because it explicitly checks, classifies and prints `M1C16_SCIENTIFIC_STATUS limitedLinear1=REJECTED upwind=ADMISSIBLE` and reports the true R2 change; **this is not a test granting high-order physical acceptance**. An initial C++ compile attempt failed because OpenFOAM v2606 `Foam::messageStream` does not have `precision()`; it was corrected using a local standard `std::ostringstream` with 16 digits so no existing logger/physical controls were modified. The first scripted native experiment therefore failed before the scientific addendum, and correction details are preserved in logs.

**Decision:** retain `Gauss upwind` as default; `Gauss limitedLinear 1` is only a rejected diagnostic counterfactual, not a recommended setting. The next substantive physical/numerical task should investigate a **provably positivity-preserving bounded transport method** (e.g., conservative FCT/MULES-based correction with signed HT/HTO activity accounting) and evaluate it under the same frozen Gaussian/Green protocol. Do not simply clip negatives (breaks tritium conservation). Real atmosphere URANS, stack thermodynamics, boundary-layer driving, soil PDE, air–soil exchange and independent external field science remain NOT IMPLEMENTED/NOT_RUN.

Software exact SHA and Ubuntu managed `full` outcomes are reported separately in the stacked Draft PR; no scientific-admission upgrade or independent non-Bot review is claimed.
