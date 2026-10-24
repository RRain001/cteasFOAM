# M1C-10 native x-only refinement findings

This report records OpenCFD v2606 *development* numerical observations following the separately committed frozen test contract `docs/m1c-10-x-refinement-source-quadrature.md`. These values were obtained with real native `wmake` C++ 3-D warped pressure-driven PISO, real `blockMesh`/`checkMesh` and real final HT/HTO fields. A subsequent exact-SHA rerun is required for official evidence; do not quote the initial development run as an independent science acceptance.

## Frozen identical physics

One-meter cubical pressure-open duct and *synthetic* Gaussian volumetric HT/HTO source only, not stack and not real-scale atmospheric dispersion. Source (0.38,0.62,0.55)m, sigma=.13m, [0.075,0.225)s, total HT 9000 Bq and HTO 1800 Bq, dt .02s, final 0.4s. **Y=20, Z=6 never vary**, X-only counts 8,12,16,24,32; strict same R3=(.88,.62,.55)m, receiver Gaussian averaging width .10m.

| NX | actual cells | R3 HT Bq/m³ | R3 HTO Bq/m³ | X-source centroid error (m) | X-source second-moment error (m²) |
|---:|---:|---:|---:|---:|---:|
| 8 | 960 | 3217.294052 | 1296.781612 | -3.64417e-4 | 6.25178e-5 |
| 12 | 1440 | 2950.025867 | 1222.241199 | -2.22015e-4 | 2.65994e-5 |
| 16 | 1920 | 2828.800637 | 1187.159237 | -1.67539e-4 | 1.34343e-5 |
| 24 | 2880 | 2718.664316 | 1154.314912 | -1.27249e-4 | 3.91456e-6 |
| 32 | 3840 | 2668.282794 | 1138.894367 | -1.12886e-4 | 5.65289e-7 |

R3 HT neighbor changes at 8→12,12→16,16→24,24→32 are -267.268, -121.225, -110.136, -50.382 Bq/m³. R3 HTO changes -74.540, -35.082, -32.844, -15.421 Bq/m³. Both satisfy frozen `TREND_DECREASED` test for the last three differences. However, relative final-step receptor sensitivity remains **1.89% HT** and **1.35% HTO** (using the NX=32 value as denominator). This is **not a grid-independence statement** or numerical-order estimate, because the physical 3-D cell aspect ratios and unrefined Y/Z errors remain, and there is no analytic independent plume concentration.

Analytic source-shape reference computed with `erf` and exact first/second Gaussian moments over the actual unit cube. Fixed Y/Z resolutions mean the **total 3-D normalized Gaussian quadrature bias does not go to zero** with X-only refinement: continuous-kernel integral relative error decreases from about 0.002275 at NX=8 to 0.001272 at NX=32, while Y/Z centroid/moment errors remain nonzero and may change non-monotonically. The source total Bq remains exactly 9000/1800 because the native source itself is **renormalized**; this does not prove its volumetric shape is exact.

No physical/source coefficient was changed. M1C-8/9 numerical/source tests retain their thresholds and reference results. The only legacy-checker code correction replaces mathematically invalid coarse-lattice inference `round(warped_x * nx)` with inversion of the same exact monotone interior deformation formula before asserting the **unchanged 2e-10m per-vertex geometric tolerance**. With NX>=16, interior shifts can exceed half a nominal x-cell and defeat naive rounding. This is **stricter correct geometry auditing, not loosening**.

## Scope and pending science

This numerical study supports a *limited* x-direction refinement trend at synthetic R3 and independently measures spatial Gaussian source quadrature error. Causality of downstream concentration variation is not resolved: numerical convection/diffusion, source discretization, mesh skewness and boundary effects all change with x refinement. No URANS, building/stack, realistic air boundary, soil PDE, interface, human receptor/dose or independently validated plume fields are present. External science and non-Bot independent acceptance **NOT_RUN**.
