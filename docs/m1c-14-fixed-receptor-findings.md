# M1C-14 findings: fixed-coordinate 3D Gaussian plume observations

Acceptance frozen in commit 846dab35ecce374ee7c2555cc22c7b892b9298d7. Initial 8³ native run scientifically FAILed because the draft observer assumed all warped OpenFOAM volume centroids had exact unwarped y/z centre coordinates. Geometry erratum committed before implementing corrective trilinear physical inversion: 0aff358. The final independent 8-native-centre interpolation uses actual physical x/y/z centroids, Newton inversion of a full nonorthogonal trilinear hexahedron, eight nonnegative weights summing to 1, and exact physical affine reproduction in 1e-9m. Grouping into logical structured columns uses a separately checked maximum 0.01m physical y/z deviation; this only locates stencil indices, never substitutes nominal positions for actual centres. No C++ source changes.

Four genuine OpenCFD v2606 warped cyclic 3-D Gaussian HT/HTO cases (8³, 12³, 16³, 24³) at dt=.005s, 80 actual solved steps, initial zero HT/HTO, fixed uniform U=(1,.05,.03)m/s, Gaussian release centre (.42,.55,.5)m sigma=.08m, HT 60000 and HTO 12000 Bq/s on [.073,.137]s; exact 3840/768 Bq injected in all four cases. Signed max species/total Bq residual 1.3573e-8 Bq on 24³, lower on coarser meshes. Native mesh, source physics, 6 cyclic boundary patches, true C/V and U, HT, HTO fields and log inventories all audited, plus 3 intended native-corruption controls.

**Fixed physical coordinates, true 3-D nonorthogonal interpolation, Green semi-analytic at SAME fixed physical coordinates**:

| N³ | R1 HT | R1 HTO | R2 HT | R2 HTO | R3 HT | R3 HTO |
|---|---:|---:|---:|---:|---:|---:|
| 8³ | 33.66% | 30.42% | 62.24% | 54.37% | 19.47% | 24.19% |
| 12³ | 19.23% | 19.25% | 48.98% | 41.97% | 0.41% | 9.54% |
| 16³ | 16.67% | 17.56% | 41.43% | 35.03% | 2.13% | 7.44% |
| 24³ | 9.35% | 11.52% | 33.07% | 27.50% | 8.64% | 1.66% |

R1/R2 both species and R3 HTO show monotone error decrease, but R3 HT minimum at 12³ followed by 16³ and 24³ increases. No all-receptor monotone convergence, no global mesh independence. The genuine per-mesh closest native centre moves; M1C-14 removes that geometric drift as the *probe evaluation position*, not as a source of mesh error. 24³ R2 fixed HT native 78975.698 versus analytic 117999.108 Bq/m³; signed error -39023.410 Bq/m³. The exact additive decomposition at R2 is ~-34504.511 Bq/m³ from the **interpolated native field minus sampled Green analytic field at the same eight anchor centres** (a combined FVM/source/reactive error, not uniquely source), plus ~-4518.899 Bq/m³ from interpolation of the analytic Green field vs Green at the precise R2 location. These close numerically to <1e-7 Bq/m³. The interpolation component is comparatively small; high R2 discrepancy persists even at 24³.

The source-shape limitation inherited from M1C-13 remains: continuous periodic *wrapped* full-Gaussian Green reference versus cell-centre-sampled and renormalized native Gaussian. Analytic Gaussian mass outside the unit box is 8.573876375539413e-8. No reference to PISO/URANS field flow, any physical stack, soil PDE or air-soil exchange. This is a scoped synthetic-geometry mathematical verification, NOT complete field model validation.

M1C-14 initial failure and documented geometry correction are preserved. Developer native cases and strict negative controls passed after that correction; the exact pinned SHA official Docker and Ubuntu full acceptance are recorded independently after commit. External science, independent non-Bot review and field validation remain NOT_RUN; stacked PR remains Draft and unmerged.
