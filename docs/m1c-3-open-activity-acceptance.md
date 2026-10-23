# M1C-3: open-boundary HT/HTO activity and signed conservation (predeclared)

Status **experimental / NOT independent scientific acceptance**. Parent is M1C-1 Draft PR #29 exact HEAD `b25cdc709fa709b104ff3cb32fbb08f3f802768c`. This acceptance specification is committed **before implementing or running** the new case.

## Physical equation and input contract

OpenFOAM OpenCFD v2606, same native 960-cell (8×20×6) 3-D *actually warped* duct as M1B-1g with fixed pressure xMin/xMax, four real no-slip walls. Use existing laminar PISO solver, nu=0.05 m²/s; pressure-drop cases +0.1 (forward), −0.1 (reverse), 0 (zero); 3 PISO × 3 nonorthogonal p solves per 20 steps to t=0.4 s, dt=0.02 s. This is not atmospheric URANS. Zero initial HT and HTO Bq/m³ (verified from actual native field), species diffusivities DHT=0.01, DHTO=0.02 m²/s, conversion kOx=0.4 s⁻¹, artificial numeric decay lambdaT=0.2 s⁻¹. These are not field-calibrated coefficients.

On xMin and xMax use OpenFOAM `inletOutlet` scalar boundary conditions with `inletValue` equal to (HT,HTO)=(1e6,2e5) Bq/m³ **only at the physically active inflow**; the opposite end inletValue must be zero, and outflow uses zeroGradient. Four wall species BCs zeroGradient. Pressure and U wall/pressure BCs unchanged from the existing verified duct. `fluxRequired` HT and HTO must be on; existing `auditActivityFlux` prints signed matrix fluxes Bq/s. For no flow and zero initial activity, inletOutlet with zero phi yields zeroGradient and must not generate activity merely because inletValue is nonzero.

For each species the finite-volume whole-domain balance is
`A_HT,n − A_HT,n−1 + Δt F_HT,n + Δt (kOx+lambdaT) A_HT,n = 0` and
`A_HTO,n − A_HTO,n−1 + Δt F_HTO,n + Δt lambdaT A_HTO,n − Δt kOx A_HT,n = 0`.
Total `A_HT+ A_HTO` excludes internal conversion but includes open boundary net flux and radioactive decay. Positive outward face flux F. Volumes are actual OpenFOAM native V, dimensional HT/HTO must be [0 −3 −1 0 0 0 0].

## Frozen fail-closed criteria

1. Real Docker official v2606, actual blockMesh + warp + checkMesh + native cell centres and volumes; 960 cells, 665 moved vertices, and pressure BC/no-slip physical faces identical for all cases. 20 actual steps / 180 p solves in each of the **three** flow cases, no substitute synthetic log.
2. Each species and total **all 20 per-step** equation residual absolute ≤ 1e−4 Bq (both logged residual and *independent* recomputation using inventories, printed flux, dt, rates). Each HT/HTO xMin+xMax signed flux sum must match total signed boundary net to 1e−3 Bq/s (walls insulating), and xMin/xMax pressure volume flux sums ≤ 1e−7 m³/s; corrected flow divergence ≤ 1e−7 s⁻¹.
3. Forward: xMin pressure volume flux negative ≤−0.003 m³/s, xMax positive ≥+0.003, and **both species'** xMin net flux negative less than −1000 Bq/s at final time; reverse swaps xMin and xMax with same thresholds; positive final total HT and HTO inventories in each. Each field finite and nonnegative ≥−1e−5 Bq/m³; no result may exceed 1.5e6 Bq/m³. Allow positive HT→HTO conversion inside the domain; total budget, not monotonicity of each species alone, is authoritative.
4. No-flow control: final HT and HTO concentration and both inventories ≤1e−8 in respective Bq units; volume/scalar boundary fluxes and pressure velocity ≤1e−8, despite nonzero `inletValue`. No invisible source term allowed.
5. Independent native **final-volume inventory** recomputation must agree with the last reported `activity_HT_Bq`, `activity_HTO_Bq` to ≤1e−4 Bq. Full 20-step time history required, not only last summary.
6. Negative controls intentionally alter one real time-series HT net flux by +100000 Bq/s and one final native HT concentration by +100000 Bq/m³; both must cause independent checker to reject. Malformed active inlet BC / missing field should also fail.
7. Existing M1B-1g-12 / M1C-0 / M1C-1 regressions must remain PASS. Run exact SHA targeted official v2606 Docker separately from managed `smoke/full`, and run exact SHA managed `smoke/full` before publishing.

The built-in `htEqn.flux()` records use the *solved* equation's discretization and are not fully independent face-flux measurement. Independent field/experimental validation, atmospheric URANS, real decay parameters, wet/dry deposition, soil and nonmatching air–soil interface, external production applicability, and independent human review are **NOT_RUN**. Keep this as a Draft research PR and never bypass scientific admission.
