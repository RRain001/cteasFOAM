# Architecture contract (planned; no solver validated)

## Physical domains

| Domain | Spatial dimensionality | Planned solver |
| --- | --- | --- |
| Atmosphere Ω_air | 3-D | independent OpenFOAM v2606 URANS / HT / HTO application |
| Air–soil interface Γ | 2-D surface with physical area | conservative mapping and bidirectional exchange |
| Porous soil Ω_soil | 3-D | independent OpenFOAM v2606 porous reactive-transport application |
| Plant/food/dose | kinetic/compartment modules | future external models |

The two applications use the same **OpenFOAM C++ FVM/MPI libraries**, not the same mesh, process, RAM, equations or timestep.

## Activity exchange contract

For component i ∈ {HT, HTO}, interface flux density J_i [Bq m⁻² s⁻¹] is positive **air → soil**. Over a coupling window:

  ΔA_i = integral(t_n,t_n+1) integral(Γ) J_i(x,t) dS dt  [Bq]

Use a **single authoritative interface exchange** per component/window:

  ΔA_air,Γ,i = −ΔA_i
  ΔA_soil,Γ,i = +ΔA_i

Nonmatching face meshes need *integrated activity-conservative* mapping, not naive concentration interpolation. Chemistry/conversion/decay are separate internal budgets. For implicit subiterations, rollback/restart both domain states at the same window start so that previously trialled transfer is not applied twice.

## Development sequence

M0 Docker/compiler scaffolding; M1 atmospheric independent verified scalar; M2 porous soil independently verified PDE; M3 equal/opposite 2-domain pure diffusion test with differing meshes/timesteps; M4 chemistry, deposition, field validation; M5 HPC and plant kinetics. No M1–M5 physics is currently implemented.

Reference papers/data are handled in RRain001/Academic-Experiments with rights checks.
