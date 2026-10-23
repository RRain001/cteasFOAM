# M1B-1d-1: Genuine three-dimensional Beltrami / ABC exact-solution verification

This increment covers only the single-rank analytic and grid-refinement
part of [Issue #11](https://github.com/RRain001/cteasFOAM/issues/11).
Unlike previous 2D-in-a-3D-mesh Taylor–Green and Poiseuille benchmarks,
**all 3 spatial dimensions contain more than one cell**, and **all 3
velocity components are nonzero and nontrivial**. This does not modify
the already accepted solver implementation.

## Independent closed-form Navier–Stokes reference

For k=2π and x,y,z∈[0,1], the ABC A=B=C=1 Beltrami field is

u0 = sin(k z) + cos(k y)
v0 = sin(k x) + cos(k z)
w0 = sin(k y) + cos(k x)

Its divergence is exactly zero, curl(U0)=k U0, Laplacian(U0)=-k² U0, and
the inertial acceleration is (U·grad)U=grad(|U|²/2). With kinematic
viscosity ν=0.01 m²/s, a genuine nontrivial nonlinear exact solution is

U(t)=U0 exp(-ν k² t)
p(t)=-|U0|²/2 exp(-2ν k² t) + arbitrary pressure gauge(t).

The domain is a unit cube periodic in **x, y and z**, with 8³ and 16³
regular finite-volume hexahedral meshes, implicit Euler dt=0.0025s,
20 steps until 0.05s, 3 PISO inner pressure correctors, no body force
and a spatially constant passive scalar c=1 with D=0.

## Independent numerical gates

- Parse real (not synthetic) OpenFOAM final U and p field files.
- Evaluate all 3 U components at every cell against the exact transient
  solution and compare p modulo its additive gauge.
- Verify actual 20-step maximum conservative face divergence
  <=1e-7 s^-1 and total six-face cyclic flux <=1e-7 m³/s.
- Verify domain-integrated constant passive scalar drift <=1e-8.
- Confirm each U component has physical dynamic range (>2.0 m/s).
- Verify finite kinetic energy in physical range 1.0..1.5 m²/s².
- Separate coarse/fine U/p error limits and require actual mesh
  refinement ratio <0.85 velocity and <0.9 pressure.

**All earlier official v2606 CI stages remain mandatory**, including
diffusion, advective refinement, 2-rank PDE reconstruction, moving-wall
Couette, Taylor–Green U/p and actual 2-rank U/p equivalence, and
no-slip Poiseuille benchmark.

## Explicitly not yet tested

Decomposed two-rank true-3D Beltrami flow, decomposition along different
axes, realistic inlet/outlet pressure and arbitrary nonorthogonal mesh,
building wakes, 3D buoyancy, LES/RANS/URANS, radioactive HT/HTO, soil
and coupled air–soil interface. These remain separate milestones, so
passing this stage must **not close Issue #7 or #11**.
