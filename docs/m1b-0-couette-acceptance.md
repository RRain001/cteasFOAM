# M1B-0 — Restrained laminar momentum benchmark

The existing atmospheric executable now supports two explicitly separate modes
through constant/transportProperties:

- momentumMode prescribed (default) — preserve all validated M1A passive-scalar
  diffusion, cyclic advection, and real two-rank MPI cases without requiring nu.
- momentumMode laminarShear — solve the unsteady kinematic-viscosity equation
  dU/dt - nu Laplacian(U) = 0 while updating phi, on a case whose velocity is
  divergence-free by construction.

For the latter, the currently accepted special-case geometry is planar
Couette: U=(u(y,t),0,0), lower plate at rest, upper plate speed 1 m/s, gap
L=1m, kinematic viscosity nu=0.05 m²/s, zero initial interior velocity.
The nonlinear convective acceleration and pressure gradient vanish in this
specific geometry only. It does **not** supply incompressibility projection
or pressure-velocity coupling for arbitrary three-dimensional flow.

The independent Fourier-series reference is:

u(y,t) = y + sum[n=1..infinity] 2*(-1)^n/(n*pi)*sin(n*pi*y)*
exp(-0.05*(n*pi)^2*t)

The 40-cell unsteady OpenFOAM case has dt=0.005s, t_end=1s.
The checker reads actual output vector fields (not internal C++ values) and
requires maximum velocity error <=0.003 m/s, negligible transverse velocity,
monotonic shear between 0 and 1m/s, and 200 actual momentum time-step records.
Previously passed M1A and MPI numerical tests must also pass in the same CI.

Remaining M1B stages: proper 3D transient incompressible momentum/pressure
coupling and divergence budget, then validated turbulence closures and
atmospheric forcing. No URANS, HT/HTO, soil or interface acceptance claim here.
