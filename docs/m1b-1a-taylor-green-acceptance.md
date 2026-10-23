# M1B-1a: first pressure-corrected laminar Navier–Stokes increment

## Scope and limits

The `momentumMode pisoLaminar` branch implements a first laminar
pressure-projection pass with inertial convection, viscosity, kinematic
pressure p in m²/s², a pressure Poisson equation, correction of face flux
phi and velocity, and propagation of the corrected phi into scalar
transport. Older `prescribed` and `laminarShear` modes remain selectable.

**This is the first increment of Issue #7, not completion of M1B-1.**
The only currently planned numerical reference is a two-dimensional
periodic Taylor–Green vortex on a 3-D finite-volume mesh. Arbitrary 3-D
boundary conditions, nonorthogonal pressure correction, wall-driven
flow, MPI pressure correction and turbulence are not yet certified.

The analytic Taylor–Green reference on the unit periodic domain is:

u = sin(2πx) cos(2πy) exp(-8π² νt)
v = -cos(2πx) sin(2πy) exp(-8π² νt)
p = +(cos(4πx)+cos(4πy))/4 exp(-16π² νt)

The vorticity evolves with molecular nu=0.01 m²/s. Both
12×12×1 and 24×24×1 OpenFOAM cases run 20 implicit Euler steps
of 0.005 s and are parsed by an independent analytic checker.
Initial p is analytical. Pressure errors are evaluated modulo
the arbitrary gauge constant, with cell averages subtracted.

The checker examines full U/p final fields, requires each time step to
log conservative face flux divergence and net boundary flux,
checks kinetic energy bounds and finite outputs, and compares coarse/fine
velocity and pressure errors. All prior M1A/Couette tests must still pass.

These bounds are initial non-regression scientific gates, not a formal
mesh-independent convergence or universal 3-D accuracy claim.

The pressure sign follows the chosen velocity orientation:
`u=sin(2πx)cos(2πy)`, `v=-cos(2πx)sin(2πy)`.
Its advective acceleration is `(+πsin(4πx), +πsin(4πy))`,
which is balanced by `-∇p` for `p=+(cos(4πx)+cos(4πy))/4`.
This sign was caught by the first full-field pressure regression;
retain both the failed CI and this correction in review evidence.
