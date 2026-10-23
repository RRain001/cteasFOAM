# M1B-1b — Real two-rank PISO pressure-correction equivalence

Follow-up to merged M1B-1a PR #8. The purpose is a scientifically
distinct check: actual two-process pressure-corrected Navier–Stokes,
not an MPI_Allreduce smoke or M1A scalar-diffusion-only verification.

The previously analytically checked 12×12×1 Taylor–Green case
is run serially and decomposed 2×1×1 along x using OpenFOAM decomposePar.
The same executable runs as 2 OpenMPI solver ranks with -parallel; then
reconstructPar produces an actual 144-cell U/p field at t=0.1s.

The checker compares the entire final U vector field, gauge-normalized
kinematic pressure field, pressure-corrected face-flux divergence logs
and final integrated kinetic energy. Requirements:
- maximum cellwise U difference <= 1e-3 m/s
- pressure difference after removal of each field mean <= 2e-3 m²/s²
- pressure-corrected face divergence <= 1e-7 s^-1 throughout MPI run
- net boundary volume flux <= 1e-7 m³/s
- final global kinetic energy difference <= 1e-4
- at least 20 time steps logged in serial and parallel
- all prior M1A, Couette and Taylor–Green analytic CI gates must remain PASS

This does not validate pressure boundary conditions for buildings,
3D complex geometries, higher MPI counts, PIMPLE, turbulence, HT/HTO or soil.

## Temporal-refinement diagnostic (2026-10-23)

After MPI pressure reference repair, the 20-step benchmark showed
conservative corrected flux (`max_div_phi=1.08e-13 s^-1`) but velocity
serial/MPI L-infinity difference `1.79374566371e-3 m/s`, exceeding
`1e-3 m/s` gate. To distinguish temporal splitting error from MPI
boundary coupling defects, this PR runs the *same* 0.1 s periodic vortex
with 40 implicit-Euler steps (`dt=0.0025s`) both serial and MPI.
The original analytic Taylor–Green tests at dt=0.005 s are unchanged.
No numeric tolerance is relaxed. If fine-time-step MPI does not meet
existing gates, leave PR Draft and diagnose pressure flux reconstruction.
