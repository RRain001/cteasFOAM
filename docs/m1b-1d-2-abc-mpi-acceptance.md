# M1B-1d-2: real two-rank 3-D ABC/Beltrami pressure PDE equivalence

Follows independently accepted genuine 3D single-rank ABC analytical
and 8³→16³ grid refinement, PR #12. This PR adds **no changes to the
laminar solver source**: it is an independent distributed-pressure scientific
gate for Issue #11. Previously, only the extruded *2D* Taylor–Green flow
was compared between serial and MPI.

## Cases

Identical 8×8×8 finite-volume periodic ABC initial U and pressure,
nu=0.01m²/s, 20 time steps dt=0.0025 through t_end=0.05s,
3 PISO pressure correctors, no body force, scalar c=1, D=0.

- A genuine serial 512-cell reference solver execution;
- decomposePar n=(2 1 1) and genuine 2-rank mpirun solver, reconstructPar;
- decomposePar n=(1 1 2) and genuine 2-rank mpirun solver, reconstructPar.

The **separate independent checker** parses each full reconstructed
0.05/U and 0.05/p field (512 real finite-volume values each) and all
three computations' actual 20-step conservative face flux/inventory logs.
It compares all U components cell by cell, gauge-normalized p fields,
max corrected |div(phi)|, all six periodic boundary face net flux,
scalar conservation, and global final kinetic energy.

Explicit predeclared gates: serial-vs-MPI and x-vs-z L-infinity U
discrepancy ≤0.015m/s; pressure modulo gauge ≤0.05m²/s²;
per-case corrected |div(phi)| ≤1e-7s^-1;
net boundary flux ≤1e-7m³/s;
passive scalar integral drift ≤1e-8;
final total kinetic energy difference ≤0.01.
All 20 global PISO/transport records are mandatory. Failed real solver,
mesh decomposition, rank count and reconstruction are fatal.

All earlier M1A scalar physics/MPI, M1B Couette, 2D Taylor–Green,
no-slip Poiseuille and 3D ABC analytic/grid-refinement benchmarks
are required on **the exact same PR HEAD** in official OpenFOAM v2606 CI.

## Limits

This verifies a **specific two-rank true-3D periodic laminar flow**,
not speedup, larger MPI topologies, arbitrary pressure gauge-reference
placement, walls/inflow/outflow, nonorthogonal meshes, buoyancy, URANS,
HT/HTO reaction transport, soil, or conservative air–soil interfaces.
Keep parent Issue #7 open even after these gates pass.
