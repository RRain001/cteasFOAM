# M1B-1f-3 — True MPI of geometrically non-orthogonal 3D ABC flow

Following merged PR #16 (physical warped mesh, OpenFOAM cell-centroid
initialization and serial ABC analytic/refinement) and #17 (distinct
non-orthogonal inner pressure iterations), this adds a **real distributed
pressure/velocity PDE verification** on the same warped 8x8x8 mesh.

## Inputs and numerical settings

For each of the serial reference, x-axis 2-rank and z-axis 2-rank
runs: create OpenFOAM 8^3 unit cube with six cyclic patches, deform
interior vertices by x'=x+0.06 sin(pi*x)sin(pi*y)sin(pi*z), preserve
each external boundary point and original mesh topology. Run native
checkMesh -allGeometry -allTopology and require Mesh OK. Native
writeCellCentres must produce 0/C, from which all initial ABC U and
p values are recomputed. The actual polyMesh points, all cell centres,
initial U and initial p must be byte-identical among the three cases.

The existing laminar PISO mode runs kinematic viscosity 0.01m^2/s,
three outer pressure-velocity correctors plus **two additional inner
nonorthogonal pressure correctors** (20 timesteps, dt=0.0025s, t=0.05s).
The three runs each must report 20 nonorth settings and **180 actual p
linear-system solves**. Two real OpenMPI processes solve after OpenFOAM
decomposePar using (2,1,1) and (1,1,2), respectively, followed by
reconstructPar -latestTime. Rank existence, solver completion and
reconstructed final fields are mandatory.

## Prespecified independent quantitative acceptance

From actual OpenFOAM full 512-cell U and p output (not synthetic arrays),
compare serial with each MPI result and compare x versus z:
- max three-component U cell-wise difference <= **0.005 m/s**
- gauge-free cell-wise kinematic p difference <= **0.02 m^2/s^2**
- serial vs MPI final global integrated kinetic energy difference <= **0.003**
- each MPI maximum corrected |div(phi)| <= **1e-7 s^-1**
- sum of periodic patch fluxes <= **1e-7 m^3/s**
- passive scalar inventory drift <= **1e-8**
- exactly 20 recorded solver steps per case.

Independently evaluate initial analytical ABC fields at actual OpenFOAM
0/C for all three runs and demand <=5e-12 discrepancy; demand centroid
deformation 0.02..0.09 m against Cartesian index coordinates.

The U/p/energy limits are set before this run and are more stringent
than the earlier orthogonal ABC MPI tolerances in #13. No relaxation
without a documented independently understood failure and new audit.

All previous official v2606 CI numerical stages must pass on the same
exact PR HEAD. This PR changes only tests/CI, no numerical solver
physics. Even when passed, it does **not** validate high skewness or
>10-degree deformation sensitivity, general pressure patches on warped
meshes, real 3D buildings, larger MPI counts, buoyancy, URANS,
radiological HT/HTO or soil transport.
