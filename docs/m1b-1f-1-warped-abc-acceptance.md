# M1B-1f-1: genuine warped 3D ABC analytic benchmark

First **serial geometry and physics** increment for Issue #15. The existing
solver remains untouched. The test builds actual 8x8x8 and 16x16x16 OpenFOAM
periodic ABC meshes, then physically moves only interior vertices using
x_new=x+0.06*sin(pi*x)*sin(pi*y)*sin(pi*z). All six boundary point sets,
periodic partner face geometry, and connectivity remain unchanged.

Quality and geometry gates on both meshes: native v2606 checkMesh
-allGeometry -allTopology must PASS; independent log parser requires
1 degree < maximum nonorthogonality < 65 degrees, max skewness < 3.5,
minimum cell volume > 0, and total volume 1 +/- 1e-6 cubic metres.
The point list must have (n+1)^3 finite vertices, unchanged surface
points and a nontrivial interior x shift (0.025-0.085 m).

Crucially, native OpenFOAM postProcess -func writeCellCentres -time 0
must write actual 0/C. Analytical ABC initial velocity and kinematic
pressure are seeded from these actual physical centroids, not Cartesian
cell-index locations. The independent checker re-evaluates initial ABC
at every centroid to <=5e-12 and requires a nontrivial 0.02-0.09 m
displacement from the undistorted index-based centre.

The existing 3D PISO laminar solver then runs 20 transient steps with
kinematic nu=.01 m2/s, dt=.0025s to t=.05s, 3 pressure correctors and
a conserved passive scalar initially equal to one. The checker reads
real 0.05/U and 0.05/p fields against exact decaying ABC at the *warped*
centres, with pressure gauge normalization, final kinetic energy and
20-step face divergence, net flux and scalar activity-independent budgets.

Predeclared coarse/fine max errors: U <0.32 / <0.22 m/s; gauge-free
p <0.80 / <0.55 m2/s2. Refinement error ratios <0.88 U, <0.92 p.
Corrected max abs face divergence <=1e-7 s^-1, net boundary flux
<=1e-7 m3/s, scalar drift <=1e-8; final kinetic energy 0.95-1.55.
No silent threshold relaxation. Retain all prior official v2606 numerical
gates on exact HEAD, including fixed-pressure open inlet/outlet.

**NOT accepted here:** explicit inner nonorthogonal pressure correction
loop, systematic warp-amplitude sensitivities, MPI on warped meshes,
nonorthogonal buildings/terrain, meteorological inflow, buoyancy, URANS,
HTO/HT transport, soil and air-soil interfaces. Do not close Issue #15.

Source of v2606 postProcess/writeCellCentres:
https://doc.openfoam.com/2606/tools/post-processing/utilities/postProcess/
https://doc.openfoam.com/2306/tools/post-processing/function-objects/field/writeCellCentres/
