# M1B-1f-2: actual inner non-orthogonal pressure correction on warped ABC

Introduces a **distinct** optional parameter nNonOrthogonalCorrectors in
transportProperties. Default 0 retains exactly one pressure equation solve
per old outer PISO corrector and preserves all existing physical modes.
Valid range 0..6; invalid values fail. The already existing independent
nPressureCorrectors parameter still controls the outer PISO passes.

OpenFOAM reference PISO pattern: inside each outer correction, recompute
the finite-volume p Laplacian equation with the previous p field for each
non-orthogonal inner iteration, while holding phiHbyA unchanged, and
update corrected phi ONLY on the final non-orthogonal pressure iteration.
All MPI ranks participate in setReference for every pressure solve,
providing only a unique valid master gauge reference cell. U correction
uses the finally corrected pressure.

The test uses the **same physical** previously scientifically accepted
8x8x8 and 16x16x16 warped ABC meshes, x displacement
0.06 sin(pi*x) sin(pi*y) sin(pi*z), native v2606 checkMesh, true
OpenFOAM writeCellCentres geometry, and analytic initial U/p at each
physical centroid. The only intended change relative to PR #16 is
nNonOrthogonalCorrectors=2 (three inner pressure solves) with three
outer PISO correctors over 20 steps, so exactly **180 actual pressure
matrix solves per grid** must appear in the solver log, plus 20
nonorthogonal parameter records.

A separate independent checker reuses all physical/geometry gates from
PR #16: actual nonorth angle, skewness, volume, true initial centroids,
full final ABC velocity and gauge-free pressure, scalar conservation and
face continuity with mesh refinement. It then checks that U/p errors
are no worse than **120%** of the accepted default-0 reference at
the same mesh and dt (prior official run #39895216715):

  n=8: U max error 0.172379028547 m/s; p 0.553776125146 m2/s2
  n=16: U max error 0.108293338904 m/s; p 0.289513430093 m2/s2

The 20% relative non-regression limit is deliberately predefined before
executing this branch. Extra nonorth pressure iterations are not
assumed to universally improve analytic L-infinity errors, which
also depend on space/time/upwind discretization.

All earlier official OpenFOAM v2606 CI gates must still PASS on exact
branch HEAD, including the M1B-1f-1 default-zero warp benchmark,
ABC 3D two-rank MPI and real static pressure boundaries.
Remaining Issue #15: warped-mesh 2-rank MPI, multiple warp amplitudes,
strong skewness and real complex building mesh; Issue #7 stays open.
No URANS/buoyancy/tritium/soil capability claimed.

Relevant source pattern:
https://api.openfoam.com/2406/incompressible_2pisoFoam_2pEqn_8H_source.html
