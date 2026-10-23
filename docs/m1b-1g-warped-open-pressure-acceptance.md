# M1B-1g: genuinely three-dimensional distorted open-pressure duct

**Status: experimental, unaccepted until actual official v2606 Docker regression and same-SHA CI plus independent scientific review.** Parent Issues #15 and #7 remain open.

## Physical problem and geometry

Extend earlier M1B-1e (8×40×1, empty z) and M1B-1f (fully periodic
warped ABC) **without changing the accepted atmosphere solver**. Use
a unit cube with (8,20,6) hexahedral finite volumes, physical
fixed kinematic pressure p_in=0.1 m²/s² and p_out=0 on its x-end
patches, four stationary no-slip walls at y/z=0,1, initially zero U,
constant nu=0.05 m²/s, no momentum body force, D=0, and a uniform
passive scalar c=1. A negative control uses p_in=p_out=0.

Three independent runs:
- **flat:** Cartesian, pressure-driven square duct;
- **warped:** same physical domain and BCs, but x'=x+0.04 sin(pi x)
  sin(pi y) sin(pi z) only at strictly interior *actual polyMesh vertices*;
- **zero:** warped geometry but zero pressure difference, initialized at rest.

All six exterior surfaces and each patch's surface mesh vertices remain
unaltered. Native v2606 blockMesh, checkMesh -allGeometry -allTopology
and postProcess -func writeCellCentres -time 0 are essential, not emulated.
Pressure is seeded from **actual OpenFOAM C centroids**, p=0.1(1-x)
or zero, avoiding unwarped index-coordinate substitution. Three
outer PISO iterations and two additional inner non-orthogonal corrections
give 20 timesteps × 3 × 3 = **180 actual pressure matrix solves per case**
(t_final=0.4s, dt=0.02s).

The independent analytic reference is the transient square-duct
unidirectional Stokes solution with zero start and Dirichlet no-slip
walls, for G=0.1 m/s²:

u(y,z,t) = sum over positive odd m,n of
[16 G / (pi² m n)] sin(m pi y) sin(n pi z)
[1-exp(-nu pi² (m²+n²)t)]/[nu pi²(m²+n²)].

u_y=u_z=0 and p(x)=G(1-x); G=0 for the negative control.
This is a pressure-driven viscous analytic reference, distinct from
the previous 2-D Poiseuille reference. The production solver remains
nonlinear incompressible laminar Navier–Stokes, so compare within
declared coarse-grid bounds rather than falsely demanding exact equality.

## Predeclared, auditable software/physics checks

These tolerances are adopted **before observing this new test's CFD outputs**:
1. Native checkMesh reports Mesh OK, positive min cell volume, unit
   total volume within 1e-6, skewness below 3.5. Flat nonorthogonality
   <=0.1°, warped >1° and <65°, actual warped-centroid shift >0.01m.
2. All three run exactly 20 complete budgets and 180 actual p solves;
   physical x patch fixedValue p is preserved; four walls remain no-slip.
   Native centroids seed initial linear p within 5e-12 m²/s².
3. Maximum cell-wise |div(phi)| <=1e-7/s, signed total boundary
   volume-flux residual <=1e-7 m³/s, inlet + outlet mismatch
   <=1e-7 m³/s and constant c inventory drift <=1e-8 m³.
4. Pressure-driven cases have q_in<-0.003 and q_out>0.003 m³/s,
   transient 3-D square-duct axial U L∞ analytic error <=0.025m/s,
   pressure linear error <=0.03m²/s² and transverse U <=0.02m/s.
   Paired physical flat/warped U L∞ difference <=0.02m/s and final
   outlet flow relative change <=20%.
5. Zero-dp warped control has U, p, q_in and q_out <=1e-8 in
   corresponding units; all measured fields and logs must be finite.
6. Existing M1A/M1B-1f regression suite must pass on the **same HEAD**
   separately; a targeted new test by itself is not full acceptance.

The checker fails closed on missing final U/p fields, missing native
centroids or checkMesh metrics, truncated pressure solve logs, missing
flux/continuity/scalar records, invalid physical pressure patches,
and unexpected wall boundaries. No synthetic external CFD or
observational data are used.

## Interpretation limits

This tests interior skew and nonperiodic *pressure* boundaries on a
full 3-D wall-bounded domain, not a realistic atmospheric inlet, not
turbulent or buoyant building wakes, not variable air density, not
HT/HTO, not deposition, and not soil/air–soil coupling. Its pressure
BCs at both ends and coarse low-Re mesh remain a narrow canonical
case. Fluid velocity and passive scalar carry no Bq units here.
Native software acceptance **does not constitute independent scientific
validation**. The source solver still has no validated URANS closure.

Run in official OpenFOAM v2606 development image with
bash tests/m1b1/run_warped_open_pressure.sh from /workspace mount.
This experimental runner is not an existing GitHub required check;
do not report GitHub Actions success as covering this new test unless
the workflow is explicitly extended and verified.

## M1B-1g-1: two-rank x/z MPI on the *same* warped open-pressure duct

This additional candidate extends M1B-1g without changing the solver, the
v2606 package, or any earlier acceptance tolerance. It specifically verifies
a **physical fixed-pressure Dirichlet gauge** together with a non-orthogonal
interior mesh and two orientations of processor boundaries. The original
M1B-1g serial flat/warped/zero suite remains a separate required test.

Run `bash tests/m1b1/run_warped_open_pressure_mpi.sh` in the official
v2606 Python-enabled Docker image. The executable is compiled from the
candidate source. Native blockMesh, checkMesh, writeCellCentres,
decomposePar, two-rank mpirun and reconstructPar are executed, not mocked.

Cases and immutable physical inputs:

- Serial warped pressure-driven duct, x-split two-rank, z-split two-rank;
  960 cells, identical polyMesh/points, native C and initial U/p/c.
- Two-rank z-split warped **zero-pressure-drop** negative control. Its
  geometry, initial U and c match the driven cases; its p is *intentionally*
  zero and must NOT match their linearly seeded initial pressure.
- Each requires 20 completed timesteps, 3 PISO correctors, 2 additional
  non-orthogonal correctors, i.e. 180 actual pressure matrix solves.
  Pressure patches remain xMin=0.1 and xMax=0 m²/s² for driven cases.

Predeclared MPI equivalence gates: all existing M1B-1g independent
transient square-duct Stokes, checkMesh, signs of q_in/q_out, no-slip,
zero-drive, continuity <=1e-7 /s, total boundary flux <=1e-7 m³/s,
and scalar inventory drift <=1e-8 m³ checks apply to **each**
reconstructed case. In addition, both MPI decompositions relative to
the serial reference must satisfy cell-wise U L∞ <=0.005 m/s,
**absolute** pressure L∞ <=0.01 m²/s² (the pressure is fixed-value,
so gauge subtraction would hide a real error), scalar c L∞ <=1e-7,
max inlet/outlet flow history discrepancy <=0.001 m³/s, relative
final signed flow difference <=5%, kinetic energy difference <=0.003
m5/s2. x/z orientations also require U L∞ <=0.005 m/s and
absolute p L∞ <=0.01 m²/s². These limits are bounded engineering
regressions, not proven asymptotic convergence. OpenFOAM reconstruction
must supply real final U/p/c on the original physical mesh.

Official GitHub CI currently exercises the earlier M1B-1f gate but **not**
this new runner; a successful existing PR check or local Docker scoped full
must not be relabeled as M1B-1g-1 scientific approval. No external
experimental observation, turbulence, real atmospheric boundaries,
HT/HTO or air–soil coupling is validated by the canonical duct.
