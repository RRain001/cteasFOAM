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

## M1B-1g-2: time-step self-convergence on the real warped duct

This is a **separate supplemental candidate**, not a replacement for the
accepted M1B-1g-1 x/z MPI check or an independent scientific review.

The official v2606 test uses one physically identical OpenFOAM 960-cell
3-D warped duct with the same actual polyMesh vertices, native C, initial
U/p/c, six boundary patches, fluid viscosity nu=0.05 m²/s, fixed
kinematic pressures 0.1 and 0 m²/s², four no-slip walls, t_final=0.4 s,
and the same PISO/nonorthogonal controls. Only deltaT changes:
0.02, 0.01 and 0.005 s, giving 20/40/80 physical time steps and
180/360/720 actual pressure solves. Native blockMesh, checkMesh,
postProcess writeCellCentres, wmake and solver are all real; no
precomputed output fields may substitute.

**A priori acceptance rules (before the first computation)**:

- Bitwise identical physical mesh and initial U, p, c, C fields; actual
  warped vertices = 665; official checkMesh Mesh OK, warped
  nonorthogonality 1–65°, positive minimum volume, skew <3.5.
- Every dt produces its complete expected number of pressure solves
  (3 outer PISO × 3 nonorthogonal per physical step), signed boundary
  fluxes, continuity/energy budgets and scalar inventories. No missing
  or nonfinite measurements.
- For *each* dt: absolute cellwise max |div(phi)| ≤1e-7/s; total
  boundary flux and inlet+outlet flux imbalance ≤1e-7 m³/s; scalar
  c/volume drift ≤1e-8 m³. Nonzero throughflow q_in <-0.003,
  q_out >0.003 m³/s, and correct fixed-pressure/no-slip BCs.
- The independently evaluated transient 3-D square-duct Stokes
  reference still has max axial velocity discrepancy ≤0.025 m/s;
  linear pressure discrepancy ≤0.03 m²/s²; transverse U ≤0.02 m/s.
- Compare full 3-D velocity at the same physical centroids and the
  final **signed outlet volume flux** at each dt. Define Δ20/10 =
  ||numerical(dt=.02) - numerical(dt=.01)|| and Δ10/05 similarly;
  the ratio Δ10/05/Δ20/10 must be ≤0.85 **for both** velocity L∞ and
  outlet volume flux, with Δ20/10 >1e-11 in the respective units to
  avoid a false zero-numerator convergence claim. Measure absolute p
  L∞ differences too (physical fixed pressure fixes the gauge), but
  do not impose a monotonic p temporal trend in this first-stage test.

The target ratio 0.85 is deliberately less strict than asymptotic
first-order Euler's ~0.5, accounting for coupled pressure iteration
effects; it is not to be weakened after observing outputs. Tests
must disclose any failure, not silently retune. This experiment
checks temporal **self-convergence at fixed spatial mesh** only;
it cannot demonstrate spatial grid convergence, continuum solution
uniqueness, turbulence accuracy, real meteorological inflows or
experimental scientific validity. URANS, HT/HTO and air–soil coupling
remain NOT_IMPLEMENTED / NOT_VALIDATED.

Runner: `bash tests/m1b1/run_warped_open_pressure_timestep.sh`.

### Local precommit observation (not independent science acceptance)

Official OpenFOAM v2606 Docker targeted run returned exit 0 on the
three dt cases after the above thresholds were written, with identical
physical mesh and initial fields. The observed L∞ successive velocity
differences were 1.79069615356e-4 and 9.32797613908e-5 m/s,
giving ratio **0.520913395638**. Signed outlet-flow differences
were 1.16976573357e-4 and 5.8875156009e-5 m³/s,
ratio **0.503307237675**. The finest-dt max |div(phi)| was
7.88750273642e-14 /s. All three cases retained 20/40/80 complete
budgets, 180/360/720 real p matrix solves and positive throughflow.

The independent Stokes **absolute** U error increased slightly from
0.00342234987151 at dt=.02 to 0.00362780220956 m/s at dt=.005,
remaining below 0.025 m/s. This is not a claim of monotonically
improved total error: changing temporal truncation at a fixed coarse
mesh can expose spatial, boundary-condition or pressure-iteration
biases. The physical pressure differences between refinements also
did not follow an imposed convergence ratio. A mesh-refinement
study and independent external review remain distinct requirements.

## M1B-1g-3: three-level *spatial* mesh refinement (predeclared gate)

This increment keeps the existing laminar OpenFOAM v2606 solver,
physical unit cube, fixed kinematic pressure 0.1 → 0 m²/s², 4 no-slip
walls, ν=0.05 m²/s, D=0 and initial U=0/c=1 unchanged.
It uses identical **physical** x-only interior deformation
x'=x+0.04 sin(πx)sin(πy)sin(πz) at three actual native hexahedral
finite-volume resolutions: **8×20×6 (960), 12×30×9 (3240),
16×40×12 (7680)**. Boundary patches are physically unchanged.
The new optional positional grid argument is limited to the existing
case generator's make/warp/seed actions; omitting it retains the original
8×20×6 regression behavior.

The temporal method is held constant: Δt=0.01 s, final time 0.4 s,
40 full time steps and 360 real p matrix solves per grid. Each
mesh is independently regenerated via v2606 `blockMesh`, genuinely
warped in its `polyMesh/points`, checked with `checkMesh`,
and seeded from OpenFOAM's actual `writeCellCentres`, rather than
reusing the coarse mesh or index-based centroids.

**Acceptance rules fixed before running the new cases:**
- Each mesh is physically matched point-by-point (order-independent)
  to the known deformed coordinate lattice, including unchanged
  physical exterior vertices; no fake mesh or cached final field.
- `checkMesh` returns Mesh OK, >1° and <65° maximum nonorthogonality,
  <3.5 skew, positive minimum cell volume and total unit volume
  within 1e-6. Every native field has exactly N cells.
- Every grid actually completes 40 PISO steps, 360 p solves,
  signed inlet/outlet fluxes, cellwise divergence and scalar budgets;
  preserved inlet/outlet pressure and all four no-slip walls,
  finite numerical output and exactly prescribed initial U/p/c.
- Each grid must have max |div(phi)|, net boundary and inlet/outlet
  imbalance ≤1e-7 in their respective units, scalar inventory/c drift
  ≤1e-8, negative inlet flow <−0.003 m³/s, positive outlet
  flow >0.003 m³/s; independent Stokes U L∞ ≤0.025 m/s,
  linear pressure L∞ ≤0.03 m²/s² and transverse U L∞ ≤0.02 m/s.
- Compute **RMS axial-velocity error** against the independently
  implemented transient square-duct Stokes series at each native
  cell centre, `sqrt(sum_i (Ux_i-u_exact(y_i,z_i,t))²/N)`.
  This is an unweighted cell-centre RMS on each logically uniform
  hexahedral mesh, not a volume-weighted/interpolated cross-mesh
  difference. Each finer resolution must reduce this RMS by at
  least **10%** (ratio <0.9), and finest/coarsest ≤**0.7**.
  Also report axial U L∞, signed outlet flow, pressure L∞,
  mesh-quality and conservation metrics; no unsupported theoretical
  second-order claim follows from this limited three-grid study.

This measures agreement with an analytic unidirectional Stokes
reference for a **finite-length nonlinear incompressible
Navier–Stokes** calculation with fixed static pressure boundaries.
There can be legitimate pressure-end and convection differences,
and changing cell-centre sample points alters the error norm.
Thus an internal grid-error decrease is evidence of a controlled
numerical refinement trend, **not a formal GCI estimate, measured
order of accuracy, a continuum error bound or external
scientific field validation**. A separate dt convergence result is
documented above. A genuinely independent science review is
still mandatory for PR merge.

Runner: `bash tests/m1b1/run_warped_open_pressure_spatial.sh`
from a read-only repository mount inside the official v2606
Ubuntu Docker image. The checker is new and cannot waive
existing serial/MPI and temporal regression gates.

### M1B-1g-3 observed initial spatial result (scope: candidate only)

Official v2606 Docker on Ubuntu, `Δt=0.01`, at `t=0.4`;
**960/3240/7680** real cells, 40 time steps and 360 pressure
matrix solves *per case*, yielded exact native-centroid RMS axial U
errors **1.59829213498e-3, 8.73701386129e-4,
5.10954721773e-4 m/s**. Ratios were **0.546646865745** and
**0.584816196798**; finest/coarsest **0.319687941017**.
All passed the separately written ≤0.9 adjacent and ≤0.7 overall
criteria. Axial U L∞ errors were **0.00355757579,
0.00218359681 and 0.00134930797 m/s**, also decreasing.

Real OpenFOAM mesh nonorthogonality: **6.64834°, 7.03556°,
7.03097°**. Finest grid max |div(phi)| **2.69e-13 /s**,
scalar inventory drift zero to printed precision, 0.4 s
final signed outlet flows **0.02655447, 0.02572988,
0.02538801 m³/s**, respectively. Since outlet flow varies with
resolution, this is not an assertion of a known continuum discharge
or a formally measured convergence order.

During initial *precommit* testing, the new independent point-set
auditor rejected real mesh geometry because raw floating-point
lexicographic sorting can reorder x-tied symmetric (y,z) lattice
vertices. The geometry check was corrected to compare unordered
multisets of coordinates rounded to 9 decimal places (≤5e-10 m
quantization per coordinate), without modifying solver, physical
mesh, scientific convergence or conservation gates. Subsequent
three-resolution Docker run returned **PASS, exit 0**. This
failure and checker correction are retained to avoid implying
first-attempt numerical acceptance.

## M1B-1g-4: separate upstream OpenFOAM icoFoam solver comparison

This increment provides a second **numerical implementation**:
the genuine unmodified `icoFoam` binary from the same pinned
official OpenCFD/Keysight v2606 distribution. It is **not** an
independent human review, external measurements, validation of URANS
or a claim of scientific signoff. The CTEAS solver remains experimental.

Compare a controlled 8×20×6 interior-warped physical fixed-pressure
unit duct and a zero-pressure-drop negative control. For **each**
pair the official icoFoam and CTEAS runs use byte-identical native
`polyMesh`, OpenFOAM-generated C centroids, U/p/c initial values,
viscosity nu=0.05 m²/s, pressure Δp=0.1 or 0 m²/s², four
no-slip walls, `fvSchemes`, `fvSolution`, timestep 0.02 s,
t_final=0.4 s and PISO counts (three pressure corrections, two
additional nonorthogonal corrections). Only the
`application` entry of the `controlDict` differs, and the CTEAS
transport dictionary includes its passive scalar settings that
stock icoFoam simply ignores. The solver runs are actual, not
replayed/fabricated fields.

**Predeclared criteria before running the cross-solver case:**
- Exact binary input equality (except `application`), with actual
  upstream `icoFoam` completing 20 steps and 180 genuine
  p solves per case; OpenFOAM mesh quality and native centres
  inspected independently.
- CTEAS retains *all* existing distorted-duct analytic, pressure
  BC, physical signed flux and scalar-conservation tests.
- Because these are distinct legitimate PISO implementations,
  not bitwise-identical algorithms, accept final field
  **componentwise U L∞ ≤0.007 m/s, U RMS ≤0.003 m/s, absolute
  fixed-gauge p L∞ ≤0.007 m²/s²**. Compare raw absolute p because
  the boundary fixedValues anchor its physical gauge.
- Upstream `icoFoam` must itself remain within independent
  transient Stokes axial-U L∞ ≤0.025 m/s and transverse-U
  L∞ ≤0.02 m/s.
- Both independent programs' zero-drive controls must end with
  maximum |U| and |p| ≤1e-8 in their relevant units;
  fields/logs must be finite and complete.

Unlike the prior analytic Stokes series, upstream `icoFoam`
solves the nonlinear incompressible pressure-corrected momentum PDE
directly. Agreement would strengthen implementation-level
verification but not give a continuum or real-site error bound.
If the cross-implementation gate fails, retain a visible
failure and investigate rather than retrofitting its thresholds.

Runner: `bash tests/m1b1/run_icofoam_cross_solver.sh`.

### M1B-1g-4 upstream solver observation (experimental scope)

On the official OpenCFD v2606 image, with identical physical mesh,
initial U/p/c, timestep, BCs and `fvSchemes/fvSolution`,
the real **unmodified upstream icoFoam binary** and CTEAS
`pisoLaminar` produced **exactly identical printed final U
and absolute kinematic p**: componentwise U L∞ = 0, RMS = 0
and p L∞ = 0 at the 960-cell warped driven case, as well
as zero U/p under Δp=0. Both completed **20 steps with 180
p matrix solves**. The upstream U analytic Stokes L∞ error
was 0.00342234987151 m/s. Separately injecting +0.02 m/s
into one native upstream final U vector caused the independently
executed checker to fail **as required**.

The first exploratory run failed with the authentic upstream
`icoFoam` error `Entry 'pFinal' not found in dictionary
system/fvSolution/solvers`. Official icoFoam requires that
final-iteration solver entry; both inputs were then given the
identical `pFinal { $p; relTol 0; }` setting. The CTEAS
solver ignores that extra entry. No physical PDE, time, mesh,
analytic/conservation threshold or solver C++ change was made.
After this **configuration correction**, the repeated real
Docker test and corruption-negative check passed, exit 0.

Matching printed numeric fields indicates the current
experimental PISO mode reproduces the same discrete upstream
OpenFOAM solver path under this canonical configuration. It
does **not** establish an independent discretization or
experimental/real-atmosphere validation; shared OpenFOAM
libraries, algorithms and conditions mean the two outputs
are not independent physical observations.

## M1B-1g-5: reversed fixed-pressure gradient sign/flow validation (pre-run gate)

The same official v2606 `run_icofoam_cross_solver.sh` now executes a
third, *reverse* run alongside the already accepted `warped` and `zero`
modes. Before observing numerical results, this increment declares:

- Identical real 960-cell warped polyMesh, true native cell centres,
  four no-slip walls, `nu=0.05`, `D=0`, initial velocity zero and
  uniform dimensionless `c=1`. Only kinematic pressure reverses from
  `p(xMin)=+0.1, p(xMax)=0` to `p(xMin)=-0.1, p(xMax)=0`
  (m²/s²). Initialise `p=-0.1(1-x)` at *native* cell centres.
- Execute both the experimental CTEAS solver and unmodified official
  `icoFoam` with byte-equal input files within each mode; require
  20 genuine transient steps, 180 actual pressure matrix solves, real
  output U/p, correct fixed pressure and zero wall velocity. The
  checker requires reverse `q_in>+0.003` and
  `q_out<-0.003 m³/s` (signs reversed relative to forward),
  full 20-step corrected-flux budgets with
  `max|div(phi)| <=1e-7 /s`, net/boundary-pair flux
  `<=1e-7 m³/s` and uniform-scalar inventory drift
  `<=1e-8 m³`. The same independent signed Stokes
  series applies with pressure-gradient parameter `G=-0.1 m/s²`.
- Reversed CTEAS vs official `icoFoam`: physical full-field
  U L∞ `<=0.007 m/s`, U RMS `<=0.003 m/s`,
  **absolute** fixed-gauge p L∞ `<=0.007 m²/s²`;
  upstream analytic signed-axial-U L∞ `<=0.025 m/s`,
  transverse U L∞ `<=0.02 m/s`.
- The forward/reverse *candidate solver* is required to have
  `max|U_forward+U_reverse| <=0.01 m/s` at identical
  physical centres and max of the two final signed inlet/outlet
  summed flow residuals `<=0.003 m³/s`. This is a
  **bounded low-Re sign-symmetry test**, not exact odd symmetry:
  nonlinear inertial terms and end effects prevent an exact Stokes
  superposition assertion. These bounds are fixed before this
  increment's first Docker run, not tuned to its observed outcomes.
- Existing positive-drive, no-drive and deliberate +0.02 m/s corrupted
  official-output negative test must remain passing; future external
  scientific or independent human review remains NOT_RUN.

This is a stronger negative-pressure boundary and flux-sign regression on a
mild warped canonical duct. It does not validate atmospheric inflow/outflow
with pressure-turbulence coupling, building walls, URANS, HT/HTO, soil or
any air–soil activity exchange. The scalar boundary condition remains
uniform-one so reverse advection is not a tracer inflow pulse test.
