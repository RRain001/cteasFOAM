# M1B-1f-4: controlled non-orthogonal mesh-amplitude sensitivity

This is a narrowly scoped **serial** geometry/physics sensitivity study
following PR #16 (actual warped OpenFOAM 3D cell centroids), #17 (two
extra nonorth inner pressure corrections) and #18 (true warped MPI x/z).
No source code change is made to the accepted laminar PISO solver.

For each amplitude a=0.04 and a=0.08 m, generate the same physical
fully periodic unit-cube 8x8x8 and 16x16x16 ABC Beltrami cases. After
real blockMesh, displace interior polyMesh points with

  x_new = x + a sin(pi*x) sin(pi*y) sin(pi*z);

all six periodic boundary sets must remain exactly fixed. Native
checkMesh -allGeometry -allTopology, native writeCellCentres, and
independent ABC U/p initialization using the *actual* warped cell C are
mandatory for all four real cases. All earlier accepted regression
tests (including a=0.06, no-extra and two-inner pressure correction
and warped 2-rank x/z MPI) remain required on the same v2606 CI HEAD.

Run real OpenFOAM v2606 transient 20-step ABC Navier–Stokes with
three outer PISO and two additional inner nonorth pressure corrections
(180 actual pressure matrix solves *per case*). Reuse the original
independent check_warped_abc.py for **each amplitude**: maximum
non-orthogonality in (1,65) degrees; skewness <3.5; positive cell
volumes and total 1m^3; correct initial ABC at C to 5e-12;
face continuity <=1e-7/s; scalar drift <=1e-8; and n8/n16 exact
U and gauge-free p errors plus actual mesh refinement ratios within
the already accepted numerical bounds.

An additional independent sensitivity checker requires at both grid
resolutions that physical max non-orthogonality grows by at least
2 degrees and a factor >1.5 when increasing a=0.04 to a=0.08, that
max skewness also grows, and that the larger warp produces no more
than 1.5 times the earlier U analytical error and no more than 1.6
times the earlier gauge-free p error. All thresholds were set before
the first run; failures must be investigated without silent weakening.

These are **mild-to-moderate controlled distorted meshes** (estimated
~7-14 degrees maximum nonorthogonality). They are not a test of
severe skewness, arbitrary curved/high-aspect-ratio body-conforming
building or terrain meshes, real meteorological inflow, URANS,
buoyancy, HT/HTO transport or soil. Parent Issues #15/#7 stay OPEN
pending those independent physical gates.
