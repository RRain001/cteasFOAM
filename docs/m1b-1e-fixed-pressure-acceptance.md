# M1B-1e: physical fixed-pressure inlet and outlet boundary verification

## Why this is distinct from M1B-1c

The previously accepted Poiseuille test prescribed a constant body
acceleration within an x-periodic channel. This milestone uses a genuinely
**open, nonperiodic** streamwise direction: xMin and xMax are patch
(not cyclic), and the specified kinematic pressures are 0.1 and 0.0 m²/s²,
respectively. The x-boundary velocity condition is
'pressureInletOutletVelocity' (OpenCFD documentation). Transverse y walls
are stationary no-slip; pressure at walls is 'fixedFluxPressure'.
The domain has 8×40×1 FV cells, unit length/gap, empty z.

A separate negative-control case has identical geometry, solver, grid,
boundary-condition classes and initial velocity U=0, but zero pressure
on **both ends**. Both cases have NO drivingAcceleration, nu=0.05 m²/s,
D=0 for scalar c, optional PISO 3 correctors, dt=0.01s, t=1s
(100 time steps). The input p field starts from an exact linear gradient;
this is not an imposed body force, and the solver must maintain pressure
boundary conditions through the Poisson corrections.

The actual solver prints per-step corrected **signed inlet and outlet
volumetric fluxes** by summing phi on xMin and xMax separately, as well as
closed net flux and maximum cell-wise |div(phi)|. The independent checker
reads the actual final 320-cell OpenFOAM U and p values, checks the
final p patch classes and fixedValue magnitudes, real inward and
outward nonzero flow, the transient Poiseuille Fourier-series profile,
pressure linearity, nonphysical transverse velocity, streamwise uniformity
and scalar uniformity/inventory. The negative control must have near-zero
flow everywhere and no inlet or outlet throughflow.

This deliberately tests **static pressure at both ends**, a configuration
that OpenFOAM identifies as potentially poorly conditioned and that
should not be generalized from one simple low-Re planar channel
(https://doc.openfoam.com/2306/tools/processing/boundary-conditions/common-combinations/).
Preserve all prior M1A/M1B-1d numerical acceptance checks on the exact
same official OpenFOAM v2606 CI job.

No generic building-wake inflow, pressure-recovery turbulence, large
Reynolds number URANS, nonorthogonal geometry, buoyancy, HT/HTO, soil,
or conservative air-soil coupling is claimed.
