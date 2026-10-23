# M1B-1c: pressure-gradient-equivalent no-slip laminar Poiseuille

The previous Taylor–Green M1B-1a/1b checks tested fully periodic
pressure–velocity correction. M1B-1c verifies actual no-slip wall boundaries
and an independently specified driving mechanism. Both walls are motionless.

The constant streamwise driving acceleration G = 0.1 m/s² is physically
equivalent to a constant streamwise kinematic pressure gradient -dp/dx = G
in an x-periodic channel. It is supplied as an optional OpenFOAM
dimensionedVector 'drivingAcceleration' in transportProperties.
It is *not* a solved imposed pressure drop between physical inlet/outlet
faces. Zero acceleration remains the default for every previously accepted
M1A/M1B benchmark and no earlier transport mode changes.

Domain: x∈[0,1] periodic, y∈[0,1] walls fixed to zero velocity, z∈[0,1]
empty, mesh 4×40×1, nu=0.05 m²/s. The scalar is initially uniform 1 with
zero diffusivity. Start from rest; implicit Euler dt=0.01 s; t_end=1 s.
Kinematic pressure has cyclic x and Neumann wall boundaries.

Exact transient solution for x velocity is:

u(y,t) = G y(1-y)/(2 nu)
        - sum_{n odd} [4G/(nu π³ n³)] sin(n π y) exp(-nu n² π² t).

Independent Python parses the actual 160-cell U field written at t=1.
Acceptance: L∞ analytic velocity error ≤0.003 m/s; transverse velocity
≤1e-8 m/s; x-wise flow homogeneity ≤1e-7 m/s; positive bounded velocity;
100 finite pressure-correction and scalar inventory records; max corrected
face divergence ≤1e-7 s^-1; wall+periodic net boundary flux ≤1e-7 m³/s;
and conserved uniform passive scalar to within 1e-8.

The above tests a no-slip boundary on an embedded 3-D FV mesh with
equivalent pressure-gradient forcing. It does NOT certify real 3-D
building wakes, pressure inlet/outlet BCs, non-orthogonal meshes,
temperature/buoyancy, PIMPLE convergence or turbulence closure.
