# M1A-2: nonzero prescribed periodic advection acceptance

Two transparent generated OpenFOAM cases with 40 and 80 cells in x,
single-cell empty y/z, cyclic x, uniform U=(0.25,0,0) m/s and D=0.
Dimensionless initial scalar: c(x,0)=1+0.2 cos(2πx).
Analytic translation after 1s: c(x,1)=1+0.2 cos(2π(x-0.25)).
Time step is 0.005s and final time 1s; FVM uses implicit Euler + upwind.

Actual OpenFOAM output is compared with the exact solution by independent
Python. Acceptance: coarse max error <=0.06; fine <=0.035, fine/coarse <0.8;
volume-mean drift <=2e-8 on both; nonzero final Fourier amplitude.
This detects inactive convection, wrong advection sign, total inventory
drift and lack of grid improvement. The generated input is reproducible.

Not validated: variable 3D wind, URANS, inflow/outflow, MPI, HT/HTO,
deposition, soil transport and air-soil coupling. See Issue #3.
