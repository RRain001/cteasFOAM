# M1C-11 source-free native three-dimensional continuum verification

This report concerns a **synthetic 1m periodic cube only**. Acceptance was committed first in docs/m1c-11-source-free-3d-analytic-acceptance.md (commit 090c8ac73dc87c2b9f850467ca7ab05dbeac0f38). Unlike M1C-10, the real v2606 C++ simulations have **no source term**, so their analytic field errors cannot originate in Gaussian source quadrature. The three-dimensional analytic eigenfunction varies in all three directions: F=cos(2πx)cos(2πy)cos(2πz). Coupled HT→HTO oxidation and synthetic decay remain active and are compared against the independent closed-form continuum solution at genuine native 3-D OpenFOAM cell centres, with cell-volume-weighted RMS. No modifications to the C++ solver.

Fixed common inputs: unit cube with six native cyclic patches, positive nonorthogonal internally x-warped mesh (0.04*sin(pi*x)*sin(pi*y)*sin(pi*z)), periodic native U, HT, HTO, kOx=.4 s^-1, synthetic lambda=.2 s^-1, DHT=.01 and DHTO=.02 m²/s, initial 1e6+2e5 F and 2e5+5e4 F Bq/m³, dt=.005s, T=.2s. ADVECT fixed vector U=(.20,.10,.07) m/s, DIFFUSE U=(0,0,0). Each mesh: genuine blockMesh/checkMesh/centres/volumes and 40 actually solved implicit Euler HT/HTO/c steps. Analytic modal translation/decay/cross-source obtained independently in Python, not from numerical solver.

| Native regime | Resolution | HT RMS Bq/m³ | HT max Bq/m³ | HTO RMS Bq/m³ | HTO max Bq/m³ |
|---|---:|---:|---:|---:|---:|
| 3D advect | 8³ | 7445.485466 | 20057.24067 | 2044.894240 | 5599.677315 |
| 3D advect | 12³ | 5305.230417 | 14878.52605 | 1493.677725 | 4251.636992 |
| 3D advect | 16³ | 4104.123265 | 11700.85680 | 1173.447625 | 3373.9922 |
| Zero-velocity diffusion | 8³ | 689.998003 | 2151.19687 | 352.081723 | 940.125914 |
| Zero-velocity diffusion | 12³ | 380.576724 | 1474.89592 | 201.078623 | 530.748743 |
| Zero-velocity diffusion | 16³ | 279.378712 | 1117.66452 | 160.744875 | 371.999075 |

Frozen adjacent volume-weighted RMS ratios (8→12, 12→16, both strictly <0.94):
- ADVECT HT 0.712543,0.773599; HTO 0.730443,0.785610.
- DIFFUSE HT 0.551562,0.734093; HTO 0.571113,0.799413.

Largest independent signed native HT/HTO/total 40-step activity budget residual ~5.49e-8 Bq; each genuine physical cyclic net activity flux below frozen 1e-5 Bq/s. All six 3-D native geometry/solver/source-free checks passed. Nonzero 3-D advection changes real final field versus zero U by ~35.5k–36.3k Bq/m³ across tested meshes. Three negative controls reject corrupted real native 4096-cell HT field (analytic/field-Bq mismatch), edited signed net flux, and changed genuine prescribed U input.

**Scientific interpretation:** these results independently verify a genuinely 3-D *source-free modal solution* under prescribed constant U and simple HT→HTO/decay. They isolate **finite-volume transport/reaction numerical error** from Gaussian injection quadrature. They do **not** establish overall mesh independence, a quantitative unique cause for M1C-10 R3, a validation of physically realistic source injection, variable 3-D PISO flow, URANS, terrain, buildings, ground exchange, soil PDE, external plume observation or radiation dose. Full external and non-Bot scientific acceptance remains NOT_RUN; keep stacked research PR Draft.

Developer-run logs are stored privately in /home/lhy/dev/cteasFOAM/reports/m1c11/dev.log; following frozen candidate SHA tests/run evidence must be logged separately. See acceptance contract for exact test thresholds and negative-control specifications.
