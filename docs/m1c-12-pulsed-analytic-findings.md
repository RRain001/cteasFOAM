# M1C-12 findings: native 3-D finite pulse with exact continuum response

Frozen before implementation in commit 047873db0939a2f7f4982cb9e563f30237e6d95b and docs/m1c-12-pulsed-analytic-acceptance.md. These results are development-stage observations until clean exact-SHA regression/full is separately verified.

Ten real native OpenCFD v2606 warped periodic 3-D FVM cases: source-free and uniform-pulse *matched twins* at 8^3 cells with dt=.01,.005,.0025s (20/40/80 steps), plus matched twins at 12^3 and 16^3 with dt=.005s. All have prescribed nonzero 3-D U=(.2,.1,.07)m/s, the M1C-11 full 3-D trigonometric HT/HTO initial profile and kOx=.4s^-1, artificial lambdaT=.2s^-1, DHT=.01 DHTO=.02m²/s. Only the pulse twin configures uniform native source HT 60000 Bq/s and HTO 12000 Bq/s in [0.073,0.137]s. Exact delivered activity is 3840 Bq HT and 768 Bq HTO on *all* grids and clocks. No C++ solver was changed, and a spatially uniform source has no physical Gaussian quadrature bias.

Independent continuum forced response: source-free M1C-11 3-D traveling wave plus a spatially constant source increment. For J_r(t)=(exp(-r*(t-min(t,b)))-exp(-r*(t-a)))/r at t>a, HT injection increment QHT*J_.6; HTO increment QHTO*J_.2+QHT*(J_.2-J_.6). Each native pulse minus matched no-source final field equals an independent Python implementation of the exact finite-step implicit-Euler HT-then-HTO source recurrence, with maximum discrepancy ~5.16e-6 Bq/m³ over all physical native cells.

| dt, s | native HT increment, Bq/m³ | native HTO increment, Bq/m³ | absolute HT continuum time error | absolute HTO continuum time error |
|---:|---:|---:|---:|---:|
| .0100 | 3617.25276225 | 899.624776697 | 10.21132592 | 5.78349836 |
| .0050 | 3622.34610434 | 896.743534573 | 5.11798384 | 2.90225624 |
| .0025 | 3624.90132178 | 895.295632266 | 2.56276639 | 1.45435393 |

HT ratios of adjacent temporal errors .501207 and .500737; HTO ratios .501817 and .501111. These verify first-order backward Euler forcing integration and exact partial-step release, with the source-free three-dimensional spatial discretization removed by genuine matched native field subtraction. Full 3-D forced 16^3 final continuum RMS remains HT 4103.931884 and HTO 1173.111593 Bq/m³, reflecting existing wave transport error. Worst signed all-step HT/HTO/total activity budget ~5.49e-8 Bq.

Native mesh/polyMesh warped vertices, real six cyclic patches, all cell centres and positive volumes, true emitted overlap/Bq logs and every 20/40/80-step activity budget are audited independently. Three negative controls: genuine 16^3 final HT field +1e5 Bq/m³, one logged emitted HT +1000 Bq, and changed sourceEnd from .137 to .147s; all fail for designated physical reasons with no relaxed numerical thresholds.

Scope: limited mathematical verification of a spatially uniform finite-duration HT/HTO pulse on a synthetic periodic 1m^3 domain. Does NOT validate 3-D spatially localized Gaussian forcing, real stacks, time-varying wind, URANS, soil PDE, air-soil coupling, external observed plume, dose, or numerical grid independence. Full independent scientific approval remains NOT_RUN, Draft PR blocked.
