#!/usr/bin/env python3
"""Independent analytic diffusion reference, not an OpenFOAM implementation.

For zero-velocity diffusion in [0,1] with Neumann walls:
c(x,0)=1+0.2*cos(pi*x), c(x,t)=1+0.2*exp(-D*pi²*t)*cos(pi*x).
The mean must remain 1 (closed domain, no sources).
"""
import math
from pathlib import Path
import re
import sys

case = Path(sys.argv[1])
solver_log = Path(sys.argv[2]).read_text()
path = case / "2" / "c"
if not path.is_file():
    raise SystemExit(f"FAIL: solver did not write final field {path}")
content = path.read_text()
m = re.search(r"internalField\\s+nonuniform\\s+List<scalar>\\s+(\\d+)\\s*\\((.*?)\\)\\s*;", content, re.S)
if not m:
    raise SystemExit("FAIL: final c field is not a nonuniform scalar list")
n = int(m.group(1))
values = [float(v) for v in m.group(2).split()]
if n != 40 or len(values) != n or not all(math.isfinite(v) for v in values):
    raise SystemExit(f"FAIL: invalid final scalar field: n={n}, entries={len(values)}")
t, diffusivity = 2.0, 0.05
expected = [
    1.0 + 0.2 * math.exp(-diffusivity * math.pi**2 * t)
    * math.cos(math.pi * (i + 0.5) / n)
    for i in range(n)
]
max_error = max(abs(a - b) for a, b in zip(values, expected))
mean_error = abs(sum(values) / n - 1.0)
initial_amplitude = 0.2
final_amplitude = (max(values) - min(values)) / 2
print(
    f"M1A analytic Neumann diffusion: max_abs_error={max_error:.9g}; "
    f"closed_domain_mean_error={mean_error:.9g}; "
    f"final_amplitude={final_amplitude:.9g}"
)
if "CTEAS_M1A initial_integral=" not in solver_log or "CTEAS_M1A integral=" not in solver_log:
    raise SystemExit("FAIL: missing conservation audit trace in solver log")
if max_error > 1.5e-3:
    raise SystemExit("FAIL: analytic diffusion error exceeds 1.5e-3")
if mean_error > 2e-8:
    raise SystemExit("FAIL: closed-domain scalar integral drift exceeds 2e-8")
if not (0.05 < final_amplitude < 0.12 and final_amplitude < initial_amplitude):
    raise SystemExit("FAIL: expected nontrivial diffusion decay was not observed")
print("PASS: analytic profile, nontrivial decay and closed-domain integral")
