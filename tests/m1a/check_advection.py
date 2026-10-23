#!/usr/bin/env python3
"""Independently verify periodic advective translation and grid improvement."""
from pathlib import Path
import math
import re
import sys

def check(case: Path, n: int):
    path = case / "1" / "c"
    if not path.is_file():
        raise SystemExit(f"FAIL: missing final actual solver field {path}")
    match = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        path.read_text(), re.S,
    )
    if not match:
        raise SystemExit("FAIL: OpenFOAM output has no scalar field list")
    values = [float(v) for v in match.group(2).split()]
    if int(match.group(1)) != n or len(values) != n:
        raise SystemExit("FAIL: wrong number of OpenFOAM finite volumes")
    if not all(math.isfinite(x) for x in values):
        raise SystemExit("FAIL: nonfinite simulated concentration")
    exact = [
        1 + 0.2 * math.cos(2 * math.pi * ((i + 0.5) / n - 0.25))
        for i in range(n)
    ]
    error = max(abs(a-b) for a,b in zip(values,exact))
    drift = abs(sum(values) / n - 1)
    amplitude = (max(values) - min(values))/2
    print(f"M1A advective n={n} max_abs_error={error:.12g} "
          f"mean_drift={drift:.12g} amplitude={amplitude:.12g}")
    if drift > 2e-8 or not 0.1 < amplitude < 0.201:
        raise SystemExit("FAIL: scalar conservation or advective amplitude")
    return error

if len(sys.argv) != 3:
    raise SystemExit("usage: check_advection.py CASE40 CASE80")
e40=check(Path(sys.argv[1]),40)
e80=check(Path(sys.argv[2]),80)
if e40>0.06 or e80>0.035 or not e80<0.8*e40:
    raise SystemExit("FAIL: advection analytic error or grid convergence")
print("PASS: nonzero cyclic advection, conservation, grid-refinement trend")
