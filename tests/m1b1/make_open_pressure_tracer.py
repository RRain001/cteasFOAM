#!/usr/bin/env python3
"""M1B-1g-6: configure nonuniform zero-start scalar on established warped duct.

Only c boundary/initial conditions and audit flags are changed.
All flow physics and actual OpenFOAM mesh, U, p remain unchanged.
"""
from pathlib import Path
import re
import sys

if len(sys.argv) != 3 or sys.argv[2] not in ("warped", "reverse", "zero"):
    raise SystemExit("usage: make_open_pressure_tracer.py CASE warped|reverse|zero")
case, mode = Path(sys.argv[1]), sys.argv[2]
forward = mode != "reverse"
in_min, in_max = (1, 0) if forward else (0, 1)
c_path = case / "0/c"
data = c_path.read_text()
start, end = data.index("internalField"), data.index("boundaryField")
if end <= start:
    raise SystemExit("FAIL: malformed initial c field")
data = data[:start] + "internalField uniform 0;\n" + data[end:]
pat = r"\bxMin\s*\{[^{}]*\}\s*\bxMax\s*\{[^{}]*\}"
replacement = (
    f"xMin {{ type inletOutlet; inletValue uniform {in_min}; value uniform 0; }}\n"
    f" xMax {{ type inletOutlet; inletValue uniform {in_max}; value uniform 0; }}"
)
data, count = re.subn(pat, replacement, data, count=1)
if count != 1 or data.count("inletOutlet") != 2:
    raise SystemExit("FAIL: expected both native scalar pressure-end patches")
c_path.write_text(data)
transport = case / "constant/transportProperties"
with transport.open("a") as f:
    f.write("auditPassiveScalarFlux true;\n")
schemes = case / "system/fvSchemes"
scheme = schemes.read_text()
old = "fluxRequired { default no; p; }"
if scheme.count(old) != 1:
    raise SystemExit("FAIL: expected native scalar fluxRequirement declaration")
schemes.write_text(scheme.replace(old, "fluxRequired { default no; p; c; }"))
print(f"M1B-1g-6 configured mode={mode} c0=0 inletValue_xMin={in_min} "
      f"inletValue_xMax={in_max} D=0 conservative cEqn flux")
