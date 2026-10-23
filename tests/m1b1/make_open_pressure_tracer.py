#!/usr/bin/env python3
"""M1B-1g-6/7: configure zero-start scalar on established warped duct.

The optional D=0.01 extension changes only scalar diffusivity in addition
to c boundary/initial conditions and audit flags. The physical mesh,
momentum equations, U and p initial fields remain unchanged.
"""
from pathlib import Path
import re
import sys

if len(sys.argv) not in (3, 4) or sys.argv[2] not in ("warped", "reverse", "zero"):
    raise SystemExit("usage: make_open_pressure_tracer.py CASE warped|reverse|zero [D=0|0.01]")
case, mode = Path(sys.argv[1]), sys.argv[2]
diffusivity = sys.argv[3] if len(sys.argv) == 4 else "0"
if diffusivity not in ("0", "0.01"):
    raise SystemExit("FAIL: only predeclared D=0 or D=0.01 is supported")
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
transport_text = transport.read_text()
transport_text, n = re.subn(
    r"(?m)^D \[0 2 -1 0 0 0 0\] 0;$",
    f"D [0 2 -1 0 0 0 0] {diffusivity};", transport_text,
)
if n != 1:
    raise SystemExit("FAIL: expected exactly one physical diffusivity entry")
transport.write_text(transport_text + "auditPassiveScalarFlux true;\n")
schemes = case / "system/fvSchemes"
scheme = schemes.read_text()
old = "fluxRequired { default no; p; }"
if scheme.count(old) != 1:
    raise SystemExit("FAIL: expected native scalar fluxRequirement declaration")
schemes.write_text(scheme.replace(old, "fluxRequired { default no; p; c; }"))
print(f"M1B-1g-6/7 configured mode={mode} c0=0 inletValue_xMin={in_min} "
      f"inletValue_xMax={in_max} D={diffusivity} conservative cEqn flux")
