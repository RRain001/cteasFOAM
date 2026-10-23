#!/usr/bin/env python3
"""Independent backward-Euler activity and budget checker (no OpenFOAM imports)."""
import math
import re
import sys
from pathlib import Path

if len(sys.argv) != 5:
    raise SystemExit("usage: check_activity_box.py CASE SOLVER_LOG K LAMBDA")
case, logpath = Path(sys.argv[1]), Path(sys.argv[2])
k, lam = float(sys.argv[3]), float(sys.argv[4])
dt, n, ht, hto = 0.05, 20, 1e6, 2e5
tol = 1e-6
log = logpath.read_text()
def fail(message):
    raise SystemExit("FAIL: M1C " + message)

if len(re.findall(r"(?m)^Time = ", log)) != n:
    fail("expected 20 real OpenFOAM time steps")
if not (case / "constant/polyMesh/owner").is_file():
    fail("real OpenFOAM mesh missing")
mesh_log = (case / "mesh-check.log").read_text()
if not re.search(r"(?m)^\s*cells:\s*24\s*$", mesh_log):
    fail("expected actual 24-cell 3D mesh")
bbox = re.search(r"Overall domain bounding box\s*\(([^)]+)\)\s*\(([^)]+)\)", mesh_log)
if not bbox:
    fail("missing mesh bounding box")
lo, hi = ([float(x) for x in v.split()] for v in bbox.groups())
if len(lo) != 3 or len(hi) != 3 or any(abs(x-y) > 1e-12 for x, y in zip(lo+hi, [0,0,0,1,1,1])):
    fail("expected unit-cube mesh")
if not re.search(r"(?m)^Time = 1(?:\.0+)?$", log):
    fail("expected final time 1.0")

def values(field):
    path = case / "1" / field
    if not path.is_file():
        fail(f"missing final {field} field")
    s = path.read_text()
    if not re.search(r"dimensions\s+\[0 -3 -1 0 0 0 0\]", s):
        fail(f"{field} incorrect Bq/m3 dimensions")
    m = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", s)
    if m:
        return [float(m.group(1))] * 24
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+24\s*\((.*?)\)\s*;", s, re.S)
    if not m:
        fail(f"{field} full field unreadable")
    result = [float(v) for v in m.group(1).split()]
    if len(result) != 24:
        fail(f"{field} wrong cell count")
    return result

p = re.compile(r"CTEAS_M1C activity_HT_Bq=([-+0-9.eE]+) activity_HTO_Bq=([-+0-9.eE]+)")
b = re.compile(r"CTEAS_M1C closed_HT_residual_Bq=([-+0-9.eE]+) closed_HTO_residual_Bq=([-+0-9.eE]+) closed_total_residual_Bq=([-+0-9.eE]+)")
inventories = [tuple(map(float, v)) for v in p.findall(log)]
budgets = [tuple(map(float, v)) for v in b.findall(log)]
if len(inventories) != n or len(budgets) != n:
    fail(f"expected 20 actual inventory/budget records, got {len(inventories)}/{len(budgets)}")
old_ht, old_hto = ht, hto
max_field = max_budget = max_total = 0.0
for i, ((a, z), (r_ht, r_hto, r_total)) in enumerate(zip(inventories, budgets), 1):
    ht = ht / (1 + dt*(k + lam))
    hto = (hto + dt*k*ht) / (1 + dt*lam)
    max_field = max(max_field, abs(a-ht), abs(z-hto))
    calc_ht = a-old_ht + dt*(k+lam)*a
    calc_hto = z-old_hto + dt*lam*z-dt*k*a
    calc_total = a+z-old_ht-old_hto + dt*lam*(a+z)
    max_budget = max(max_budget, abs(calc_ht), abs(calc_hto), abs(r_ht), abs(r_hto),
                     abs(calc_ht-r_ht), abs(calc_hto-r_hto))
    max_total = max(max_total, abs(calc_total), abs(r_total), abs(calc_total-r_total))
    if not all(math.isfinite(x) and x >= -1e-12 for x in (a, z)):
        fail(f"negative or nonfinite inventory at step {i}")
    old_ht, old_hto = a, z

for field, expected in (("HT", ht), ("HTO", hto)):
    actual = values(field)
    if not all(math.isfinite(v) and v >= -1e-12 for v in actual):
        fail(f"negative or nonfinite {field} cell")
    max_field = max(max_field, max(abs(v-expected) for v in actual))
    if abs(sum(actual)/24 - inventories[-1][0 if field == "HT" else 1]) > tol:
        fail(f"{field} full-field vs volume inventory mismatch")

if max_field > tol:
    fail(f"analytical backward-Euler field error {max_field:.9g} exceeds {tol}")
if max_budget > tol:
    fail(f"independent species budget error {max_budget:.9g} exceeds {tol}")
if max_total > 2*tol:
    fail(f"combined activity budget error {max_total:.9g} exceeds {2*tol}")
print(f"PASS: M1C k={k} lambda={lam} steps={n} cells=24 "
      f"max_BE_error={max_field:.12g} Bq/m3 "
      f"max_species_budget={max_budget:.12g} Bq "
      f"max_total_budget={max_total:.12g} Bq "
      f"final_HT={ht:.12g} final_HTO={hto:.12g}")
