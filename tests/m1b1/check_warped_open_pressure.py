#!/usr/bin/env python3
"""M1B-1g independently audit actual 3-D nonperiodic, distorted duct PDE.

The analytic reference solves transient square-duct laminar Stokes flow with
constant pressure gradient, initial rest and no-slip y/z walls. The numerical
Navier-Stokes solver includes convection; a relatively coarse 8x20x6 grid
and static-p two-end BC mean this is a bounded verification, not a full
building/terrain validation.
"""
from pathlib import Path
import math
import re
import sys

N = 8*20*6
VEC = re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", re.S)
SCA = re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", re.S)
VEC3 = re.compile(r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)")
UNIFORM_V = re.compile(r"internalField\s+uniform\s+\(([^()]*)\)\s*;")
UNIFORM_S = re.compile(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;")
FLOW = re.compile(r"CTEAS_M1B1 pressure_inlet_flux=([-+0-9.eE]+) pressure_outlet_flux=([-+0-9.eE]+)")
BUDGET = re.compile(r"CTEAS_M1B1 max_div_phi=([-+0-9.eE]+) boundary_net_flux=([-+0-9.eE]+) kinetic_energy=([-+0-9.eE]+)")
SCALAR = re.compile(r"CTEAS_M1A integral=([-+0-9.eE]+)")


def read_field(case, time, name, vector):
    file = Path(case) / time / name
    if not file.is_file():
        raise SystemExit(f"FAIL: missing actual OpenFOAM {file}")
    txt = file.read_text()
    found = (VEC if vector else SCA).search(txt)
    if found:
        if int(found.group(1)) != N:
            raise SystemExit(f"FAIL: wrong field cell count in {file}")
        values = ([tuple(map(float, x)) for x in VEC3.findall(found.group(2))]
                  if vector else [float(x) for x in found.group(2).split()])
    else:
        uni = (UNIFORM_V if vector else UNIFORM_S).search(txt)
        if not uni:
            raise SystemExit(f"FAIL: no actual internalField in {file}")
        v = (tuple(map(float, uni.group(1).split())) if vector else float(uni.group(1)))
        values = [v]*N
    if len(values) != N or (vector and any(len(v) != 3 for v in values)):
        raise SystemExit(f"FAIL: malformed actual field in {file}")
    if not all(math.isfinite(x) for v in values for x in (v if vector else [v])):
        raise SystemExit(f"FAIL: nonfinite actual values in {file}")
    return values


def audit_physical_points(case, warped):
    """Independently verify real polyMesh points and six unchanged physical faces."""
    file = Path(case)/"constant"/"polyMesh"/"points"
    data = file.read_text()
    m = re.search(r"(?m)^\s*(\d+)\s*\n\s*\(\s*\n", data)
    if not m:
        raise SystemExit("FAIL: actual counted OpenFOAM mesh vertices missing")
    tail = data[m.end():]
    terminator = re.search(r"(?m)^\s*\)\s*$", tail)
    if not terminator:
        raise SystemExit("FAIL: invalid real polyMesh vertices terminator")
    pts = [tuple(map(float, v)) for v in VEC3.findall(tail[:terminator.start()])]
    if len(pts) != int(m.group(1)) or len(pts) != 9*21*7:
        raise SystemExit("FAIL: actual 3D mesh point count mismatch")
    moved = 0
    for x,y,z in pts:
        ix,iy,iz = round(8*x), round(20*y), round(6*z)
        # The x deformation is < 0.04m and Δx/2 = 0.0625m, hence
        # nearest original x grid plane is unique even after distortion.
        base = (ix/8, iy/20, iz/6)
        if any(abs(v-t) > 1e-9 for v,t in zip((y,z),base[1:])):
            raise SystemExit("FAIL: actual mesh y/z positions changed")
        if not (0 <= ix <= 8 and 0 <= iy <= 20 and 0 <= iz <= 6):
            raise SystemExit("FAIL: deformed mesh vertex outside physical cube")
        boundary = ix in (0,8) or iy in (0,20) or iz in (0,6)
        expected_dx = (0 if boundary or not warped else
                       .04*math.sin(math.pi*base[0])
                       *math.sin(math.pi*base[1])*math.sin(math.pi*base[2]))
        if abs(x-base[0]-expected_dx) > 1e-9:
            raise SystemExit("FAIL: physical wall/inlet/outlet altered or wrong real warp amplitude")
        moved += abs(x-base[0]) > 1e-10
    expected_moved = 7*19*5 if warped else 0
    if moved != expected_moved:
        raise SystemExit("FAIL: wrong number of truly deformed interior polyMesh points")
    return moved


def quality(log, warped):
    data = Path(log).read_text()
    if "Mesh OK." not in data:
        raise SystemExit("FAIL: official v2606 checkMesh did not accept topology")
    def number(pattern):
        m = re.search(pattern, data, re.I)
        if not m:
            raise SystemExit(f"FAIL: missing genuine checkMesh metric {pattern}")
        val = float(m.group(1).rstrip("."))
        if not math.isfinite(val):
            raise SystemExit("FAIL: nonfinite checkMesh metric")
        return val
    nonorth = number(r"Mesh non-orthogonality Max:\s*([-+0-9.eE]+)")
    skew = number(r"Max skewness\s*=\s*([-+0-9.eE]+)")
    minvol = number(r"Min volume\s*=\s*([-+0-9.eE]+)")
    total = number(r"Total volume\s*=\s*([-+0-9.eE]+)")
    if not minvol > 0 or abs(total-1) > 1e-6 or not 0 <= skew < 3.5:
        raise SystemExit("FAIL: non-positive cells / wrong volume / severe skew")
    if warped and not 1 < nonorth < 65:
        raise SystemExit("FAIL: controlled internal warp did not create acceptable nonorthogonality")
    if not warped and nonorth > 0.1:
        raise SystemExit("FAIL: nominally Cartesian control is not orthogonal")
    return nonorth, skew, minvol


def duct_velocity(y, z, t=0.4, nu=0.05, gradient=0.1):
    # u_t = nu*(u_yy+u_zz) + G; u(0)=0 and four no-slip walls.
    # Double sine-series coefficient of uniform drive is 16G/(pi^2*m*n).
    result = 0.0
    for m in range(1, 40, 2):
        sy = math.sin(m*math.pi*y)
        for n in range(1, 40, 2):
            lam = nu*math.pi**2*(m*m+n*n)
            result += (16*gradient/(math.pi**2*m*n)
                       *sy*math.sin(n*math.pi*z)
                       *(-math.expm1(-lam*t))/lam)
    return result


def audit(case, meshlog, solverlog, mode):
    warped = mode != "flat"
    moved = audit_physical_points(case, warped)
    nonorth, skew, minvol = quality(meshlog, warped)
    C = read_field(case, "0", "C", True)
    U = read_field(case, "0.4", "U", True)
    p = read_field(case, "0.4", "p", False)
    initial_p = read_field(case, "0", "p", False)
    logfile = Path(solverlog).read_text()
    pressure = [tuple(map(float, x)) for x in FLOW.findall(logfile)]
    budget = [tuple(map(float, x)) for x in BUDGET.findall(logfile)]
    scalar = [float(x) for x in SCALAR.findall(logfile)]
    if len(pressure) != 20 or len(budget) != 20 or len(scalar) != 20:
        raise SystemExit("FAIL: missing complete 20-step PDE/flux/scalar history")
    if ("mode=pisoLaminar" not in logfile or
        len(re.findall(r"^CTEAS_M1B1 nonorthogonal_correctors=2$", logfile, re.M)) != 20 or
        len(re.findall(r"Solving for p,", logfile)) != 180):
        raise SystemExit("FAIL: not 20 real PDE steps × 3 PISO × 3 nonorth pressure solves")
    if not all(math.isfinite(a) for row in pressure+budget for a in row):
        raise SystemExit("FAIL: nonfinite pressure-flow or continuity log")
    dp = 0 if mode == "zero" else (-.1 if mode == "reverse" else .1)
    pressuretext = (Path(case)/"0.4"/"p").read_text()
    for name, expected in (("xMin", dp), ("xMax", 0)):
        pat = rf"\b{name}\s*\{{\s*type\s+fixedValue;\s*value\s+uniform\s+([-+0-9.eE]+)\s*;"
        m = re.search(pat, pressuretext)
        if not m or abs(float(m.group(1))-expected) > 1e-12:
            raise SystemExit(f"FAIL: nonperiodic fixed kinematic pressure BC lost on {name}")
    walltext = (Path(case)/"0.4"/"U").read_text()
    for name in ("yMin", "yMax", "zMin", "zMax"):
        if not re.search(rf"\b{name}\s*\{{\s*type\s+fixedValue;\s*value\s+uniform\s+\(0 0 0\)", walltext):
            raise SystemExit(f"FAIL: no-slip wall {name} not preserved")
    dx = max(abs(x-((i%8)+.5)/8) for i,(x,y,z) in enumerate(C))
    if warped and dx < .01:
        raise SystemExit("FAIL: actual native OpenFOAM centres not warped")
    if not warped and dx > 1e-12:
        raise SystemExit("FAIL: flat control centres changed")
    initial_mismatch = max(abs(a-dp*(1-x)) for a,(x,y,z) in zip(initial_p,C))
    if initial_mismatch > 5e-12:
        raise SystemExit("FAIL: initial p not seeded at native OpenFOAM C")
    maxdiv = max(abs(row[0]) for row in budget)
    maxnet = max(abs(row[1]) for row in budget)
    scalar_drift = max(abs(x-1) for x in scalar)
    inq, outq = pressure[-1]
    flux_error = max(abs(row[0]+row[1]) for row in pressure)
    transverse = max(max(abs(v[1]), abs(v[2])) for v in U)
    xpressure = max(abs(pv-dp*(1-C[i][0])) for i,pv in enumerate(p))
    axialerror = max(abs(ux[0]-duct_velocity(y,z, gradient=dp))
                     for ux, (x,y,z) in zip(U,C))
    print(f"M1B-1g mode={mode} real3d=960 moved_interior_vertices={moved} nonorth_deg={nonorth:.12g} "
          f"skew={skew:.12g} minvol={minvol:.12g} max_center_dx={dx:.12g} "
          f"seed_p_err={initial_mismatch:.12g} q_in={inq:.12g} q_out={outq:.12g} "
          f"max_pair_flux_mismatch={flux_error:.12g} max_div_phi={maxdiv:.12g} "
          f"max_boundary_net_flux={maxnet:.12g} scalar_drift={scalar_drift:.12g} "
          f"analytical_duct_U_err={axialerror:.12g} linear_p_err={xpressure:.12g} "
          f"transverse_U={transverse:.12g}")
    if maxdiv > 1e-7 or maxnet > 1e-7 or flux_error > 1e-7 or scalar_drift > 1e-8:
        raise SystemExit("FAIL: actual 3D pressure-corrected conservative flux/scalar")
    if dp:
        if (dp > 0 and not (inq < -.003 and outq > .003)) or (dp < 0 and not (inq > .003 and outq < -.003)):
            raise SystemExit("FAIL: fixed-pressure gradient did not produce correctly signed throughflow")
        if axialerror > .025 or xpressure > .03 or transverse > .02:
            raise SystemExit("FAIL: actual 3D transient laminar duct analytical agreement")
    elif (max(abs(v) for u in U for v in u) > 1e-8 or
          max(abs(v) for v in p) > 1e-8 or max(abs(v) for pair in pressure for v in pair) > 1e-8):
        raise SystemExit("FAIL: zero-Δp warped 3D duct spontaneously generated flow")
    return U, pressure[-1]


if __name__ == "__main__":
    if len(sys.argv) != 10:
        raise SystemExit("usage: check_warped_open_pressure.py FLAT_CASE FLAT_CHECK FLAT_SOLVE "
                         "WARP_CASE WARP_CHECK WARP_SOLVE ZERO_CASE ZERO_CHECK ZERO_SOLVE")
    flat, fq = audit(*sys.argv[1:4], "flat")
    warped, wq = audit(*sys.argv[4:7], "warped")
    audit(*sys.argv[7:10], "zero")
    diff = max(abs(a-b) for f,w in zip(flat,warped) for a,b in zip(f,w))
    flow_rel = abs(fq[1]-wq[1])/max(abs(fq[1]),1e-12)
    print(f"M1B-1g flat_warp_U_Linf={diff:.12g} outlet_Q_relative_change={flow_rel:.12g}")
    if diff > .02 or flow_rel > .2:
        raise SystemExit("FAIL: interior mesh warp changed physically identical duct flow too much")
    print("PASS: true 3D nonperiodic fixed-pressure/4-wall duct, warped mesh, analytical transient Stokes, zero-drive control")
