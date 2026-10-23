#!/usr/bin/env python3
"""M1B-1g-3: native-geometry, three-resolution spatial convergence audit.

All three meshes preserve identical physical boundaries, pressure drop,
viscosity, timestep and stop time. Error versus the independent analytic
square-duct Stokes reference is evaluated on each mesh's actual OpenFOAM
finite-volume centroids, not index-based or interpolated coordinates.
"""
import math
from pathlib import Path
import re
import sys

from check_warped_open_pressure import (
    BUDGET, FLOW, SCALAR, VEC, SCA, VEC3, UNIFORM_V, UNIFORM_S,
    duct_velocity, quality,
)

SPECS = ((8, 20, 6), (12, 30, 9), (16, 40, 12))
DT, END, NSTEPS, PSOLVES = 0.01, 0.4, 40, 360


def field(case, time, name, vector, expected_cells):
    path = Path(case) / time / name
    if not path.is_file():
        raise SystemExit(f"FAIL: missing actual OpenFOAM field {path}")
    data = path.read_text()
    match = (VEC if vector else SCA).search(data)
    if match:
        if int(match.group(1)) != expected_cells:
            raise SystemExit(f"FAIL: wrong native field cell count in {path}")
        values = ([tuple(map(float, v)) for v in VEC3.findall(match.group(2))]
                  if vector else [float(x) for x in match.group(2).split()])
    else:
        match = (UNIFORM_V if vector else UNIFORM_S).search(data)
        if not match:
            raise SystemExit(f"FAIL: missing field internalField in {path}")
        scalar = (tuple(map(float, match.group(1).split()))
                  if vector else float(match.group(1)))
        values = [scalar] * expected_cells
    if len(values) != expected_cells:
        raise SystemExit(f"FAIL: truncated native field in {path}")
    if any(not math.isfinite(x) for v in values for x in (v if vector else (v,))):
        raise SystemExit(f"FAIL: nonfinite native field in {path}")
    return values


def true_geometry(case, nx, ny, nz):
    """Match each actual vertex against the analytic warped lattice as a set."""
    path = Path(case) / "constant/polyMesh/points"
    data = path.read_text()
    m = re.search(r"(?m)^\s*(\d+)\s*\n\s*\(\s*\n", data)
    if not m:
        raise SystemExit("FAIL: missing OpenFOAM physical vertices")
    tail = data[m.end():]
    end = re.search(r"(?m)^\s*\)\s*$", tail)
    if not end:
        raise SystemExit("FAIL: unterminated physical points file")
    actual = [tuple(map(float, v)) for v in VEC3.findall(tail[:end.start()])]
    expected = []
    for iz in range(nz+1):
        z = iz / nz
        for iy in range(ny+1):
            y = iy / ny
            for ix in range(nx+1):
                x = ix / nx
                boundary = ix in (0,nx) or iy in (0,ny) or iz in (0,nz)
                warp = (0 if boundary else 0.04*math.sin(math.pi*x)
                        *math.sin(math.pi*y)*math.sin(math.pi*z))
                expected.append((x+warp, y, z))
    if len(actual) != len(expected) or len(actual) != int(m.group(1)):
        raise SystemExit("FAIL: actual warped vertex count is wrong")
    # Check an unordered geometric point *multiset*. Sorting raw floats by
    # x can permute mirror-symmetric (y,z) points when their x values
    # differ by ~1e-16 in the native OpenFOAM ASCII representation.
    # Quantizing every coordinate to 9 decimals is a stricter ≤5e-10 m
    # mesh comparison than the original 1e-9 m coordinate tolerance.
    from collections import Counter
    key = lambda xyz: tuple(round(value, 9) for value in xyz)
    if Counter(map(key, actual)) != Counter(map(key, expected)):
        raise SystemExit("FAIL: actual warped geometry differs from expected physical lattice")
    moved = (nx-1)*(ny-1)*(nz-1)
    if moved < 100:
        raise SystemExit("FAIL: insufficient genuine 3-D interior deformation")
    return moved


def run_one(case, meshlog, solverlog, dims):
    nx,ny,nz = dims
    n = nx*ny*nz
    moved = true_geometry(case,nx,ny,nz)
    nonorth, skew, vol = quality(meshlog, warped=True)
    text = Path(solverlog).read_text()
    pressure = [tuple(map(float, q)) for q in FLOW.findall(text)]
    budget = [tuple(map(float, q)) for q in BUDGET.findall(text)]
    inventories = [float(x) for x in SCALAR.findall(text)]
    nsolves = len(re.findall(r"Solving for p,", text))
    ncorr = len(re.findall(r"^CTEAS_M1B1 nonorthogonal_correctors=2$",text,re.M))
    if not ("mode=pisoLaminar" in text and "Completed mode=pisoLaminar" in text):
        raise SystemExit(f"FAIL: real PISO laminar solve incomplete at n={n}")
    if (len(pressure)!=NSTEPS or len(budget)!=NSTEPS
            or len(inventories)!=NSTEPS or ncorr!=NSTEPS or nsolves!=PSOLVES):
        raise SystemExit(f"FAIL: partial 40-step and 360-solve history at n={n}")
    if any(not math.isfinite(q) for rows in (pressure,budget)
           for row in rows for q in row):
        raise SystemExit("FAIL: nonfinite signed flux or cellwise conservation")
    C = field(case,"0","C",True,n)
    U0 = field(case,"0","U",True,n)
    p0 = field(case,"0","p",False,n)
    c0 = field(case,"0","c",False,n)
    U = field(case,"0.4","U",True,n)
    p = field(case,"0.4","p",False,n)
    c = field(case,"0.4","c",False,n)
    if max(abs(v) for u in U0 for v in u) > 1e-12:
        raise SystemExit("FAIL: initial flow was not rest")
    if max(abs(v-0.1*(1-xyz[0])) for v,xyz in zip(p0,C)) > 5e-12:
        raise SystemExit("FAIL: initial pressure not seeded at native OpenFOAM centroids")
    if max(abs(v-1) for v in c0) > 1e-12:
        raise SystemExit("FAIL: initial scalar not unit uniform")
    control = (Path(case)/"system/controlDict").read_text()
    if not re.search(r"\bdeltaT\s+0\.01\s*;",control) or not re.search(r"\bendTime\s+0\.4\s*;",control):
        raise SystemExit("FAIL: spatial study changed physical integration interval")
    tp = (Path(case)/"0.4/p").read_text()
    tu = (Path(case)/"0.4/U").read_text()
    for patch,dp in (("xMin",0.1),("xMax",0.0)):
        m = re.search(rf"\b{patch}\s*\{{\s*type\s+fixedValue;\s*value\s+uniform\s+([-+0-9.eE]+)\s*;",tp)
        if not m or abs(float(m.group(1))-dp)>1e-12:
            raise SystemExit(f"FAIL: changed Dirichlet kinematic pressure at {patch}")
    for patch in ("yMin","yMax","zMin","zMax"):
        if not re.search(rf"\b{patch}\s*\{{\s*type\s+fixedValue;\s*value\s+uniform\s+\(0 0 0\)",tu):
            raise SystemExit(f"FAIL: missing no-slip wall {patch}")
    errs = [u[0]-duct_velocity(xyz[1],xyz[2]) for u,xyz in zip(U,C)]
    l2 = math.sqrt(math.fsum(e*e for e in errs)/n)
    linf = max(abs(e) for e in errs)
    transverse = max(abs(v) for u in U for v in u[1:])
    p_err = max(abs(v-0.1*(1-xyz[0])) for v,xyz in zip(p,C))
    divergence = max(abs(row[0]) for row in budget)
    boundary = max(abs(row[1]) for row in budget)
    fluxpair = max(abs(a+b) for a,b in pressure)
    qin,qout = pressure[-1]
    drift = max(abs(v-1) for v in inventories)
    c_err = max(abs(v-1) for v in c)
    print(f"M1B-1g-3 dims={nx}x{ny}x{nz} cells={n} moved={moved} "
          f"nonorth_deg={nonorth:.12g} skew={skew:.12g} minvol={vol:.12g} "
          f"p_solves={nsolves} U_error_L2={l2:.12g} U_error_Linf={linf:.12g} "
          f"pressure_error_Linf={p_err:.12g} transverse_U={transverse:.12g} "
          f"q_in={qin:.12g} q_out={qout:.12g} "
          f"max_div_phi={divergence:.12g} net_boundary={boundary:.12g} "
          f"max_pair_flux={fluxpair:.12g} inventory_drift={drift:.12g}")
    if (divergence>1e-7 or boundary>1e-7 or fluxpair>1e-7 or
            drift>1e-8 or c_err>1e-8 or qin>=-0.003 or qout<=0.003 or
            linf>0.025 or p_err>0.03 or transverse>0.02):
        raise SystemExit(f"FAIL: finite-volume analytical/flux gate at n={n}")
    return l2,linf,qout


if len(sys.argv)!=10:
    raise SystemExit("usage: check_warped_open_pressure_spatial.py CASE1 CHECK1 LOG1 CASE2 CHECK2 LOG2 CASE3 CHECK3 LOG3")
runs = [run_one(*sys.argv[1+3*i:4+3*i],SPECS[i]) for i in range(3)]
errs = [r[0] for r in runs]
ratios = [errs[i+1]/errs[i] if errs[i]>1e-12 else math.inf for i in range(2)]
fluxes = [r[2] for r in runs]
print(f"M1B-1g-3 native_L2_errors={errs} successive_ratios={ratios} "
      f"coarse_to_fine={errs[-1]/errs[0]:.12g} outlet_Q={fluxes}")
# Predeclared deliberately modest gate: decrease at least 10% each level
# and at least 30% over two doublings; no claim of measured 2nd order.
if any(r>=0.9 for r in ratios) or errs[-1]/errs[0]>0.7:
    raise SystemExit("FAIL: genuine warped 3-D spatial refinement did not decrease Stokes L2 error")
print("PASS: native-centroid physical 3D warped duct spatial L2 error convergence")
