#!/usr/bin/env python3
"""Independent continuum traveling-wave and discrete FVM budget audit."""
from pathlib import Path
import math
import re
import sys

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VEC = re.compile(r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)")
SCALAR = re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", re.S)
VECTOR = re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", re.S)
FLUX = re.compile(r"CTEAS_M1A scalar_flux_xMin=(" + NUM + r") scalar_flux_xMax=(" + NUM + r") scalar_flux_net=(" + NUM + r") scalar_budget_residual=(" + NUM + r")")
MASS = re.compile(r"CTEAS_M1A integral=(" + NUM + r")")
INITIAL = re.compile(r"CTEAS_M1A initial_integral=(" + NUM + r")")


def field(path, n, vector=False):
    data = path.read_text()
    m = (VECTOR if vector else SCALAR).search(data)
    if m:
        if int(m.group(1)) != n:
            raise SystemExit(f"FAIL: wrong native field length {path}")
        values = ([tuple(map(float, x)) for x in VEC.findall(m.group(2))] if vector
                  else [float(x) for x in m.group(2).split()])
    else:
        pattern = (r"internalField\s+uniform\s+\(([^()]*)\)\s*;" if vector else
                   r"internalField\s+uniform\s+(" + NUM + r")\s*;")
        m = re.search(pattern, data)
        if not m:
            raise SystemExit(f"FAIL: native field not found {path}")
        val = tuple(map(float, m.group(1).split())) if vector else float(m.group(1))
        values = [val]*n
    if len(values) != n or (vector and any(len(v)!=3 for v in values)):
        raise SystemExit(f"FAIL: native field malformed {path}")
    if not all(math.isfinite(x) for v in values for x in (v if vector else (v,))):
        raise SystemExit(f"FAIL: nonfinite native field {path}")
    return values


def mesh_audit(case, meshlog, dims):
    nx,ny,nz = dims
    n = nx*ny*nz
    points_data = (case/"constant/polyMesh/points").read_text()
    m = re.search(r"(?m)^\s*(\d+)\s*\n\s*\(\s*\n", points_data)
    if not m:
        raise SystemExit("FAIL: actual native mesh points missing")
    tail = points_data[m.end():]
    end = re.search(r"(?m)^\s*\)\s*$",tail)
    if not end:
        raise SystemExit("FAIL: actual native mesh terminator missing")
    points = [tuple(map(float,p)) for p in VEC.findall(tail[:end.start()])]
    if len(points) != int(m.group(1)) or len(points) != (nx+1)*(ny+1)*(nz+1):
        raise SystemExit("FAIL: actual native mesh vertex count")
    moved = 0
    for x,y,z in points:
        ix,iy,iz = round(nx*x),round(ny*y),round(nz*z)
        if not (0<=ix<=nx and 0<=iy<=ny and 0<=iz<=nz):
            raise SystemExit("FAIL: vertex outside cube")
        bx,by,bz = ix/nx,iy/ny,iz/nz
        wall = ix in (0,nx) or iy in (0,ny) or iz in (0,nz)
        shift = 0 if wall else 0.04*math.sin(math.pi*bx)*math.sin(math.pi*by)*math.sin(math.pi*bz)
        if abs(x-bx-shift)>1e-9 or abs(y-by)>1e-9 or abs(z-bz)>1e-9:
            raise SystemExit("FAIL: wrong physical mesh warp or changed outer boundary")
        moved += abs(shift)>1e-12
    if moved != (nx-1)*(ny-1)*(nz-1):
        raise SystemExit("FAIL: interior vertices not all warped")
    data = meshlog.read_text()
    if "Mesh OK." not in data:
        raise SystemExit("FAIL: checkMesh not OK")
    def metric(pat):
        found = re.search(pat,data,re.I)
        if not found:
            raise SystemExit(f"FAIL: missing checkMesh metric {pat}")
        return float(found.group(1).rstrip('.'))
    nonorth = metric(r"Mesh non-orthogonality Max:\s*("+NUM+")")
    skew = metric(r"Max skewness\s*=\s*("+NUM+")")
    volume = metric(r"Total volume\s*=\s*("+NUM+")")
    minvol = metric(r"Min volume\s*=\s*("+NUM+")")
    if not (1<nonorth<65 and 0<=skew<3.5 and abs(volume-1)<1e-6 and minvol>0):
        raise SystemExit("FAIL: native mesh geometry metrics")
    boundary = (case/"constant/polyMesh/boundary").read_text()
    for patch in ("xMin","xMax"):
        if not re.search(rf"\b{patch}\s*\{{\s*type\s+cyclic;",boundary):
            raise SystemExit(f"FAIL: {patch} is not a real cyclic polyMesh patch")
    for patch in ("yMin","yMax","zMin","zMax"):
        if not re.search(rf"\b{patch}\s*\{{\s*type\s+wall;",boundary):
            raise SystemExit(f"FAIL: {patch} not a real physical wall")
    return n,moved,nonorth


PISO = re.compile(r"CTEAS_M1B1 max_div_phi=(" + NUM + r") boundary_net_flux=(" + NUM + r") kinetic_energy=(" + NUM + r")")


def audit(case, meshlog, solverlog, dims, max_linf, max_rms):
    n, moved, nonorth = mesh_audit(case, meshlog, dims)
    props = (case / "constant/transportProperties").read_text()
    scheme = (case / "system/fvSchemes").read_text()
    control = (case / "system/controlDict").read_text()
    required = (
        "D [0 2 -1 0 0 0 0] 0.01;",
        "nu [0 2 -1 0 0 0 0] 0.05;",
        "momentumMode pisoLaminar;",
        "drivingAcceleration [0 1 -2 0 0 0 0] (0.5 0 0);",
        "nPressureCorrectors 3;",
        "nNonOrthogonalCorrectors 2;",
        "auditPassiveScalarFlux true;",
    )
    if any(x not in props for x in required):
        raise SystemExit("FAIL: PISO acceleration/diffusion coefficient contract")
    if not all(x in scheme for x in (
        "gradSchemes { default leastSquares; }",
        "div(phi,U) Gauss upwind;",
        "div(phi,c) Gauss upwind;",
        "laplacianSchemes { default Gauss linear corrected; }",
        "fluxRequired { default no; p; c; }",
    )):
        raise SystemExit("FAIL: unexpected PISO/scalar discretization")
    if not all(x in control for x in ("deltaT 0.01;", "endTime 0.4;", "writeInterval 40;")):
        raise SystemExit("FAIL: PISO time integration contract")
    for name in ("U", "p", "c"):
        data = (case / "0" / name).read_text()
        for patch in ("xMin", "xMax"):
            if not re.search(rf"\b{patch}\s*\{{\s*type\s+cyclic;", data):
                raise SystemExit(f"FAIL: {name}/{patch} is not cyclic")
    for patch in ("yMin", "yMax", "zMin", "zMax"):
        if not re.search(rf"\b{patch}\s*\{{\s*type\s+slip;", (case/"0/U").read_text()):
            raise SystemExit("FAIL: expected physical slip walls for U")
        if not re.search(rf"\b{patch}\s*\{{\s*type\s+zeroGradient;", (case/"0/c").read_text()):
            raise SystemExit("FAIL: expected insulating scalar walls")
    centres = field(case / "0/C", n, True)
    c0 = field(case / "0/c", n)
    cf = field(case / "0.4/c", n)
    u0 = field(case / "0/U", n, True)
    uf = field(case / "0.4/U", n, True)
    pf = field(case / "0.4/p", n)
    if max(abs(v) for row in u0 for v in row) > 1e-12:
        raise SystemExit("FAIL: initial PISO velocity not zero")
    uerr = max(abs(value - target) for row in uf for value, target in zip(row, (0.2, 0., 0.)))
    if uerr > 2e-5:
        raise SystemExit(f"FAIL: PISO-computed final U error {uerr:.12g}")
    prange = max(pf) - min(pf)
    if prange > 2e-5:
        raise SystemExit(f"FAIL: PISO pressure not spatially constant {prange:.12g}")
    init_ref = [0.5 + 0.25*math.cos(2*math.pi*x) for x,y,z in centres]
    final_ref = [0.5 + 0.25*math.exp(-0.01*(2*math.pi)**2*0.4)*
                 math.cos(2*math.pi*(x-0.25*0.4**2)) for x,y,z in centres]
    init_error = max(abs(a-b) for a,b in zip(c0, init_ref))
    errors = [abs(a-b) for a,b in zip(cf, final_ref)]
    linf = max(errors)
    rms = math.sqrt(sum(e*e for e in errors)/n)
    if init_error > 1e-12:
        raise SystemExit(f"FAIL: PISO scalar native-centre initialization {init_error:.12g}")
    if linf > max_linf or rms > max_rms:
        raise SystemExit(f"FAIL: PISO scalar analytical error Linf={linf:.12g} RMS={rms:.12g}")
    if not all(0.23 <= x <= 0.77 for x in cf):
        raise SystemExit("FAIL: PISO scalar maximum principle")
    log = solverlog.read_text()
    if "mode=pisoLaminar" not in log or len(re.findall(r"(?m)^Time = ", log)) != 40:
        raise SystemExit("FAIL: missing 40 real PISO time steps")
    if (any(len(re.findall(rf"Solving for U{axis},", log)) != 40
            for axis in ("x", "y", "z")) or
        len(re.findall(r"Solving for c,", log)) != 40 or
        len(re.findall(r"Solving for p,", log)) != 360):
        raise SystemExit("FAIL: incomplete native momentum/pressure/scalar solves")
    p_records = [tuple(map(float, x)) for x in PISO.findall(log)]
    if len(p_records) != 40:
        raise SystemExit("FAIL: incomplete PISO pressure/flux/energy history")
    maxdiv = max(abs(row[0]) for row in p_records)
    maxnet = max(abs(row[1]) for row in p_records)
    maxenergy = max(abs(row[2] - 0.5*(0.5*(i*0.01))**2)
                    for i, row in enumerate(p_records, 1))
    if maxdiv > 1e-8 or maxnet > 1e-8 or maxenergy > 2e-5:
        raise SystemExit(f"FAIL: PISO momentum energy/continuity maxdiv={maxdiv:.12g} "
                         f"maxnet={maxnet:.12g} energy={maxenergy:.12g}")
    fluxes = [tuple(map(float, row)) for row in FLUX.findall(log)]
    masses = [float(x) for x in MASS.findall(log)]
    starts = [float(x) for x in INITIAL.findall(log)]
    if len(fluxes) != 40 or len(masses) != 40 or len(starts) != 1:
        raise SystemExit("FAIL: incomplete PISO scalar inventory history")
    if not all(math.isfinite(x) for x in [*starts,*masses,
                *(v for row in fluxes for v in row),
                *(v for row in p_records for v in row)]):
        raise SystemExit("FAIL: nonfinite native PISO or scalar log")
    if abs(starts[0] - 0.5) > 0.02:
        raise SystemExit("FAIL: unexpected initial PISO scalar inventory")
    previous = starts[0]
    worst = 0.
    mincyclic = float('inf')
    for i, (mass, (left,right,net,residual)) in enumerate(zip(masses,fluxes), 1):
        error = abs(mass - previous + 0.01*net)
        worst = max(worst,error,abs(residual))
        mincyclic = min(mincyclic,abs(left),abs(right))
        if abs(left+right) > 1e-8 or abs(net) > 1e-8:
            raise SystemExit(f"FAIL: PISO cyclic scalar flux mismatch step={i}")
        if min(abs(left),abs(right)) <= 1e-4:
            raise SystemExit(f"FAIL: PISO cyclic scalar flux missing step={i}")
        if error > 2e-6 or abs(residual) > 1e-7 or abs(mass-starts[0]) > 1e-6:
            raise SystemExit(f"FAIL: PISO scalar conservation step={i}")
        previous = mass
    print(f"M1B-1g-11 cells={n} moved={moved} nonorth_deg={nonorth:.9g} "
          f"final_U_Linf={uerr:.12g} p_range={prange:.12g} "
          f"max_div_phi={maxdiv:.12g} max_energy_error={maxenergy:.12g} "
          f"scalar_Linf={linf:.12g} scalar_RMS={rms:.12g} "
          f"min_cyclic_flux={mincyclic:.12g} max_scalar_budget={worst:.12g}")
    return rms


if len(sys.argv) != 7:
    raise SystemExit("usage: check_warped_piso_accelerating_advdiff.py COARSE_CASE COARSE_MESH COARSE_SOLVE FINE_CASE FINE_MESH FINE_SOLVE")
coarse = audit(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), (8,20,6), 0.055, 0.035)
fine = audit(Path(sys.argv[4]), Path(sys.argv[5]), Path(sys.argv[6]), (12,30,9), 0.040, 0.025)
ratio = fine/coarse if coarse > 0 else float('inf')
print(f"M1B-1g-11 refinement_RMS_ratio={ratio:.12g}")
if not (0 < fine and 0 < coarse and ratio <= 0.85):
    raise SystemExit("FAIL: PISO scalar spatial refinement")
print("PASS: M1B-1g-11 PISO-computed accelerating velocity and scalar analytic reference")
