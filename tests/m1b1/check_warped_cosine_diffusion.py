#!/usr/bin/env python3
"""Independent non-affine transient Neumann reference and FVM conservation audit."""
from pathlib import Path
import math
import re
import sys

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VECTOR = re.compile(r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)")
FLUX = re.compile(
    r"CTEAS_M1A scalar_flux_xMin=(" + NUM + r") scalar_flux_xMax=(" + NUM +
    r") scalar_flux_net=(" + NUM + r") scalar_budget_residual=(" + NUM + r")"
)
MASS = re.compile(r"CTEAS_M1A integral=(" + NUM + r")")
START = re.compile(r"CTEAS_M1A initial_integral=(" + NUM + r")")
FLOW = re.compile(r"CTEAS_M1B1 pressure_inlet_flux=(" + NUM + r") pressure_outlet_flux=(" + NUM + r")")
BUDGET = re.compile(r"CTEAS_M1B1 max_div_phi=(" + NUM + r") boundary_net_flux=(" + NUM + r") kinetic_energy=(" + NUM + r")")
SCA = re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", re.S)
VEC = re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", re.S)


def field(path, n, vector=False):
    data = path.read_text()
    match = (VEC if vector else SCA).search(data)
    if match:
        if int(match.group(1)) != n:
            raise SystemExit(f"FAIL: wrong native field count in {path}")
        vals = ([tuple(map(float, v)) for v in VECTOR.findall(match.group(2))]
                if vector else [float(x) for x in match.group(2).split()])
    else:
        pattern = (r"internalField\s+uniform\s+\(([^()]*)\)\s*;" if vector else
                   r"internalField\s+uniform\s+(" + NUM + r")\s*;")
        match = re.search(pattern, data)
        if not match:
            raise SystemExit(f"FAIL: native field missing {path}")
        one = (tuple(map(float, match.group(1).split())) if vector else float(match.group(1)))
        vals = [one] * n
    if len(vals) != n or (vector and any(len(v) != 3 for v in vals)):
        raise SystemExit(f"FAIL: native field malformed {path}")
    if not all(math.isfinite(x) for v in vals for x in (v if vector else (v,))):
        raise SystemExit(f"FAIL: nonfinite native field {path}")
    return vals


def geometry(case, dims, meshlog):
    nx, ny, nz = dims
    n = nx * ny * nz
    data = (case / "constant/polyMesh/points").read_text()
    m = re.search(r"(?m)^\s*(\d+)\s*\n\s*\(\s*\n", data)
    if not m:
        raise SystemExit("FAIL: real polyMesh points missing")
    tail = data[m.end():]
    close = re.search(r"(?m)^\s*\)\s*$", tail)
    if not close:
        raise SystemExit("FAIL: native points terminator missing")
    points = [tuple(map(float, p)) for p in VECTOR.findall(tail[:close.start()])]
    if len(points) != int(m.group(1)) or len(points) != (nx+1)*(ny+1)*(nz+1):
        raise SystemExit("FAIL: actual 3D point count")
    moved = 0
    for x,y,z in points:
        ix,iy,iz = round(nx*x),round(ny*y),round(nz*z)
        if not (0 <= ix <= nx and 0 <= iy <= ny and 0 <= iz <= nz):
            raise SystemExit("FAIL: physical point outside cube")
        bx,by,bz = ix/nx,iy/ny,iz/nz
        wall = ix in (0,nx) or iy in (0,ny) or iz in (0,nz)
        shift = 0 if wall else 0.04*math.sin(math.pi*bx)*math.sin(math.pi*by)*math.sin(math.pi*bz)
        if abs(x-bx-shift)>1e-9 or abs(y-by)>1e-9 or abs(z-bz)>1e-9:
            raise SystemExit("FAIL: unexpected native warp or moved physical boundary")
        moved += abs(shift)>1e-12
    if moved != (nx-1)*(ny-1)*(nz-1):
        raise SystemExit("FAIL: not all interior vertices warped")
    log = meshlog.read_text()
    if "Mesh OK." not in log:
        raise SystemExit("FAIL: official checkMesh not OK")
    def metric(pattern):
        m = re.search(pattern, log, re.I)
        if not m:
            raise SystemExit("FAIL: checkMesh metric missing")
        return float(m.group(1).rstrip("."))
    nonorth = metric(r"Mesh non-orthogonality Max:\s*(" + NUM + ")")
    skew = metric(r"Max skewness\s*=\s*(" + NUM + ")")
    minvol = metric(r"Min volume\s*=\s*(" + NUM + ")")
    total = metric(r"Total volume\s*=\s*(" + NUM + ")")
    if not (1 < nonorth < 65 and 0 <= skew < 3.5 and minvol > 0 and abs(total-1)<1e-6):
        raise SystemExit("FAIL: mesh quality or physical volume")
    return n, moved, nonorth


def audit(case, meshlog, solverlog, dims):
    n, moved, nonorth = geometry(case, dims, meshlog)
    c0 = (case / "0/c").read_text()
    props = (case / "constant/transportProperties").read_text()
    schemes = (case / "system/fvSchemes").read_text()
    if len(re.findall(r"(?m)^D \[0 2 -1 0 0 0 0\] 0\.01;$", props)) != 1:
        raise SystemExit("FAIL: actual diffusivity not D=0.01")
    if props.count("auditPassiveScalarFlux true;") != 1:
        raise SystemExit("FAIL: actual finite-volume flux audit disabled")
    if ("gradSchemes { default leastSquares; }" not in schemes or
        "laplacianSchemes { default Gauss linear corrected; }" not in schemes or
        "fluxRequired { default no; p; c; }" not in schemes):
        raise SystemExit("FAIL: unexpected diffusion discretization")
    for patch in ("xMin","xMax","yMin","yMax","zMin","zMax"):
        if not re.search(rf"\b{patch}\s*\{{\s*type\s+zeroGradient;", c0):
            raise SystemExit(f"FAIL: {patch} not insulating")
    for patch in ("xMin","xMax"):
        if not re.search(rf"\b{patch}\s*\{{\s*type\s+fixedValue;\s*value\s+uniform\s+0;", (case/"0/p").read_text()):
            raise SystemExit(f"FAIL: {patch} has nonzero pressure")
    centres = field(case/"0/C", n, True)
    seed = field(case/"0/c", n)
    final = field(case/"0.4/c", n)
    U = field(case/"0.4/U", n, True)
    p = field(case/"0.4/p", n)
    initial_ref = [0.5+0.25*math.cos(math.pi*x) for x,y,z in centres]
    ref = [0.5+0.25*math.exp(-0.01*math.pi**2*0.4)*math.cos(math.pi*x) for x,y,z in centres]
    init_err = max(abs(a-b) for a,b in zip(seed,initial_ref))
    errors = [abs(a-b) for a,b in zip(final,ref)]
    linf, rms = max(errors), math.sqrt(sum(e*e for e in errors)/n)
    if init_err>1e-12:
        raise SystemExit(f"FAIL: wrong native-centre initial cosine {init_err}")
    if linf>0.004 or rms>0.002:
        raise SystemExit(f"FAIL: eigenmode analytical error Linf={linf} RMS={rms}")
    if not all(0.24 <= c <= 0.76 for c in final):
        raise SystemExit("FAIL: cosine scalar bounds violated")
    if max(abs(v) for u in U for v in u)>1e-8 or max(abs(v) for v in p)>1e-8:
        raise SystemExit("FAIL: spurious zero-drive momentum or pressure")
    log = solverlog.read_text()
    if "mode=pisoLaminar" not in log or len(re.findall(r"(?m)^CTEAS_M1B1 nonorthogonal_correctors=2$",log))!=20:
        raise SystemExit("FAIL: incomplete 20 PISO steps")
    if len(re.findall(r"Solving for p,",log)) != 180:
        raise SystemExit("FAIL: incomplete 180 real pressure solves")
    f = [tuple(map(float,row)) for row in FLUX.findall(log)]
    q = [tuple(map(float,row)) for row in FLOW.findall(log)]
    b = [tuple(map(float,row)) for row in BUDGET.findall(log)]
    mass = [float(v) for v in MASS.findall(log)]
    initial = [float(v) for v in START.findall(log)]
    if any(len(v)!=20 for v in (f,q,b,mass)) or len(initial)!=1:
        raise SystemExit("FAIL: incomplete scalar/flow finite-volume budget history")
    if not all(math.isfinite(x) for x in [*initial,*mass,*(x for rows in (f,q,b) for row in rows for x in row)]):
        raise SystemExit("FAIL: nonfinite logged FVM budgets")
    if abs(initial[0]-0.5)>0.01:
        raise SystemExit("FAIL: wrong initial volume inventory")
    old=initial[0]
    worst=0
    for step,(new,flux,flow,budget) in enumerate(zip(mass,f,q,b),1):
        left,right,net,residual=flux
        err=abs(new-old+0.02*net)
        worst=max(worst,err,abs(residual))
        if max(abs(left),abs(right),abs(net))>1e-8:
            raise SystemExit(f"FAIL: nonzero Neumann scalar boundary flux at step {step}")
        if err>2e-6 or abs(residual)>1e-7 or abs(new-initial[0])>1e-6:
            raise SystemExit(f"FAIL: nonconservative transient scalar inventory at step {step}")
        if max(abs(v) for v in flow)>1e-8 or max(abs(budget[0]),abs(budget[1]))>1e-7:
            raise SystemExit(f"FAIL: spurious flow/continuity at step {step}")
        old=new
    print(f"M1B-1g-9 cells={n} moved={moved} nonorth_deg={nonorth:.9g} "
          f"initial_Linf={init_err:.12g} final_Linf={linf:.12g} "
          f"final_RMS={rms:.12g} initial_inventory={initial[0]:.12g} "
          f"max_mass_budget={worst:.12g}")
    return rms


if len(sys.argv)!=7:
    raise SystemExit("usage: check_warped_cosine_diffusion.py COARSE_CASE COARSE_MESH COARSE_SOLVE FINE_CASE FINE_MESH FINE_SOLVE")
coarse = audit(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),(8,20,6))
fine = audit(Path(sys.argv[4]),Path(sys.argv[5]),Path(sys.argv[6]),(12,30,9))
ratio = fine/coarse if coarse>0 else float("inf")
print(f"M1B-1g-9 refinement_RMS_ratio={ratio:.12g}")
if not (0 < fine and 0 < coarse and ratio<=0.85):
    raise SystemExit("FAIL: non-affine transient spatial refinement criterion")
print("PASS: official warped 3D non-affine transient Neumann diffusion, conservation and refinement")
