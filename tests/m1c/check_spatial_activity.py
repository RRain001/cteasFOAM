#!/usr/bin/env python3
"""Independent continuum two-species traveling-wave and discrete FVM budget audit."""
import math
import re
import sys
from pathlib import Path

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
VECTOR = re.compile(r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)")
SCALAR = re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", re.S)
VECTORS = re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", re.S)
INVENTORY = re.compile(r"CTEAS_M1C activity_HT_Bq=(" + NUM + r") activity_HTO_Bq=(" + NUM + r")")
HT_FLUX = re.compile(r"CTEAS_M1C flux_HT_xMin_Bq_s=(" + NUM + r") flux_HT_xMax_Bq_s=(" + NUM + r") flux_HT_net_Bq_s=(" + NUM + r") flux_HT_budget_residual_Bq=(" + NUM + r")")
HTO_FLUX = re.compile(r"CTEAS_M1C flux_HTO_xMin_Bq_s=(" + NUM + r") flux_HTO_xMax_Bq_s=(" + NUM + r") flux_HTO_net_Bq_s=(" + NUM + r") flux_HTO_budget_residual_Bq=(" + NUM + r") flux_total_budget_residual_Bq=(" + NUM + r")")

def fail(msg):
    raise SystemExit("FAIL: M1C-1 " + msg)

def scalar_field(path, n, dims):
    data = path.read_text()
    if not re.search(r"dimensions\s+\[" + re.escape(dims) + r"\]", data):
        fail(f"wrong physical dimensions: {path}")
    m = SCALAR.search(data)
    if m:
        if int(m.group(1)) != n:
            fail(f"wrong native field size: {path}")
        values = [float(v) for v in m.group(2).split()]
    else:
        m = re.search(r"internalField\s+uniform\s+(" + NUM + r")\s*;", data)
        if not m:
            fail(f"unreadable native field: {path}")
        values = [float(m.group(1))]*n
    if len(values) != n or not all(math.isfinite(v) for v in values):
        fail(f"nonfinite/incomplete field: {path}")
    return values

def audit(case, meshlog, logpath, dims, limits):
    nx, ny, nz = dims
    n = nx*ny*nz
    mesh = meshlog.read_text()
    if "Mesh OK." not in mesh or not re.search(rf"(?m)^\s*cells:\s*{n}\s*$", mesh):
        fail("real checkMesh cell count or topology")
    nonorth = re.search(r"Mesh non-orthogonality Max:\s*(" + NUM + r")", mesh)
    if not nonorth or not 1 < float(nonorth.group(1)) < 65:
        fail("real mesh nonorthogonality")
    points = (case/"constant/polyMesh/points").read_text()
    m = re.search(r"(?m)^\s*(\d+)\s*\n\s*\(\s*\n", points)
    if not m:
        fail("native polyMesh points missing")
    end = re.search(r"(?m)^\s*\)\s*$", points[m.end():])
    if not end:
        fail("native polyMesh points malformed")
    xyz = [tuple(map(float, row)) for row in VECTOR.findall(points[m.end():m.end()+end.start()])]
    if len(xyz) != int(m.group(1)) or len(xyz) != (nx+1)*(ny+1)*(nz+1):
        fail("native warped vertex count")
    moved = 0
    for x,y,z in xyz:
        ix,iy,iz = round(nx*x),round(ny*y),round(nz*z)
        if not (0<=ix<=nx and 0<=iy<=ny and 0<=iz<=nz):
            fail("native vertex outside grid")
        bx,by,bz = ix/nx,iy/ny,iz/nz
        boundary = ix in (0,nx) or iy in (0,ny) or iz in (0,nz)
        shift = 0 if boundary else 0.04*math.sin(math.pi*bx)*math.sin(math.pi*by)*math.sin(math.pi*bz)
        if abs(x-bx-shift)>1e-9 or abs(y-by)>1e-9 or abs(z-bz)>1e-9:
            fail("wrong real 3D vertex warp")
        moved += abs(shift)>1e-12
    if moved != (nx-1)*(ny-1)*(nz-1):
        fail("warped vertices not all present")
    boundary = (case/"constant/polyMesh/boundary").read_text()
    for patch in ("xMin","xMax"):
        if not re.search(rf"\b{patch}\s*\{{\s*type\s+cyclic;", boundary):
            fail("native x boundary not cyclic")
    props = (case/"constant/transportProperties").read_text()
    for required in ("activityMode htHto;", "auditActivityFlux true;",
                     "kOx [0 0 -1 0 0 0 0] 0.4;",
                     "lambdaT [0 0 -1 0 0 0 0] 0.2;",
                     "DHT [0 2 -1 0 0 0 0] 0.01;",
                     "DHTO [0 2 -1 0 0 0 0] 0.02;",
                     "momentumMode prescribed;"):
        if required not in props:
            fail(f"wrong physical coefficient {required}")
    schemes = (case/"system/fvSchemes").read_text()
    for required in ("div(phi,HT) Gauss upwind;", "div(phi,HTO) Gauss upwind;",
                     "gradSchemes { default leastSquares; }",
                     "fluxRequired { default no; p; c; HT; HTO; }"):
        if required not in schemes:
            fail(f"wrong scheme {required}")
    control = (case/"system/controlDict").read_text()
    for required in ("deltaT 0.01;", "endTime 0.4;", "writeInterval 40;"):
        if required not in control:
            fail(f"wrong integration {required}")
    centres_data = (case/"0/C").read_text()
    centres_match = VECTORS.search(centres_data)
    if not centres_match or int(centres_match.group(1)) != n:
        fail("native cell centres")
    centres = [tuple(map(float,row)) for row in VECTOR.findall(centres_match.group(2))]
    if len(centres) != n:
        fail("native cell centres length")
    volumes = scalar_field(case/"0/V", n, "0 3 0 0 0 0 0")
    if min(volumes) <= 0 or abs(sum(volumes)-1)>1e-8:
        fail("native cell volumes")
    log = logpath.read_text()
    if len(re.findall(r"(?m)^Time = ",log)) != 40 or len(re.findall(r"Solving for HT,",log)) != 40 or len(re.findall(r"Solving for HTO,",log)) != 40 or "Solving for p," in log:
        fail("wrong number/type of real FVM solves")
    if "mode=prescribed" not in log or "Time = 0.4\n" not in log:
        fail("wrong prescribed-flow mode or final time")
    inventories = [tuple(map(float,v)) for v in INVENTORY.findall(log)]
    hflux = [tuple(map(float,v)) for v in HT_FLUX.findall(log)]
    zflux = [tuple(map(float,v)) for v in HTO_FLUX.findall(log)]
    if len(inventories)!=40 or len(hflux)!=40 or len(zflux)!=40:
        fail("incomplete 40-step activity inventory/flux history")
    start = []
    final = []
    errors = []
    rms = []
    t, u, k, lam, dht, dz = 0.4, 0.25, 0.4, 0.2, 0.01, 0.02
    q = 2*math.pi
    a, b = dht*q*q+k+lam, dz*q*q+lam
    h_mean = 1e6*math.exp(-(k+lam)*t)
    z_mean = math.exp(-lam*t)*(2e5+1e6*(1-math.exp(-k*t)))
    h_amp = 2e5*math.exp(-a*t)
    z_amp = 5e4*math.exp(-b*t) + k*2e5*(math.exp(-a*t)-math.exp(-b*t))/(b-a)
    for name, mean, amp, initial_mean, initial_amp, maxerr in (
        ("HT",h_mean,h_amp,1e6,2e5,limits[0]),
        ("HTO",z_mean,z_amp,2e5,5e4,limits[1]),
    ):
        init = scalar_field(case/"0"/name,n,"0 -3 -1 0 0 0 0")
        out = scalar_field(case/"0.4"/name,n,"0 -3 -1 0 0 0 0")
        for patch in ("xMin","xMax"):
            if not re.search(rf"\b{patch}\s*\{{\s*type\s+cyclic;",(case/"0"/name).read_text()):
                fail(f"{name} cyclic boundary missing")
        expected0 = [initial_mean+initial_amp*math.cos(q*x) for x,y,z in centres]
        expected = [mean+amp*math.cos(q*(x-u*t)) for x,y,z in centres]
        if max(abs(v-e) for v,e in zip(init,expected0)) > 1e-6:
            fail(f"{name} wrong native-centre initialization")
        if min(out)<0:
            fail(f"{name} negative activity concentration")
        err = [abs(v-e) for v,e in zip(out,expected)]
        mx = max(err)
        r = math.sqrt(sum(e*e for e in err)/n)
        if mx > maxerr:
            fail(f"{name} continuum error Linf={mx:.12g} > {maxerr}")
        if abs(sum(v*vol for v,vol in zip(out,volumes)) - inventories[-1][0 if name=="HT" else 1]) > 1e-5:
            fail(f"{name} actual-volume integral disagrees with native log")
        errors.append(mx)
        rms.append(r)
        start.append(sum(v*vol for v,vol in zip(init,volumes)))
        final.append(sum(v*vol for v,vol in zip(out,volumes)))
    old_h, old_z = start
    max_budget, min_cyclic = 0.,float("inf")
    for i,((h,z),hf,zf) in enumerate(zip(inventories,hflux,zflux),1):
        hl,hr,hn,hreported = hf
        zl,zr,zn,zreported,totalreported = zf
        if not all(math.isfinite(v) for v in (h,z,*hf,*zf)) or min(h,z)<0:
            fail(f"nonfinite/negative inventory or flux step {i}")
        if max(abs(hl+hr),abs(hn),abs(zl+zr),abs(zn))>1e-5:
            fail(f"cyclic signed activity flux not balanced step {i}")
        min_cyclic = min(min_cyclic,abs(hl),abs(hr),abs(zl),abs(zr))
        if min(abs(hl),abs(hr),abs(zl),abs(zr))<=1e3:
            fail(f"no genuine cyclic species transport step {i}")
        rh = h-old_h+0.01*hn+0.01*(k+lam)*h
        rz = z-old_z+0.01*zn+0.01*lam*z-0.01*k*h
        rt = h+z-old_h-old_z+0.01*(hn+zn)+0.01*lam*(h+z)
        max_budget=max(max_budget,*map(abs,(rh,rz,rt,hreported,zreported,totalreported,
                   rh-hreported,rz-zreported,rt-totalreported)))
        if max(map(abs,(rh,rz,rt,hreported,zreported,totalreported)))>1e-5:
            fail(f"independent species/total activity budget step {i}")
        old_h,old_z=h,z
    print(f"M1C-1 cells={n} warped_vertices={moved} nonorth_deg={float(nonorth.group(1)):.8g} "
          f"HT_Linf={errors[0]:.12g} HT_RMS={rms[0]:.12g} "
          f"HTO_Linf={errors[1]:.12g} HTO_RMS={rms[1]:.12g} "
          f"min_cyclic_flux_Bq_s={min_cyclic:.12g} max_budget_Bq={max_budget:.12g} "
          f"final_HT_Bq={final[0]:.12g} final_HTO_Bq={final[1]:.12g}")
    return rms

if len(sys.argv)!=7:
    raise SystemExit("usage: check_spatial_activity.py COARSE_CASE COARSE_MESH COARSE_SOLVE FINE_CASE FINE_MESH FINE_SOLVE")
coarse = audit(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),(8,20,6),(45000,18000))
fine = audit(Path(sys.argv[4]),Path(sys.argv[5]),Path(sys.argv[6]),(12,30,9),(33000,14000))
for name,c,f in zip(("HT","HTO"),coarse,fine):
    ratio = f/c if c>0 else float("inf")
    print(f"M1C-1 {name} refinement_RMS_ratio={ratio:.12g}")
    if not (0<f and 0<c and ratio<=0.85):
        fail(f"{name} spatial refinement criterion")
print("PASS: M1C-1 official OpenFOAM v2606 spatial HT/HTO activity verification")
