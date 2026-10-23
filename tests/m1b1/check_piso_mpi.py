#!/usr/bin/env python3
"""Compare real two-rank pressure-corrected Taylor-Green with serial solve."""
from pathlib import Path
import math
import re
import sys

if len(sys.argv)!=5:
    raise SystemExit("usage: check_piso_mpi.py SERIAL_CASE MPI_CASE SERIAL_LOG MPI_LOG")
vecpat=re.compile(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
scapat=re.compile(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",re.S)
logpat=re.compile(r"CTEAS_M1B1 max_div_phi=([-+0-9.eE]+) boundary_net_flux=([-+0-9.eE]+) kinetic_energy=([-+0-9.eE]+)")
def field(folder,name,kind):
    f=Path(folder)/"0.1"/name
    if not f.is_file():
        raise SystemExit(f"FAIL: missing real reconstructed output {f}")
    text=f.read_text()
    m=(vecpat if kind=="vector" else scapat).search(text)
    if not m or int(m.group(1))!=144:
        raise SystemExit(f"FAIL: malformed {name} field at {f}")
    if kind=="vector":
        x=[tuple(float(w) for w in z.split()) for z in re.findall(r"\(([^()]*)\)",m.group(2))]
        if len(x)!=144 or any(len(v)!=3 for v in x):
            raise SystemExit("FAIL: malformed vector field")
        flat=[z for v in x for z in v]
    else:
        x=[float(z) for z in m.group(2).split()]
        flat=x
    if len(x)!=144 or not all(math.isfinite(z) for z in flat):
        raise SystemExit("FAIL: nonfinite field or wrong count")
    return x
def log_check(path):
    txt=Path(path).read_text()
    data=[tuple(map(float,z)) for z in logpat.findall(txt)]
    if len(data)!=40:
        raise SystemExit(f"FAIL: pressure PDE did not complete 40 timesteps in {path}")
    if not all(math.isfinite(z) for row in data for z in row):
        raise SystemExit("FAIL: nonfinite MPI flux or energy")
    return data

serial,parallel,slog,plog=sys.argv[1:]
sU,pU=field(serial,"U","vector"),field(parallel,"U","vector")
sp,pp=field(serial,"p","scalar"),field(parallel,"p","scalar")
sd,pd=log_check(slog),log_check(plog)
udiff=max(abs(x-y) for v,w in zip(sU,pU) for x,y in zip(v,w))
sgauge=sum(sp)/len(sp)
pgauge=sum(pp)/len(pp)
pdiff=max(abs((x-sgauge)-(y-pgauge)) for x,y in zip(sp,pp))
div=max(abs(z[0]) for z in pd)
net=max(abs(z[1]) for z in pd)
energy_diff=abs(sd[-1][2]-pd[-1][2])
print(
    f"M1B-1b MPI max_U_difference={udiff:.12g} "
    f"gauge_free_p_difference={pdiff:.12g} "
    f"max_div_phi={div:.12g} max_net_boundary_flux={net:.12g} "
    f"final_kinetic_energy_difference={energy_diff:.12g}"
)
if div>1e-7 or net>1e-7:
    raise SystemExit("FAIL: parallel corrected flux continuity")
if udiff>1e-3 or pdiff>2e-3 or energy_diff>1e-4:
    raise SystemExit("FAIL: 2-rank pressure-corrected PDE differs from serial")
print("PASS: pressure-corrected MPI velocity, gauge-free pressure, and continuity")
