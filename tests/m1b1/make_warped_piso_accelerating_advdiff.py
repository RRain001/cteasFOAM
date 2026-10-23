#!/usr/bin/env python3
"""Prepare real warped 3-D cyclic-x PISO acceleration and passive-scalar case."""
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE / "make_warped_periodic_advection_diffusion.py"


def replace_once(path, before, after):
    source = path.read_text()
    if source.count(before) != 1:
        raise SystemExit(f"FAIL: expected unique {before!r} in {path}")
    path.write_text(source.replace(before, after))


def prepare(case):
    subprocess.run([sys.executable, str(BASE), "prepare", str(case)], check=True)
    props = case / "constant/transportProperties"
    replace_once(props, "momentumMode prescribed;", """momentumMode pisoLaminar;
nu [0 2 -1 0 0 0 0] 0.05;
drivingAcceleration [0 1 -2 0 0 0 0] (0.5 0 0);
nPressureCorrectors 3;
nNonOrthogonalCorrectors 2;""")
    u = case / "0/U"
    u.write_text("""FoamFile { version 2.0; format ascii; class volVectorField; object U; }
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
 xMin { type cyclic; }
 xMax { type cyclic; }
 yMin { type slip; }
 yMax { type slip; }
 zMin { type slip; }
 zMax { type slip; }
}
""")
    (case / "0/p").write_text("""FoamFile { version 2.0; format ascii; class volScalarField; object p; }
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
 xMin { type cyclic; }
 xMax { type cyclic; }
 yMin { type fixedFluxPressure; value uniform 0; }
 yMax { type fixedFluxPressure; value uniform 0; }
 zMin { type fixedFluxPressure; value uniform 0; }
 zMax { type fixedFluxPressure; value uniform 0; }
}
""")
    print("M1B-1g-11 configured genuine PISO U0=0 a=(0.5,0,0) D=0.01 dt=0.01")


def seed(case, count):
    subprocess.run([sys.executable, str(BASE), "seed", str(case), str(count)], check=True)


if len(sys.argv) not in (3, 4) or sys.argv[1] not in ("prepare", "seed"):
    raise SystemExit("usage: make_warped_piso_accelerating_advdiff.py prepare CASE | seed CASE COUNT")
case = Path(sys.argv[2])
if sys.argv[1] == "prepare" and len(sys.argv) == 3:
    prepare(case)
elif sys.argv[1] == "seed" and len(sys.argv) == 4:
    seed(case, int(sys.argv[3]))
else:
    raise SystemExit("FAIL: wrong generator arguments")
