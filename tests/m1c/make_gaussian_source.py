#!/usr/bin/env python3
"""Prepare real OpenFOAM M1C-6 Gaussian injection cases from M1C-5."""
from pathlib import Path
import subprocess
import sys

VALID = ("coarse", "fine", "unknown_shape", "missing_center", "negative_sigma",
         "bad_center_dims", "bad_sigma_dims", "orphan_shape", "no_support")
if len(sys.argv) != 3 or sys.argv[2] not in VALID:
    raise SystemExit("usage: make_gaussian_source.py CASE " + "|".join(VALID))
case, mode = Path(sys.argv[1]), sys.argv[2]
subprocess.run([sys.executable, str(Path(__file__).with_name("make_timed_emission.py")),
                str(case), "pulse"], check=True)
if mode == "fine":
    path = case / "system/blockMeshDict"
    text = path.read_text()
    before = "(4 3 2) simpleGrading"
    if text.count(before) != 1:
        raise SystemExit("FAIL: source mesh 4x3x2 template missing")
    path.write_text(text.replace(before, "(8 6 4) simpleGrading"))
props = case / "constant/transportProperties"
if mode == "orphan_shape":
    data = props.read_text()
    for line in data.splitlines(keepends=True):
        if line.startswith("sourceHT_Bq_s ") or line.startswith("sourceHTO_Bq_s "):
            data = data.replace(line, "")
    props.write_text(data)
with props.open("a") as f:
    shape = "unsupported" if mode == "unknown_shape" else "gaussian"
    f.write(f"sourceShape {shape};\n")
    if mode != "missing_center":
        dims = "[0 0 0 0 0 0 0]" if mode == "bad_center_dims" else "[0 1 0 0 0 0 0]"
        center = "(100 100 100)" if mode == "no_support" else "(0.38 0.62 0.55)"
        f.write(f"sourceCenter {dims} {center};\n")
    dims = "[0 0 0 0 0 0 0]" if mode == "bad_sigma_dims" else "[0 1 0 0 0 0 0]"
    sigma = "-0.22" if mode == "negative_sigma" else "0.22"
    f.write(f"sourceSigma {dims} {sigma};\n")
print(f"M1C-6 mode={mode} sourceShape={shape} sourceCenter=(0.38 0.62 0.55) sigma={sigma} native volume mesh")
