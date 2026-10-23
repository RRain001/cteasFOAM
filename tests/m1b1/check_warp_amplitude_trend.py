#!/usr/bin/env python3
"""Independent amplitude-response check on already-verified warped ABC output.

check_warped_abc.py does full official-geometry/centroid/analytic U,p,
conservation and grid refinement for each amplitude. This checker tests
the systematic geometry perturbation and checks for material accuracy loss.
"""
from pathlib import Path
import math
import re
import sys

pattern=re.compile(
    r"M1B-1f nonorth_ABC n=(8|16) "
    r"max_nonorth_deg=([-+0-9.eE]+) "
    r"max_skewness=([-+0-9.eE]+).*?"
    r"max_U_error=([-+0-9.eE]+) "
    r"gauge_free_p_error=([-+0-9.eE]+)"
)

def read_log(path):
    values={}
    for line in Path(path).read_text().splitlines():
        match=pattern.search(line)
        if match:
            n=int(match.group(1))
            values[n]=tuple(map(float,match.groups()[1:]))
    if set(values)!={8,16}:
        raise SystemExit(f"FAIL: missing real mesh/ABC fields in {path}")
    if any(not math.isfinite(x) for row in values.values() for x in row):
        raise SystemExit("FAIL: nonfinite physical sensitivity metric")
    return values

if len(sys.argv)!=3:
    raise SystemExit("usage: check_warp_amplitude_trend.py METRICS_A04 METRICS_A08")

small=read_log(sys.argv[1])
large=read_log(sys.argv[2])
for n in (8,16):
    nonorth_s,skew_s,u_s,p_s=small[n]
    nonorth_l,skew_l,u_l,p_l=large[n]
    ur=u_l/u_s
    pr=p_l/p_s
    print(f"M1B-1f4 warpSensitivity n={n} "
          f"nonorth04={nonorth_s:.12g} nonorth08={nonorth_l:.12g} "
          f"skew04={skew_s:.12g} skew08={skew_l:.12g} "
          f"U_error_ratio_08_to_04={ur:.12g} p_error_ratio_08_to_04={pr:.12g}")
    if not (nonorth_l>nonorth_s+2 and nonorth_l/nonorth_s>1.5):
        raise SystemExit("FAIL: physical maximum nonorthogonality fails amplitude-response test")
    if not skew_l>skew_s:
        raise SystemExit("FAIL: actual OpenFOAM maximum skewness does not respond to amplitude")
    if not (0.04<ur<=1.5 and 0.04<pr<=1.6):
        raise SystemExit("FAIL: numerical solution degrades or suspiciously collapses")
print("PASS: controlled mild/moderate nonorthogonality trend and bounded ABC error")
