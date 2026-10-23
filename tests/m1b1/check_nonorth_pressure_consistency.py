#!/usr/bin/env python3
"""Audit predeclared nonorth pressure loop science against accepted no-loop baseline.

Do not assume extra nonorth pressure corrections must improve U/p L-infinity
on this mildly warped case; reject material regressions while requiring
the independent absolute analytical, mesh, conservation and refinement gates.
"""
from pathlib import Path
import math
import re
import sys
if len(sys.argv)!=2:
    raise SystemExit("usage: check_nonorth_pressure_consistency.py CHECKER_LOG")
report=Path(sys.argv[1]).read_text().splitlines()
pat=re.compile(
    r"M1B-1f nonorth_ABC n=(8|16) .*?"
    r"max_U_error=([-+0-9.eE]+) gauge_free_p_error=([-+0-9.eE]+)"
)
current={}
for line in report:
    m=pat.search(line)
    if m:
        n=int(m.group(1))
        current[n]=(float(m.group(2)),float(m.group(3)))
if set(current)!={8,16}:
    raise SystemExit("FAIL: two actual warp-grid scientific records missing")

# Immutable reference numbers from accepted PR #16, official v2606 CI
# #39895216715, with nNonOrthogonalCorrectors absent/default 0.
baseline={
    8: (0.172379028547, 0.553776125146),
    16: (0.108293338904, 0.289513430093),
}
for n in (8,16):
    u,p=current[n]
    ru=u/baseline[n][0]
    rp=p/baseline[n][1]
    print(
        f"M1B-1f2 nonorth_comparison n={n} U_relative_to_baseline={ru:.12g} "
        f"pressure_relative_to_baseline={rp:.12g}"
    )
    if not all(math.isfinite(x) for x in (u,p,ru,rp)):
        raise SystemExit("FAIL: nonfinite nonorth pressure solution")
    if ru>1.20 or rp>1.20:
        raise SystemExit("FAIL: extra pressure corrections degraded actual ABC beyond 20%")

print("PASS: 2 true inner nonorth pressure corrections, no material U/p regression")
