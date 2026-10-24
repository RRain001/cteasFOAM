#!/usr/bin/env python3
"""Prepare an isolated M1C-5 closed 3-D activity-source acceptance case."""
from pathlib import Path
import subprocess
import sys

if len(sys.argv) != 3 or sys.argv[2] not in (
    "pulse", "ht_only", "zero_rate", "negative_rate", "reverse_window",
    "missing_window", "wrong_dimensions"
):
    raise SystemExit("usage: make_timed_emission.py CASE pulse|ht_only|zero_rate|negative_rate|reverse_window|missing_window|wrong_dimensions")
case, mode = Path(sys.argv[1]), sys.argv[2]
subprocess.run(
    [sys.executable, str(Path(__file__).with_name("make_activity_box.py")),
     str(case), "0.4", "0.2", "valid"],
    check=True
)
for name, initial in (("HT", "1000000"), ("HTO", "200000")):
    path = case / "0" / name
    src = path.read_text()
    before = "internalField uniform " + initial + ";"
    if src.count(before) != 1:
        raise SystemExit("FAIL: unexpected HT/HTO initial field")
    path.write_text(src.replace(before, "internalField uniform 0;"))
props = case / "constant" / "transportProperties"
q_ht, q_hto = 60000, 12000
if mode == "ht_only":
    q_hto = 0
if mode == "zero_rate":
    q_ht = q_hto = 0
if mode == "negative_rate":
    q_ht = -60000
start, end = (0.225, 0.075) if mode == "reverse_window" else (0.075, 0.225)
rate_dim = "[0 0 -1 0 0 0 0]" if mode == "wrong_dimensions" else "[0 0 -2 0 0 0 0]"
with props.open("a") as f:
    f.write(f"sourceHT_Bq_s {rate_dim} {q_ht};\n")
    f.write(f"sourceHTO_Bq_s [0 0 -2 0 0 0 0] {q_hto};\n")
    f.write(f"sourceStart [0 0 1 0 0 0 0] {start:.17g};\n")
    if mode != "missing_window":
        f.write(f"sourceEnd [0 0 1 0 0 0 0] {end:.17g};\n")
print(f"M1C-5 mode={mode} HT_rate={q_ht} HTO_rate={q_hto} window=[{start},{end})")
