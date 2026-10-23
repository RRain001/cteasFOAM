#!/usr/bin/env bash
# M1B-1g-4: real native v2606 icoFoam vs cteasAtmosphereFoam.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
rc=$?
set -e
echo "M1B-1g-4 official OpenFOAM setup_rc=$rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for app in python3 icoFoam blockMesh checkMesh postProcess wmake wclean timeout; do
 command -v "$app" >/dev/null || { echo "FAIL: required official executable absent: $app"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-icofoam-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1B-1g-4 upstream comparison rc=$rc";
 for f in "$work"/build.log "$work"/mesh.log "$work"/check.log "$work"/centres.log "$work"/*solve*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -70 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || { cat "$work/build.log"; exit 21; }
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then
 echo "FAIL: hidden compiler failure"; exit 22;
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for mode in warped zero; do
 case_dir="$work/cteas-$mode"
 ref_dir="$work/ico-$mode"
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" "$mode"
 if [[ "$mode" == warped ]]; then
  blockMesh -case "$case_dir" > "$work/mesh.log" 2>&1 || { tail -80 "$work/mesh.log"; exit 23; }
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" "$mode"
  checkMesh -case "$case_dir" -allGeometry -allTopology > "$work/check.log" 2>&1 || { tail -100 "$work/check.log"; exit 24; }
  grep -q 'Mesh OK.' "$work/check.log" || exit 24
  postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$work/centres.log" 2>&1 || {
   tail -90 "$work/centres.log"; exit 25;
  }
 else
  cp -a "$work/cteas-warped/constant/polyMesh" "$case_dir/constant/"
  cp -a "$work/cteas-warped/0/C" "$case_dir/0/"
 fi
 python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$case_dir" "$mode"
 cat >> "$case_dir/system/fvSolution" <<'ENDPISO'
PISO
{
    momentumPredictor yes;
    nCorrectors 3;
    nNonOrthogonalCorrectors 2;
    pRefCell 0;
    pRefValue 0;
}
ENDPISO
 # Official v2606 icoFoam selects pFinal on its final PISO iteration.
 # Add the same solver setting to both candidates, preserving equality.
 python3 - "$case_dir/system/fvSolution" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); data=p.read_text()
start = data.index("solvers")
opening = data.index("{", start)
depth, end = 0, None
for idx in range(opening,len(data)):
    if data[idx]=="{": depth += 1
    elif data[idx]=="}":
        depth -= 1
        if depth==0: end=idx; break
if end is None or "pFinal" in data:
    raise SystemExit("FAIL: invalid pressure solver dictionary")
p.write_text(data[:end] + 'pFinal { $p; relTol 0; }\n' + data[end:])
PY
 cp -a "$case_dir" "$ref_dir"
 python3 - "$ref_dir/system/controlDict" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]);data=p.read_text()
old="application cteasAtmosphereFoam;"
if data.count(old)!=1: raise SystemExit("FAIL: mismatching application field")
p.write_text(data.replace(old, "application icoFoam;"))
PY
 for rel in constant/polyMesh/points constant/polyMesh/faces constant/polyMesh/owner constant/polyMesh/neighbour constant/polyMesh/boundary 0/C 0/U 0/p 0/c constant/transportProperties system/fvSchemes system/fvSolution; do
  cmp -s "$case_dir/$rel" "$ref_dir/$rel" || {
   echo "FAIL: independent case initial inputs differ $rel"; exit 26;
  }
 done
 echo "M1B-1g-4 exact identical native initial mesh/U/p/c for mode=$mode"
 timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$work/cteas-$mode-solve.log" 2>&1 || {
  tail -100 "$work/cteas-$mode-solve.log"; exit 31;
 }
 echo "M1B-1g-4 completed CTEAS mode=$mode"
 timeout 180s icoFoam -case "$ref_dir" > "$work/ico-$mode-solve.log" 2>&1 || {
  tail -100 "$work/ico-$mode-solve.log"; exit 32;
 }
 echo "M1B-1g-4 completed official icoFoam mode=$mode"
done
python3 /workspace/tests/m1b1/check_icofoam_cross_solver.py \
 "$work/cteas-warped" "$work/ico-warped" "$work/check.log" \
 "$work/cteas-warped-solve.log" "$work/ico-warped-solve.log" \
 "$work/cteas-zero" "$work/ico-zero" \
 "$work/cteas-zero-solve.log" "$work/ico-zero-solve.log"
# Deliberately contaminate a reference final field and ensure the
# evaluator flags the discrepancy, ruling out a silent same-file
# comparison. The corrupted output is confined to the disposable case.
python3 - "$work/ico-warped/0.4/U" <<'PY'
from pathlib import Path
import re, sys
p=Path(sys.argv[1])
txt=p.read_text()
pat=(r"(internalField\s+nonuniform\s+List<vector>\s+\d+\s*\(\s*)"
     r"\(([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\)")
m=re.search(pat,txt)
if m is None:
    raise SystemExit("FAIL: unable to locate actual upstream final U vector for negative check")
changed="(" + f"{float(m.group(2))+0.02:.16g}" + " " + m.group(3) + " " + m.group(4) + ")"
p.write_text(txt[:m.start(2)-1] + changed + txt[m.end():])
PY
if python3 /workspace/tests/m1b1/check_icofoam_cross_solver.py \
 "$work/cteas-warped" "$work/ico-warped" "$work/check.log" \
 "$work/cteas-warped-solve.log" "$work/ico-warped-solve.log" \
 "$work/cteas-zero" "$work/ico-zero" \
 "$work/cteas-zero-solve.log" "$work/ico-zero-solve.log" \
 > "$work/injected-check.log" 2>&1; then
 echo "FAIL: corrupted official final velocity went undetected"
 exit 35
fi
grep -q 'FAIL: official icoFoam and experimental solver diverge beyond predeclared bounds' "$work/injected-check.log" || {
 cat "$work/injected-check.log"; echo "FAIL: wrong negative-check rejection"; exit 36;
}
echo "PASS: injected official upstream final U corruption rejected by independent comparison"
echo "PASS: M1B-1g-4 real official v2606 upstream solver cross-reference"
