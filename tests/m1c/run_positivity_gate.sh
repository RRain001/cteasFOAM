#!/usr/bin/env bash
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c17-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then echo "FAIL: M1C17 exit=$rc"; for f in "$work"/build.log "$work"/*/solve.log; do if [[ -f "$f" ]]; then echo "== $f =="; tail -55 "$f"; fi; done; fi; rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere" && wclean && wmake > "$work/build.log" 2>&1) || { tail -95 "$work/build.log"; exit 21; }
exe="$FOAM_USER_APPBIN/cteasAtmosphereFoam"
test -x "$exe" || exit 22
for kind in base guarded unsafe badseed; do
    case_dir="$work/$kind"
    python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" zero 8x8x8
    python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$case_dir" advect
    python3 /workspace/tests/m1c/make_gaussian_green.py "$case_dir" normal
    python3 - "$case_dir" "$kind" <<'PYCASE'
from pathlib import Path
import sys
p = Path(sys.argv[1]); kind = sys.argv[2]
props = p/'constant/transportProperties'
props.write_text(props.read_text()+'auditActivityExtrema true;\n'+('enforceNonnegativeActivity true;\n' if kind != 'base' else ''))
if kind == 'unsafe':
    f = p/'system/fvSchemes'; txt = f.read_text()
    for sp in ('HT','HTO'):
        old = f'div(phi,{sp}) Gauss upwind;'
        assert txt.count(old) == 1
        txt = txt.replace(old, f'div(phi,{sp}) Gauss limitedLinear 1;')
    f.write_text(txt)
if kind == 'badseed':
    f = p/'0/HT'; txt = f.read_text(); old = 'internalField uniform 0;'
    assert txt.count(old) == 1
    f.write_text(txt.replace(old, 'internalField uniform -1;'))
print(f'M1C17_NATIVE_INPUT kind={kind} opt_in={kind != "base"}')
PYCASE
    blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || exit 23
    python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" zero 8x8x8
    checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || exit 24
    grep -q 'Mesh OK.' "$case_dir/mesh-check.log" || exit 25
    postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$case_dir/centres.log" 2>&1 || exit 26
    postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$case_dir/volumes.log" 2>&1 || exit 27
    if [[ "$kind" == base || "$kind" == guarded ]]; then
        timeout 480s "$exe" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || { echo "FAIL: expected positive native completion $kind"; exit 31; }
    else
        if timeout 480s "$exe" -case "$case_dir" > "$case_dir/solve.log" 2>&1; then
            echo "FAIL: nonnegative guard did not abort unsafe $kind"; exit 32
        fi
        python3 /workspace/tests/m1c/check_positivity_gate.py --reject "$case_dir" "$([[ "$kind" == unsafe ]] && echo postSolve || echo initial)"
    fi
done
python3 /workspace/tests/m1c/check_positivity_gate.py --compare "$work/base" "$work/guarded"
# Negative parser control: real native unsafe log with marker removed must be rejected.
cp "$work/unsafe/solve.log" "$work/unsafe/native-original.log"
sed -i 's/CTEAS_M1C_POSITIVITY_FAILURE/M1C17_TAMPERED/' "$work/unsafe/solve.log"
if python3 /workspace/tests/m1c/check_positivity_gate.py --reject "$work/unsafe" postSolve > "$work/tamper.log" 2>&1; then
    echo 'FAIL: forged success accepted after genuine failure marker removed'; exit 40
fi
grep -q 'expected exactly one genuine native positivity failure' "$work/tamper.log" || { cat "$work/tamper.log"; exit 41; }
echo 'PASS: M1C-17 real v2606 two upwind field-preserving solves, real limitedLinear intermediate fail, initial-field fail, tamper negative control'
