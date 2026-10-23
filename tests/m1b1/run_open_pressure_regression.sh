#!/usr/bin/env bash
# Actual static-pressure inlet and outlet (NOT periodic equivalent body force).
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B-1e fixed pressure setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
set -u
work=$(mktemp -d /tmp/cteasfoam-fixedpressure-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "M1B-1e failed rc=$rc";
 for f in "$work"/build.log "$work"/mesh*.log "$work"/solve*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -150 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(
 cd "$work/atmosphere"
 wclean
 wmake > "$work/build.log" 2>&1 || { cat "$work/build.log"; exit 21; }
)
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then
 echo "FAIL: hidden solver compiler error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for kind in driven zero; do
 python3 /workspace/tests/m1b1/make_open_pressure_case.py "$work/$kind" "$kind"
 blockMesh -case "$work/$kind" > "$work/mesh$kind.log" 2>&1 || {
  cat "$work/mesh$kind.log"; exit 23;
 }
 echo "[M1B-1e] actual p-fixedValue inlet/outlet, case=$kind"
 timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" \
   -case "$work/$kind" > "$work/solve$kind.log" 2>&1 || {
   cat "$work/solve$kind.log"; exit 24;
 }
done
python3 /workspace/tests/m1b1/check_open_pressure.py \
 "$work/driven" "$work/solvedriven.log" "$work/zero" "$work/solvezero.log"
echo "PASS: actual two-end static-pressure-inlet/outlet Poiseuille and zero-dp control"
