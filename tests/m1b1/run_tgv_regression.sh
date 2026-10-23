#!/usr/bin/env bash
# Official OpenFOAM-v2606 true PISO-style pressure-corrected TGV benchmark.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B1 setup_rc=$setup_rc WM_PROJECT_VERSION=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
set -u
work=$(mktemp -d /tmp/cteasfoam-tgv-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
  echo "M1B1 TGV FAILED, rc=$rc";
  for f in "$work"/compile.log "$work"/mesh*.log "$work"/solve*.log; do
    if [[ -f "$f" ]]; then echo "=== $f ==="; tail -140 "$f"; fi
  done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(
 cd "$work/atmosphere"
 wclean
 wmake > "$work/compile.log" 2>&1 || { cat "$work/compile.log"; exit 21; }
)
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/compile.log"; then
  echo "FAIL: wmake reported hidden compiler error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for n in 12 24; do
  python3 /workspace/tests/m1b1/make_tgv_case.py "$work/case$n" "$n"
  blockMesh -case "$work/case$n" > "$work/mesh$n.log" 2>&1 || {
      cat "$work/mesh$n.log"; exit 23;
  }
  echo "[M1B1] actual pressure-corrected TGV n=$n"
  timeout 120s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" \
     -case "$work/case$n" > "$work/solve$n.log" 2>&1 || {
      cat "$work/solve$n.log"; exit 24;
  }
  grep -q 'CTEAS_M1B1 max_div_phi=' "$work/solve$n.log" || {
    echo "FAIL: no actual pressure continuity records"; exit 25;
  }
done
python3 /workspace/tests/m1b1/check_tgv.py \
 "$work/case12" "$work/solve12.log" "$work/case24" "$work/solve24.log"
echo "PASS: M1B1 pressure-corrected TGV reference"
