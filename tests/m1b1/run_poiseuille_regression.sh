#!/usr/bin/env bash
# Verify finite-volume wall-driven-equivalent Poiseuille in official v2606.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B1 Poiseuille setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
set -u
work=$(mktemp -d /tmp/cteasfoam-poiseuille-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "M1B1 Poiseuille regression failed rc=$rc";
 for f in "$work"/build.log "$work"/mesh.log "$work"/solver.log; do
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
 echo "FAIL: hidden build error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
python3 /workspace/tests/m1b1/make_poiseuille_case.py "$work/case"
blockMesh -case "$work/case" > "$work/mesh.log" 2>&1 || {
 cat "$work/mesh.log"; exit 23;
}
echo "[M1B-1c] Solve 100-step transient no-slip channel with equivalent dp/dx"
timeout 120s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$work/case" \
  > "$work/solver.log" 2>&1 || { cat "$work/solver.log"; exit 24; }
python3 /workspace/tests/m1b1/check_poiseuille.py "$work/case" "$work/solver.log"
echo "PASS: M1B-1c 3-D mesh with nonperiodic no-slip wall reference"
