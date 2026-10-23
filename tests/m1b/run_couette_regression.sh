#!/usr/bin/env bash
# M1B-0 deliberately restricted momentum PDE: unsteady planar Couette.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B setup_rc=$setup_rc WM_PROJECT_VERSION=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for command_name in wmake wclean blockMesh python3 timeout; do
    command -v "$command_name" >/dev/null || { echo "FAIL: missing $command_name"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-couette-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
    for f in "$work"/compile.log "$work"/mesh.log "$work"/solver.log; do
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
    wmake > "$work/compile.log" 2>&1 || { cat "$work/compile.log"; exit 21; }
)
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/compile.log"; then
    echo "FAIL: hidden C++ compiler diagnostics"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
python3 /workspace/tests/m1b/make_couette_case.py "$work/case"
blockMesh -case "$work/case" > "$work/mesh.log" 2>&1 || {
    cat "$work/mesh.log"; exit 23;
}
echo "[M1B-0] solving transient viscous Couette momentum PDE"
timeout 120s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$work/case" \
    > "$work/solver.log" 2>&1 || {
    cat "$work/solver.log"; exit 24;
}
echo "[M1B-0] checking exact transient Couette eigenfunction"
python3 /workspace/tests/m1b/check_couette.py "$work/case" "$work/solver.log"
echo "PASS: M1B-0 restricted laminar momentum PDE verified"
