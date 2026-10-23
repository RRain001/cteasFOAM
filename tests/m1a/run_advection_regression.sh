#!/usr/bin/env bash
# Independent real OpenFOAM nonzero-velocity scalar transport benchmark.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1A setup_rc=$setup_rc WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
[[ "${WM_PROJECT_VERSION:-}" == v2606 || "${WM_PROJECT_VERSION:-}" == 2606 ]] || exit 18
set -u
work=$(mktemp -d /tmp/cteasfoam-advect-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
    for f in "$work"/compile.log "$work"/mesh*.log "$work"/solve*.log; do
        if [[ -f "$f" ]]; then echo "=== $f ==="; tail -100 "$f"; fi
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
    echo "FAIL: hidden wmake error"; exit 22
fi
for n in 40 80; do
    python3 /workspace/tests/m1a/make_advection_case.py "$work/case$n" "$n"
    (
        cd "$work/case$n"
        blockMesh > "$work/mesh$n.log" 2>&1 || { cat "$work/mesh$n.log"; exit 23; }
        "$FOAM_USER_APPBIN/cteasAtmosphereFoam" > "$work/solve$n.log" 2>&1 || {
            cat "$work/solve$n.log"; exit 24;
        }
    )
    grep -q 'CTEAS_M1A integral=' "$work/solve$n.log" || {
        echo "FAIL: solver wrote no integral audit"; exit 25;
    }
done
python3 /workspace/tests/m1a/check_advection.py "$work/case40" "$work/case80"
echo "PASS: M1A nonzero periodic advection"
