#!/usr/bin/env bash
# M1B-1g-2: time-step self-convergence, identical physical 3D warped duct.
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
setup_rc=$?
set -e
echo "M1B-1g-2 official OpenFOAM setup_rc=$setup_rc version=$WM_PROJECT_VERSION"
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
for cmd in python3 wmake wclean blockMesh checkMesh postProcess timeout; do
 command -v "$cmd" >/dev/null || { echo "FAIL: missing $cmd"; exit 18; }
done
set -u
work=$(mktemp -d /tmp/cteasfoam-duct-dt-XXXXXX)
trap 'rc=$?; if ((rc != 0)); then
 echo "FAIL: M1B-1g-2 time-step study rc=$rc";
 for f in "$work"/build.log "$work"/mesh.log "$work"/check.log "$work"/centres.log "$work"/solve*.log; do
  if [[ -f "$f" ]]; then echo "=== $f ==="; tail -70 "$f"; fi
 done
fi
rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere"; wclean; wmake > "$work/build.log" 2>&1) || {
 cat "$work/build.log"; exit 21;
}
if grep -E '(^|[[:space:]])(fatal )?error:|wmkdepend: could not open' "$work/build.log"; then
 echo "FAIL: hidden compiler error"; exit 22
fi
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
base="$work/dt20"
python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$base" warped
blockMesh -case "$base" > "$work/mesh.log" 2>&1 || {
 tail -90 "$work/mesh.log"; exit 23;
}
python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$base" warped || exit 24
checkMesh -case "$base" -allGeometry -allTopology > "$work/check.log" 2>&1 || {
 tail -120 "$work/check.log"; exit 25;
}
grep -q 'Mesh OK.' "$work/check.log" || exit 25
postProcess -case "$base" -func writeCellCentres -time 0 > "$work/centres.log" 2>&1 || {
 tail -100 "$work/centres.log"; exit 26;
}
test -s "$base/0/C" || exit 26
python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py seed "$base" warped || exit 27
for tag in 10 05; do
 cp -a "$base" "$work/dt$tag"
done
python3 - "$work" <<'PY'
from pathlib import Path
import sys
folder = Path(sys.argv[1])
for tag, dt in (("10","0.01"), ("05","0.005")):
    path = folder / ("dt"+tag) / "system" / "controlDict"
    data = path.read_text()
    if data.count("deltaT 0.02;") != 1:
        raise SystemExit("FAIL: expected exactly one baseline deltaT")
    path.write_text(data.replace("deltaT 0.02;", "deltaT "+dt+";"))
PY
for tag in 10 05; do
 for relative in constant/polyMesh/points constant/polyMesh/faces constant/polyMesh/owner constant/polyMesh/neighbour constant/polyMesh/boundary 0/C 0/U 0/p 0/c; do
  cmp -s "$base/$relative" "$work/dt$tag/$relative" || {
   echo "FAIL: changed input $tag/$relative"; exit 28;
  }
 done
done
echo 'PASS: bitwise identical actual mesh, centroids, U/p/c and pressure/wall BCs across time steps'
for tag in 20 10 05; do
 case "$tag" in 20) dt="0.02";; 10) dt="0.01";; 05) dt="0.005";; esac
 echo "M1B-1g-2 running dt=$dt"
 timeout 180s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$work/dt$tag" > "$work/solve$tag.log" 2>&1 || {
  tail -100 "$work/solve$tag.log"; exit 31;
 }
done
python3 /workspace/tests/m1b1/check_warped_open_pressure_timestep.py \
 "$work/check.log" \
 "$work/dt20" "$work/solve20.log" \
 "$work/dt10" "$work/solve10.log" \
 "$work/dt05" "$work/solve05.log"
echo 'PASS: M1B-1g-2 full OpenFOAM v2606 warped duct temporal convergence'
