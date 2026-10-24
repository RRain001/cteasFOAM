#!/usr/bin/env bash
set -eo pipefail
set +e
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
[[ "$WM_PROJECT_VERSION" == v2606 || "$WM_PROJECT_VERSION" == 2606 ]] || exit 18
set -u
work=$(mktemp -d /tmp/cteasfoam-m1c16-XXXXXX)
trap 'rc=$?; if ((rc!=0)); then echo "FAIL: M1C-16 exit=$rc"; for f in "$work"/build.log "$work"/*/solve.log "$work"/negative-*.log; do if [[ -f "$f" ]]; then echo "== $f =="; tail -42 "$f"; fi; done; fi; rm -rf "$work"' EXIT
export FOAM_USER_APPBIN="$work/bin"
mkdir -p "$FOAM_USER_APPBIN"
cp -a /workspace/src/atmosphere "$work/atmosphere"
(cd "$work/atmosphere" && wclean && wmake > "$work/build.log" 2>&1) || { tail -95 "$work/build.log";exit 21; }
test -x "$FOAM_USER_APPBIN/cteasAtmosphereFoam" || exit 22
for n in 8 16 24; do
 for scheme in upwind limited; do
  case_dir="$work/$n"_"$scheme"
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py make "$case_dir" zero "$n"x"$n"x"$n"
  python3 /workspace/tests/m1c/make_3d_analytic.py prepare "$case_dir" advect
  python3 /workspace/tests/m1c/make_gaussian_green.py "$case_dir" normal
  python3 - "$case_dir" "$scheme" <<'PYCASE'
from pathlib import Path
import sys
case=Path(sys.argv[1])
scheme=sys.argv[2]
p=case/"constant/transportProperties"
s=p.read_text()
assert "auditActivityExtrema" not in s
p.write_text(s+"auditActivityExtrema true;\n")
fv=case/"system/fvSchemes"
s=fv.read_text()
if scheme=="limited":
    for species in ("HT","HTO"):
        old=f"div(phi,{species}) Gauss upwind;"
        assert s.count(old)==1
        s=s.replace(old,f"div(phi,{species}) Gauss limitedLinear 1;")
fv.write_text(s)
print(f"M1C16_NATIVE_CASE N={sys.argv[1]} scheme={scheme}")
PYCASE
  blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1 || { tail -60 "$case_dir/blockMesh.log";exit 23; }
  python3 /workspace/tests/m1b1/make_warped_open_pressure_case.py warp "$case_dir" zero "$n"x"$n"x"$n"
  checkMesh -case "$case_dir" -allGeometry -allTopology > "$case_dir/mesh-check.log" 2>&1 || { tail -60 "$case_dir/mesh-check.log";exit 24; }
  grep -q 'Mesh OK.' "$case_dir/mesh-check.log" || exit 24
  postProcess -case "$case_dir" -func writeCellCentres -time 0 > "$case_dir/centres.log" 2>&1 || { tail -45 "$case_dir/centres.log";exit 25; }
  postProcess -case "$case_dir" -func writeCellVolumes -time 0 > "$case_dir/volumes.log" 2>&1 || { tail -45 "$case_dir/volumes.log";exit 26; }
  timeout 480s "$FOAM_USER_APPBIN/cteasAtmosphereFoam" -case "$case_dir" > "$case_dir/solve.log" 2>&1 || { tail -85 "$case_dir/solve.log";exit 31; }
  if [[ "$scheme" == limited ]]; then
    python3 /workspace/tests/m1c/check_convection_scheme.py --single-diagnostic "$case_dir" "$n" "$scheme"
    if [[ "$n" == 8 ]]; then
      if python3 /workspace/tests/m1c/check_convection_scheme.py --single "$case_dir" "$n" "$scheme" > "$work/negative-physical-limiter.log" 2>&1; then
        echo "FAIL: M1C16 physically unsafe 8^3 limitedLinear passed strict acceptance"; exit 35
      fi
      grep -q 'FAIL: M1C-16 native real 80-step HT negative concentration' "$work/negative-physical-limiter.log" || { cat "$work/negative-physical-limiter.log"; exit 36; }
      echo "M1C16_SCIENTIFIC_REJECTED real 8^3 limitedLinear1 native HT intermediate concentration is negative"
    fi
  else
    python3 /workspace/tests/m1c/check_convection_scheme.py --single "$case_dir" "$n" "$scheme"
  fi
 done
done
checker=/workspace/tests/m1c/check_convection_scheme.py
python3 "$checker" --matrix-diagnostic "$work"
native="$work/16_upwind"
cp "$native/0.4/HT" "$work/pristine-final-HT"
python3 - "$native/0.4/HT" <<'PYHT'
from pathlib import Path
import sys,re
p=Path(sys.argv[1]);s=p.read_text()
m=re.search(r"(internalField\s+nonuniform\s+List<scalar>\s+4096\s*\(\s*)([-+0-9.eE]+)",s)
assert m
p.write_text(s[:m.start(2)]+f"{float(m.group(2))+1000000:.17g}"+s[m.end(2):])
PYHT
if python3 "$checker" --single "$native" 16 upwind > "$work/negative-field.log" 2>&1; then echo "FAIL: M1C16 real final cell tamper accepted";exit 40;fi
grep -q 'FAIL: M1C-16 real final HT field differs from logged native extrema' "$work/negative-field.log" || { cat "$work/negative-field.log";exit 41; }
mv "$work/pristine-final-HT" "$native/0.4/HT"
echo "PASS: M1C16 native HT final physical cell corruption rejected"
cp "$native/solve.log" "$work/pristine-solve"
python3 - "$native/solve.log" <<'PYEXT'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
s,n=re.subn(r"(extrema_HT_min_Bq_m3=)([-+0-9.eE]+)",r"\g<1>-1",s,count=1)
if n!=1:raise SystemExit("FAIL: native real extrema missing")
p.write_text(s)
PYEXT
if python3 "$checker" --single "$native" 16 upwind > "$work/negative-extrema.log" 2>&1; then echo "FAIL: M1C16 negative native extrema tamper accepted";exit 42;fi
grep -q 'FAIL: M1C-16 native real 80-step HT negative concentration' "$work/negative-extrema.log" || { cat "$work/negative-extrema.log";exit 43; }
mv "$work/pristine-solve" "$native/solve.log"
echo "PASS: M1C16 native per-step negative concentration record rejected"
cp "$native/system/fvSchemes" "$work/pristine-fv"
sed -i 's/div(phi,HT) Gauss upwind;/div(phi,HT) Gauss bogusScheme;/' "$native/system/fvSchemes"
if python3 "$checker" --single "$native" 16 upwind > "$work/negative-scheme.log" 2>&1; then echo "FAIL: M1C16 bad input scheme tamper accepted";exit 44;fi
grep -q 'FAIL: M1C-16 frozen physical advection scheme changed HT' "$work/negative-scheme.log" || { cat "$work/negative-scheme.log";exit 45; }
mv "$work/pristine-fv" "$native/system/fvSchemes"
echo "PASS: M1C16 unsupported actual HT scheme in dictionary rejected before output interpretation"
cp "$native/solve.log" "$work/pristine-solve"
python3 - "$native/solve.log" <<'PYEM'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
s,n=re.subn(r"(emitted_HT_Bq=)([-+0-9.eE]+)",lambda m:m.group(1)+str(float(m.group(2))+1000),s,count=1)
assert n==1
p.write_text(s)
PYEM
if python3 "$checker" --single "$native" 16 upwind > "$work/negative-emitted.log" 2>&1; then echo "FAIL: M1C16 emitted native Bq tamper accepted";exit 46;fi
grep -q 'FAIL: M1C-14 genuine Gaussian finite-window emission mismatch' "$work/negative-emitted.log" || { cat "$work/negative-emitted.log";exit 47; }
mv "$work/pristine-solve" "$native/solve.log"
echo "PASS: M1C16 source window signed Bq audit catches native log tamper"
echo "PASS: M1C-16 6 native 3D Gaussian upwind/limitedLinear diagnostic simulations; physically negative high-order candidate explicitly REJECTED; 4 strict negatives"
