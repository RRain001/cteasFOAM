#!/usr/bin/env bash
# Genuine native C++ negative controls using the same compiled solver and real case.
for testCase in bad-low missing-low missing-sngrad wrong-sngrad invalid-mode incompatible-mode; do
  c="$work/fail-$testCase"
  cp -a "$work/normal" "$c"
  python3 - "$c" "$testCase" <<'PYERROR'
from pathlib import Path
import sys
c=Path(sys.argv[1]);mode=sys.argv[2]
p=c/"system/fvSchemes";s=p.read_text()
props=c/"constant/transportProperties";q=props.read_text()
if mode=="bad-low":
    old="laplacian(DHT,HTLowDiff) Gauss linear uncorrected;"
    assert s.count(old)==1;s=s.replace(old,old.replace("uncorrected","corrected"))
if mode=="missing-low":
    old="laplacian(DHT,HTLowDiff) Gauss linear uncorrected;"
    assert s.count(old)==1;s=s.replace(old,"")
if mode=="missing-sngrad":
    old="snGrad(HTHighDiff) corrected;"
    assert s.count(old)==1;s=s.replace(old,"")
if mode=="wrong-sngrad":
    old="snGrad(HTHighDiff) corrected;"
    assert s.count(old)==1;s=s.replace(old,"snGrad(HTHighDiff) uncorrected;")
if mode=="invalid-mode":
    assert q.count("activityDiffusionCorrection faceMules;")==1
    q=q.replace("activityDiffusionCorrection faceMules;","activityDiffusionCorrection nonexistent;")
if mode=="incompatible-mode":
    assert q.count("activityConvectionBlend faceMules;")==1
    q=q.replace("activityConvectionBlend faceMules;","activityConvectionBlend none;")
p.write_text(s);props.write_text(q)
PYERROR
  case "$testCase" in
    bad-low) marker=M1C21_REQUIRED_POSITIVE_LOW_SCHEME_WRONG;;
    missing-low) marker=M1C21_REQUIRED_POSITIVE_LOW_SCHEME_MISSING;;
    missing-sngrad) marker=M1C21_REQUIRED_NAMED_SNGRAD_MISSING;;
    wrong-sngrad) marker=M1C21_REQUIRED_NAMED_SNGRAD_WRONG;;
    invalid-mode) marker=M1C21_UNSUPPORTED_DIFFUSION_MODE;;
    incompatible-mode) marker=M1C21_REQUIRES_FACE_MULES_HT_HTO;;
  esac
  if timeout 45s "$exe" -case "$c" > "$work/fail-$testCase.log" 2>&1;then
    echo "FAIL: M1C21 malformed $testCase configuration was accepted";exit 101
  fi
  grep -q "$marker" "$work/fail-$testCase.log" || {
    echo "FAIL: wrong reason for M1C21 $testCase";tail -45 "$work/fail-$testCase.log";exit 102;
  }
  if grep -q '^Time = ' "$work/fail-$testCase.log";then
    echo "FAIL: M1C21 $testCase advanced past configuration gate";exit 103
  fi
  echo "PASS: M1C21 authentic native fail-closed $testCase rejected before time step"
done
# Independently reject doctored native correction evidence.
cp "$work/normal/solve.log" "$work/m1c21-pristine.log"
python3 - "$work/normal/solve.log" <<'PYLOG'
import re,sys
from pathlib import Path
p=Path(sys.argv[1]);s=p.read_text()
s,n=re.subn(
    r'(CTEAS_M1C21_DIFFUSION species=HT max_abs_high_minus_low_Bq_s=)([-+0-9.eE]+)',
    r'\g<1>-5.25',s,count=1)
assert n==1
p.write_text(s)
PYLOG
if python3 /workspace/tests/m1c/check_bounded_diffusion.py "$work/normal" 8 80 0.005 > "$work/tampered.log" 2>&1;then
 echo "FAIL: M1C21 negative signed diffusive face amplitude doctored accepted";exit 104
fi
grep -q "FAIL: M1C21 invalid correction face amplitude/count" "$work/tampered.log" || {
 tail -50 "$work/tampered.log";exit 105;
}
mv "$work/m1c21-pristine.log" "$work/normal/solve.log"
echo "PASS: M1C21 independent correction-log negative-amplitude tamper rejected"
