#!/usr/bin/env bash
# Run after 8/16/24 and MPI/open checks; preserve actual native evidence.
cp "$work/8_face/solve.log" "$work/8-face-pristine.log"
python3 - "$work/8_face/solve.log" <<'PY'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
s,n=re.subn(r'(CTEAS_M1C19 species=HT face_lambda_min=)([-+0-9.eE]+)',r'\g<1>-0.75',s,count=1)
assert n==1
p.write_text(s)
PY
if python3 /workspace/tests/m1c/check_local_fct.py --compare "$work" 8 > "$work/forged-lambda.log" 2>&1; then
 echo 'FAIL: M1C19 forged negative local face lambda accepted';exit 71
fi
grep -q 'FAIL: M1C19 invalid real limiter' "$work/forged-lambda.log" || { tail -45 "$work/forged-lambda.log";exit 72; }
mv "$work/8-face-pristine.log" "$work/8_face/solve.log"
echo 'PASS: M1C-19 independently rejected forged negative real-face lambda record'

bad="$work/8_face_negative_initial"
cp -a "$work/8_face" "$bad"
python3 - "$bad" <<'PY'
from pathlib import Path
import sys
c=Path(sys.argv[1]); p=c/'0/HT';t=p.read_text()
assert t.count('internalField uniform 0;')==1
p.write_text(t.replace('internalField uniform 0;','internalField uniform -1;'))
props=c/'constant/transportProperties';s=props.read_text()
assert s.count('enforceNonnegativeActivity true;')==1
props.write_text(s.replace('enforceNonnegativeActivity true;','enforceNonnegativeActivity false;'))
PY
if timeout 40s "$exe" -case "$bad" > "$work/negative-initial.log" 2>&1; then
 echo 'FAIL: M1C19 negative initial HT admitted with optional M1C17 guard disabled';exit 73
fi
grep -q 'CTEAS_M1C_POSITIVITY_FAILURE species=HT phase=initial' "$work/negative-initial.log" || { tail -50 "$work/negative-initial.log";exit 74; }
if grep -q '^Time = ' "$work/negative-initial.log"; then
 echo 'FAIL: M1C19 negative input advanced native solver timestep'; exit 75
fi
echo 'PASS: M1C-19 rejects genuine negative initial HT before time loop even with optional diagnostic guard disabled'
