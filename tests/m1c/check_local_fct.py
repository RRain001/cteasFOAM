#!/usr/bin/env python3
"""Independent M1C-19 true-cell facewise MULES source/Bq/positivity audit."""
import math, re, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_fixed_receptors import audit
from check_convection_scheme import EXTREMA

NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
PAT=re.compile(r'CTEAS_M1C19 species=(HT|HTO) face_lambda_min=('+NUM+r') face_lambda_max=('+NUM+r') limited_faces=(\d+) activity_min_Bq_m3=('+NUM+r')')

def checked(case,n,mode):
    log=(case/'solve.log').read_text()
    data=audit(case,n,solves_per_species_step=1)
    rows=PAT.findall(log)
    if mode=='baseline':
        if rows:raise SystemExit('FAIL: M1C19 default upwind unexpectedly used facewise FCT')
    else:
        if len(rows)!=160:raise SystemExit(f'FAIL: M1C19 requires 160 true limiter records, got {len(rows)}')
        limited={'HT':0,'HTO':0}
        minima={'HT':1.0,'HTO':1.0}
        for j,(sp,mn,mx,count,cmin) in enumerate(rows):
            if sp!=('HT' if j%2==0 else 'HTO'):raise SystemExit('FAIL: M1C19 wrong native species order')
            mn,mx,cmin=float(mn),float(mx),float(cmin);count=int(count)
            if not all(map(math.isfinite,(mn,mx,cmin))) or not (0<=mn<=mx<=1) or count<0 or cmin < -1e-8:
                raise SystemExit('FAIL: M1C19 invalid real limiter/negative native activity')
            minima[sp]=min(minima[sp],mn)
            limited[sp]+=int(count>0)
        for sp in ('HT','HTO'):
            print(f'M1C19_FACE_LIMITER N={n} species={sp} min_lambda={minima[sp]:.16g} limited_steps={limited[sp]}/80')
    extrema=[tuple(map(float,m)) for m in EXTREMA.findall(log)]
    if len(extrema)!=80 or any(not all(map(math.isfinite,row)) or row[0]<-1e-8 or row[2]<-1e-8 for row in extrema):
        raise SystemExit('FAIL: M1C19 real per-step HT/HTO extrema invalid')
    print(f'M1C19_NATIVE N={n} mode={mode} max_signed_Bq={data["max_budget_Bq"]:.16g}')
    return data

def compare(root,n):
    baseline=checked(root/f'{n}_base',n,'baseline')
    face=checked(root/f'{n}_face',n,'face')
    for probe in ('R1','R2','R3'):
        for sp in ('HT','HTO'):
            x=baseline['receptors'][probe]['species'][sp]['relative_error']
            y=face['receptors'][probe]['species'][sp]['relative_error']
            print(f'M1C19_RECEIVER N={n} {probe} {sp} upwind_relative={x:.12g} face_fct_relative={y:.12g}')
    print(f'PASS: M1C19 native {n}^3 face-local MULES source/Green/activity/positivity scope')

if __name__=='__main__':
    if len(sys.argv)==4 and sys.argv[1]=='--compare':compare(Path(sys.argv[2]),int(sys.argv[3]))
    else:raise SystemExit('usage: --compare ROOT 8|16|24')
