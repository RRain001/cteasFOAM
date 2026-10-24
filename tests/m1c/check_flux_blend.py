#!/usr/bin/env python3
"""Independent M1C-18 global convex flux-blend evidence; no implicit scientific acceptance."""
import math
import re
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_convection_scheme import single, EXTREMA
from check_fixed_receptors import audit
from check_gaussian_green import field, DIM
ALPHA=re.compile(r'CTEAS_M1C18 species=(HTO|HT) global_alpha=([-+0-9.eE]+) limited=([01])')

def audit_case(case,n,mode):
    log=(case/'solve.log').read_text()
    data=(single(case,n,'upwind') if mode=='base'
          else audit(case,n,solves_per_species_step=2))
    extrema=[tuple(map(float,m)) for m in EXTREMA.findall(log)]
    if len(extrema)!=80 or any(not all(math.isfinite(x) for x in row) for row in extrema):
        raise SystemExit('FAIL: M1C18 missing/nonfinite genuine 80-step species extrema')
    for species,col in (('HT',0),('HTO',2)):
        if min(row[col] for row in extrema)<-1e-8:
            raise SystemExit('FAIL: M1C18 actual negative intermediate species concentration')
        final=field(case/'0.4'/species,n**3,'scalar',DIM)
        if max(abs(min(final)-extrema[-1][col]),abs(max(final)-extrema[-1][col+1]))>2e-5:
            raise SystemExit('FAIL: M1C18 recorded native extrema disagree with physical final cells')
    rows=ALPHA.findall(log)
    if mode=='blend':
        if len(rows)!=160:raise SystemExit(f'FAIL: M1C18 expected genuine 160 stepwise blend records, got {len(rows)}')
        steps={s:[] for s in ('HT','HTO')}
        for j,(species,a,limited) in enumerate(rows):
            if species!= ('HT' if j%2==0 else 'HTO'):
                raise SystemExit('FAIL: M1C18 wrong species sequence')
            val=float(a)
            if not math.isfinite(val) or not (0<=val<=1):
                raise SystemExit('FAIL: M1C18 invalid native alpha')
            if bool(int(limited)) != (val < 1-1e-12):
                raise SystemExit('FAIL: M1C18 alpha / limiter status inconsistent')
            steps[species].append(val)
        for species,values in steps.items():
            print(f'M1C18_ALPHA N={n} {species} min={min(values):.16g} max={max(values):.16g} constrained={sum(v<1-1e-12 for v in values)}/80')
    elif rows:
        raise SystemExit('FAIL: M1C18 default upwind unexpectedly invoked blending')
    print(f'M1C18_CASE N={n} mode={mode} max_signed_Bq_residual={data.get("native_activity_budget_max_Bq",data.get("max_budget_Bq")):.12g}')
    return data

def compare(root,n):
    base=audit_case(root/f'{n}_base',n,'base')
    blend=audit_case(root/f'{n}_blend',n,'blend')
    for name in ('R1','R2','R3'):
        for species in ('HT','HTO'):
            a=base['receptors'][name]['species'][species]['relative_error']
            b=blend['receptors'][name]['species'][species]['relative_error']
            print(f'M1C18_RECEIVER N={n} {name} {species} upwind_error={a:.12g} global_blend_error={b:.12g} delta={b-a:.12g}')
    print(f'PASS: M1C-18 {n}^3 real high/low finite-volume convex blending, strict native positivity + 80-step Bq budgets (NOT local FCT/MULES)')

if __name__=='__main__':
    if len(sys.argv)==4 and sys.argv[1]=='--compare':compare(Path(sys.argv[2]),int(sys.argv[3]))
    else:raise SystemExit('usage: --compare ROOT 8|16|24')
