#!/usr/bin/env python3
"""M1C-17: classify native OpenFOAM positivity failures; never certify unsafe output."""
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_gaussian_green import field, DIM
from check_convection_scheme import single

NUM = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
FAIL = re.compile(r'CTEAS_M1C_POSITIVITY_FAILURE species=(HT|HTO) phase=(initial|postSolve) minimum_Bq_m3=('+NUM+r') threshold_Bq_m3=('+NUM+r')')


def verify_rejection(case, phase):
    log = (case / 'solve.log').read_text()
    matches = FAIL.findall(log)
    if len(matches) != 1:
        raise SystemExit(f'FAIL: M1C-17 expected exactly one genuine native positivity failure, got {len(matches)}')
    species, observed_phase, minimum, threshold = matches[0]
    minimum, threshold = float(minimum), float(threshold)
    if observed_phase != phase or species not in ('HT', 'HTO') or not math.isfinite(minimum) or not math.isfinite(threshold):
        raise SystemExit('FAIL: M1C-17 wrong species/phase/nonfinite minimum')
    if abs(threshold + 1e-8) > 1e-15 or minimum >= threshold:
        raise SystemExit('FAIL: M1C-17 actual minimum was not below frozen positivity tolerance')
    if 'Completed mode=' in log:
        raise SystemExit('FAIL: M1C-17 failed-positive native solver nevertheless completed')
    steps = len(re.findall(r'^Time = [0-9]', log, re.M))
    if phase == 'initial' and steps:
        raise SystemExit('FAIL: M1C-17 bad initial field admitted to time loop')
    if phase == 'postSolve' and not (0 < steps < 80):
        raise SystemExit(f'FAIL: M1C-17 expected real transient negative before step 80, got {steps}')
    print(f'M1C17_EXPECTED_NATIVE_REJECTION species={species} phase={phase} min={minimum:.16g} step_count={steps}')


def verify_pair(first, second):
    for case in (first, second):
        log = (case / 'solve.log').read_text()
        if FAIL.search(log) or 'Completed mode=' not in log:
            raise SystemExit('FAIL: M1C-17 upwind run failed or was incomplete')
        single(case, 8, 'upwind')  # preserves all existing Green, source, budget, extrema checks
    for species in ('HT', 'HTO'):
        a = field(first / '0.4' / species, 512, 'scalar', DIM)
        b = field(second / '0.4' / species, 512, 'scalar', DIM)
        mismatch = max(abs(x-y) for x,y in zip(a,b))
        if mismatch > 1e-11:
            raise SystemExit(f'FAIL: M1C-17 guarding changed physical {species} field by {mismatch}')
        print(f'M1C17_UPWIND_PARITY species={species} max_delta_Bq_m3={mismatch:.12g}')
    print('PASS: M1C-17 nonnegative opt-in guard preserves physical upwind fields, signed budgets and source totals')


if __name__ == '__main__':
    if len(sys.argv) == 4 and sys.argv[1] == '--reject':
        verify_rejection(Path(sys.argv[2]), sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == '--compare':
        verify_pair(Path(sys.argv[2]), Path(sys.argv[3]))
    else:
        raise SystemExit('usage: --reject CASE initial|postSolve | --compare BASE GUARDED')
