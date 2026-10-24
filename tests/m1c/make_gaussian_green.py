#!/usr/bin/env python3
"""M1C-13 physical native warped periodic 3D localized Gaussian pulse configuration."""
from pathlib import Path
import sys
from make_3d_analytic import change

CTR=(.42,.55,.50)
SIG=.08
U=(1,.05,.03)
S=(.073,.137)
RATES=(60000,12000)

def configure(case,step):
    if step not in ("normal","coarse"):raise SystemExit("FAIL: M1C-13 unknown native clock")
    u=case/"0/U"
    change(u,"internalField uniform (0.2 0.1 0.07);",
           "internalField uniform (1 0.05 0.03);")
    for name,value in (("HT","1000000"),("HTO","200000")):
        change(case/"0"/name,f"internalField uniform {value};","internalField uniform 0;")
    c=case/"system/controlDict"
    change(c,"endTime 0.2;","endTime 0.4;")
    change(c,"writeInterval 40;",f"writeInterval {80 if step=='normal' else 40};")
    if step=="coarse":change(c,"deltaT 0.005;","deltaT 0.01;")
    props=case/"constant/transportProperties"
    x=props.read_text()
    if any(v in x for v in ("sourceShape","sourceHT_Bq_s","sourceEnd")):
        raise SystemExit("FAIL: M1C-13 inherited source pollution")
    props.write_text(x+
        "sourceShape gaussian;\n"
        "sourceCenter [0 1 0 0 0 0 0] (0.42 0.55 0.5);\n"
        "sourceSigma [0 1 0 0 0 0 0] 0.08;\n"
        "sourceHT_Bq_s [0 0 -2 0 0 0 0] 60000;\n"
        "sourceHTO_Bq_s [0 0 -2 0 0 0 0] 12000;\n"
        "sourceStart [0 0 1 0 0 0 0] 0.073;\n"
        "sourceEnd [0 0 1 0 0 0 0] 0.137;\n")
    print(f"M1C13_CONFIG native source shape Gaussian center={CTR} sigma={SIG} U={U} "
          f"window={S} release_HT={RATES[0]*(S[1]-S[0]):.12g} "
          f"release_HTO={RATES[1]*(S[1]-S[0]):.12g} step={step}")

if __name__=="__main__":
    if len(sys.argv)!=3:raise SystemExit("usage: make_gaussian_green.py CASE coarse|normal")
    configure(Path(sys.argv[1]),sys.argv[2])
