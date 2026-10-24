#!/usr/bin/env python3
"""M1C-12: configure matched source-free and pulse 3D eigenmode twins."""
from pathlib import Path
import sys

DT={"coarse":("0.01",20),"normal":("0.005",40),"fine":("0.0025",80)}
START=.073
END=.137
Q_HT=60000
Q_HTO=12000

def apply(case,source,step):
    if source not in ("control","pulse") or step not in DT:
        raise SystemExit("FAIL: M1C-12 unsupported source/timestep")
    raw,count=DT[step]
    ctrl=case/"system/controlDict"
    s=ctrl.read_text()
    if step!="normal":
        if s.count("deltaT 0.005;")!=1 or s.count("writeInterval 40;")!=1:
            raise SystemExit("FAIL: M1C-12 no frozen base clock")
        ctrl.write_text(s.replace("deltaT 0.005;",f"deltaT {raw};")
                         .replace("writeInterval 40;",f"writeInterval {count};"))
    p=case/"constant/transportProperties"
    s=p.read_text()
    if any(v in s for v in ("sourceHT_Bq_s","sourceHTO_Bq_s","sourceStart","sourceEnd","sourceShape")):
        raise SystemExit("FAIL: M1C-12 inherited source is not clean")
    if source=="pulse":
        with p.open("a") as f:
            f.write(
                "sourceShape uniform;\n"
                f"sourceHT_Bq_s [0 0 -2 0 0 0 0] {Q_HT};\n"
                f"sourceHTO_Bq_s [0 0 -2 0 0 0 0] {Q_HTO};\n"
                f"sourceStart [0 0 1 0 0 0 0] {START:.17g};\n"
                f"sourceEnd [0 0 1 0 0 0 0] {END:.17g};\n"
            )
    print(f"M1C-12 prepared {source}/{step} steps={count} dt={raw} "
          f"window=[{START},{END}]s uniform_release_HT={Q_HT*(END-START) if source=='pulse' else 0:.10g}"
          f" HTO={Q_HTO*(END-START) if source=='pulse' else 0:.10g}")

if __name__=="__main__":
    if len(sys.argv)!=4:raise SystemExit("usage: make_pulsed_analytic.py CASE control|pulse coarse|normal|fine")
    apply(Path(sys.argv[1]),sys.argv[2],sys.argv[3])
