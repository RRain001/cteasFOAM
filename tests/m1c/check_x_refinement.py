#!/usr/bin/env python3
"""Independent M1C-10 x-only native source/quadrature versus downstream receptors."""
from pathlib import Path
import json
import math
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_receiver_diagnostics import CASES, probe
from check_gaussian_source import field

MODES={
    "baseline":(8,20,6,.02,20),
    "space_x":(12,20,6,.02,20),
    "x16":(16,20,6,.02,20),
    "x24":(24,20,6,.02,20),
    "x32":(32,20,6,.02,20)
}
CASES.update(MODES)
C=(.38,.62,.55)
SIGMA=.13
SOURCE_SHA="7b402a4420aa255db9c9183c416917ed5015058a"
R=Path(__file__).with_name("m1c9_xref_exact_observations.json")

def fail(s):raise SystemExit("FAIL: M1C-10 "+s)

def exact_unit_truncated_normal(m,s):
    """Analytic integrals of exp(-(x-m)^2/(2s^2)), x in [0,1]."""
    aa=-m
    bb=1-m
    a=aa/(math.sqrt(2)*s)
    b=bb/(math.sqrt(2)*s)
    q0=s*math.sqrt(math.pi/2)*(math.erf(b)-math.erf(a))
    ea=math.exp(-aa*aa/(2*s*s))
    eb=math.exp(-bb*bb/(2*s*s))
    if q0<=0:fail("continuous Gaussian integral vanishes")
    mean_offset=s*s*(ea-eb)/q0
    second=s*s*(1+(aa*ea-bb*eb)/q0)
    if not (math.isfinite(mean_offset) and second>0):fail("invalid analytic source moments")
    return q0, m+mean_offset, second

def native_geometry(case,n):
    c=field(case/"0/C",n,"vector","0 1 0 0 0 0 0")
    v=field(case/"0/V",n,"scalar","0 3 0 0 0 0 0")
    if len(c)!=n or len(v)!=n:fail("native geometry count incorrect")
    if min(v)<=0 or abs(sum(v)-1)>1e-8:fail("native mesh physical volume error")
    if any(any(not (0<=q<=1) for q in coords) for coords in c):
        fail("native cell centers outside unit cube")
    return c,v

def moments(case,n):
    c,v=native_geometry(case,n)
    source_weights=[vol*math.exp(-sum((q-m)**2 for q,m in zip(center,C))
                                      /(2*SIGMA**2))
                    for center,vol in zip(c,v)]
    z=sum(source_weights)
    if not (math.isfinite(z) and z>0):fail("source weights missing/nonfinite")
    w=[t/z for t in source_weights]
    if abs(sum(w)-1)>1e-12:fail("Gaussian cell-volume weight not normalized")
    centroid=[sum(p*pos[k] for p,pos in zip(w,c)) for k in range(3)]
    second=[[sum(p*(pos[i]-C[i])*(pos[j]-C[j]) for p,pos in zip(w,c))
             for j in range(3)] for i in range(3)]
    exact=[exact_unit_truncated_normal(m,SIGMA) for m in C]
    exact_q=math.prod(a[0] for a in exact)
    expected_centroid=[a[1] for a in exact]
    mean_delta=[a[1]-m for m,a in zip(C,exact)]
    expected_second=[
      [(exact[i][2] if i==j else mean_delta[i]*mean_delta[j])
       for j in range(3)] for i in range(3)
    ]
    return {
      "normalizer_native_m3":z,
      "normalizer_analytic_m3":exact_q,
      "source_normalizer_relative_error":(z-exact_q)/exact_q,
      "discrete_centroid_m":centroid,
      "continuous_centroid_m":expected_centroid,
      "centroid_error_m":[a-b for a,b in zip(centroid,expected_centroid)],
      "discrete_second_moment_m2":second,
      "continuous_second_moment_m2":expected_second,
      "second_moment_error_m2":[
         [second[i][j]-expected_second[i][j] for j in range(3)]
         for i in range(3)
      ],
      "source_weight_sum":sum(w)
    }

def single(case,name):
    if name not in MODES:fail("unknown fixed native resolution")
    native=probe(case,name)  # Frozen M1C-8 audit; real fields, budgets and geometry
    nm=native["native_qa"]["cells"]
    if nm!=math.prod(MODES[name][:3]) or native["native_qa"]["dt_s"]!=.02:
        fail("native physical mesh/time changed")
    mom=moments(case,nm)
    if max(abs(a-b) for a,b in zip(native["centroid_m"],
                                 mom["discrete_centroid_m"]))>2e-10:
        fail("independent native Gaussian centroid audit disagrees with M1C-9")
    return {"native_qa":native["native_qa"],
            "source_quadrature":mom,
            "receptor_Bq_m3":native["receptor_Bq_m3"],
            "nearest_receiver_m":native["nearest_receiver_m"]}

def compare(root):
    ref=json.loads(R.read_text())
    if ref["frozen_parent_sha"]!=SOURCE_SHA or set(ref["reference"])!={"baseline","space_x"}:
        fail("missing or mismatched independently frozen M1C-9 reference")
    outcomes={name:single(root/name,name) for name in MODES}
    refdiff=0
    for name,old in ref["reference"].items():
        case=outcomes[name]
        for metric,expected in old["receptor_Bq_m3"].items():
            diff=abs(case["receptor_Bq_m3"][metric]-expected)
            refdiff=max(diff,refdiff)
            if diff>1e-5:fail(f"{name}: prior independent R1–R3 native concentration changed {metric} delta={diff}")
        for a,b in zip(case["source_quadrature"]["discrete_centroid_m"],old["native_centroid_m"]):
            if abs(a-b)>2e-10:fail("parent source centroid changed")
    trend={}
    for sp in ("HT","HTO"):
        key=f"R3_{sp}_w0.10_dx+0.000"
        values=[outcomes[name]["receptor_Bq_m3"][key] for name in MODES]
        differences=[values[i+1]-values[i] for i in range(4)]
        last_decreases=abs(differences[3])<abs(differences[2])<abs(differences[1])
        last_change_relative=abs(differences[3])/max(abs(values[-1]),1e-9)
        status="TREND_DECREASED" if last_decreases else "NOT_DEMONSTRATED"
        trend[sp]={"cells_by_nx":[MODES[name][0] for name in MODES],
                   "native_R3_Bq_m3":values,
                   "signed_neighbor_deltas_Bq_m3":differences,
                   "fine_step_relative_change":last_change_relative,
                   "trend":status}
        print(f"M1C10_R3_{sp} native="+",".join(f"{v:.10g}" for v in values)
              +" signed_deltas="+",".join(f"{v:.10g}" for v in differences)
              +f" relative_last={last_change_relative:.9g} {status}")
    for name,case in outcomes.items():
        m=case["source_quadrature"]
        print(f"M1C10_SOURCE {name} NX={MODES[name][0]}"
              f" integrated_HT_Bq={case['native_qa']['HT_released_Bq']:.10g}"
              f" integrated_HTO_Bq={case['native_qa']['HTO_released_Bq']:.10g}"
              f" quadrature_normalizer_relative={m['source_normalizer_relative_error']:.10g}"
              f" centroid_error_m={m['centroid_error_m']}"
              f" second_diag_error_m2={[m['second_moment_error_m2'][i][i] for i in range(3)]}")
    report={"scope":"M1C-10 synthetic native one-dimensional mesh refinement; not a real-world atmospheric domain",
            "parent_source_sha":SOURCE_SHA,"prior_max_abs_error_Bq_m3":refdiff,
            "source_center_m":C,"source_sigma_m":SIGMA,
            "x_only_meshes":[MODES[name][:3] for name in MODES],
            "outcomes":outcomes,"R3_refinement":trend,
            "conclusion":"LIMITED_SOURCE_QUADRATURE_AND_X_SENSITIVITY_AUDIT_ONLY"}
    (root/"m1c10-report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("M1C10_REPORT_JSON="+json.dumps(report,sort_keys=True,separators=(",",":")))
    print(f"PASS: M1C-10 5 independent exact-native x-only cases, source analytic moments and parent reference maxError={refdiff:.10g} Bq/m3; R3 trends classified separately")
    return report

if __name__=="__main__":
    if len(sys.argv)==4 and sys.argv[1]=="--single":
        single(Path(sys.argv[2]),sys.argv[3])
        print("PASS: M1C-10 native x-only signed source and boundary Bq validation")
    elif len(sys.argv)==3 and sys.argv[1]=="--matrix":
        compare(Path(sys.argv[2]))
    else:
        raise SystemExit("usage: check_x_refinement.py --single CASE MODE | --matrix ROOT")
