#!/usr/bin/env python3
"""M1C-16: exact native finite-volume scheme, 80-step extrema and fixed Green probes."""
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_fixed_receptors import audit
from check_gaussian_green import field, DIM, POS

NUM=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
EXTREMA=re.compile(
    r"CTEAS_M1C extrema_HT_min_Bq_m3=("+NUM+r")"
    r" extrema_HT_max_Bq_m3=("+NUM+r")"
    r" extrema_HTO_min_Bq_m3=("+NUM+r")"
    r" extrema_HTO_max_Bq_m3=("+NUM+r")")

SCHEMES={"upwind":"Gauss upwind",
         "limited":"Gauss limitedLinear 1"}
NAMES=("HT","HTO")

def fail(s):
    raise SystemExit("FAIL: M1C-16 "+s)

def single(case,n,scheme,diagnostic=False):
    if n not in (8,16,24) or scheme not in SCHEMES:
        fail("unsupported frozen native mesh or advection scheme")
    props=(case/"constant/transportProperties").read_text()
    if "auditActivityExtrema true;" not in props:
        fail("native per-step activity extrema switch not enabled")
    if "sourceSpatialQuadrature" in props or "writeSourceProfile" in props:
        fail("legacy Gaussian cellCentre spatial quadrature was altered")
    fv=(case/"system/fvSchemes").read_text()
    for sp in NAMES:
        requested=f"div(phi,{sp}) {SCHEMES[scheme]};"
        if requested not in fv:
            fail("frozen physical advection scheme changed "+sp)
    if "div(phi,c) Gauss upwind;" not in fv:
        fail("passive scalar convection scheme was modified")
    log=(case/"solve.log").read_text()
    moments=[tuple(map(float,m)) for m in EXTREMA.findall(log)]
    if len(moments)!=80:
        fail("genuine per-step native species extrema absent or incomplete")
    worst=0.
    maxvals={}
    minvals={}
    scientific_rejections=[]
    for sp,idx in (("HT",0),("HTO",2)):
        minimums=[step[idx] for step in moments]
        maximums=[step[idx+1] for step in moments]
        if not all(math.isfinite(v) for v in minimums+maximums):
            fail("native activity extrema nonfinite")
        if any(mn>mx for mn,mx in zip(minimums,maximums)):
            fail("native min/max activity ordering impossible")
        worst=min(minimums)
        if worst < -1e-8:
            reason=f"native real 80-step {sp} negative concentration: {worst:.10g}"
            if not diagnostic or scheme=="upwind":
                fail(reason)
            scientific_rejections.append(reason)
            print("M1C16_SCIENTIFIC_REJECTED "+reason)
        final=field(case/"0.4"/sp,n**3,"scalar",DIM)
        native_mn=min(final)
        native_mx=max(final)
        if abs(minimums[-1]-native_mn)>2e-5 or abs(maximums[-1]-native_mx)>2e-5:
            fail(f"real final {sp} field differs from logged native extrema")
        maxvals[sp]=max(maximums)
        minvals[sp]=worst
    base=audit(case,n)
    result={"N":n,"scheme":scheme,"native_min_by_species":minvals,
            "native_max_by_species":maxvals,
            "native_activity_budget_max_Bq":base["max_budget_Bq"],
            "emitted_Bq":base["emitted_Bq"],
            "receptors":base["receptors"],
            "real_extrema_samples_per_species":len(moments),
            "scientific_scheme_admissible":not bool(scientific_rejections),
            "scientific_rejection_reasons":scientific_rejections}
    for name in POS:
        for sp in NAMES:
            r=result["receptors"][name]["species"][sp]
            print(f"M1C16_RECEIVER {scheme} N={n} {name} {sp} native={r['interpolated_native_Bq_m3']:.10g} "
                  f"green={r['fixed_analytic_Bq_m3']:.10g} signed_error={r['signed_error_Bq_m3']:.10g} "
                  f"relative={r['relative_error']:.10g}")
    print(f"M1C16_CASE scheme={scheme} N={n} HTmin={minvals['HT']:.10g} "
          f"HTOmin={minvals['HTO']:.10g} maxBq={base['max_budget_Bq']:.10g}")
    return result

def matrix(root,diagnostic=False):
    records={}
    for n in (8,16,24):
        records[n]={scheme:single(root/f"{n}_{scheme}",n,scheme,diagnostic=diagnostic) for scheme in ("upwind","limited")}
    for n in (8,16,24):
        for name in POS:
            for sp in NAMES:
                base=records[n]["upwind"]["receptors"][name]["species"][sp]
                trial=records[n]["limited"]["receptors"][name]["species"][sp]
                print(f"M1C16_COMPARISON N={n} {name} {sp} "
                      f"upwind_rel={base['relative_error']:.12g} limited_rel={trial['relative_error']:.12g} "
                      f"delta_rel={trial['relative_error']-base['relative_error']:.12g}")
    rejection_flags=[]
    for sp in NAMES:
        base=records[24]["upwind"]["receptors"]["R2"]["species"][sp]["relative_error"]
        trial=records[24]["limited"]["receptors"]["R2"]["species"][sp]["relative_error"]
        if trial>1.1*base:
            fail("limitedLinear 1 R2 relative error > frozen 10-percent-relative degradation "+sp)
    for sp in NAMES:
        deltas=[records[n]["limited"]["receptors"]["R2"]["species"][sp]["relative_error"] -
                records[n]["upwind"]["receptors"]["R2"]["species"][sp]["relative_error"]
                for n in (8,16,24)]
        print(f"M1C16_R2_DIAGNOSTIC {sp} relative_error_deltas_8_16_24={deltas} "
              f"improves_all={all(v<0 for v in deltas)}")
    rejected=any(not records[n]["limited"]["scientific_scheme_admissible"] for n in (8,16,24))
    rejected=rejected or bool(rejection_flags)
    if diagnostic and not rejected:
        fail("diagnostic specifically intended to classify the authentic previously observed unsafe limitedLinear 1 run, but no candidate rejected")
    record={"scope":"prescribed-U warped periodic synthetic 3D Gaussian source; cell-centre source unchanged, compare Gauss upwind and Gauss limitedLinear 1 only for HT/HTO",
            "parent_SHA":"320a0a0485a7f1af5c0bc7bdac36a3d06e784176",
            "cases":records,
            "scheme_candidate_scientific_status":"REJECTED" if rejected else "ADMISSIBLE_IN_SYNTHETIC_TESTS",
            "scheme_candidate_additional_rejections":rejection_flags,
            "conclusion":"SCOPED_NATIVE_LIMITED_CONVECTION_DIAGNOSTIC_HIGH_ORDER_REJECTED" if rejected else "SCOPED_NATIVE_LIMITED_CONVECTION_SCHEME_NUMERICAL_COMPARISON"}
    (root/"m1c16-report.json").write_text(json.dumps(record,indent=2,sort_keys=True)+"\n")
    print("M1C16_REPORT_JSON="+json.dumps(record,sort_keys=True,separators=(",",":")))
    if rejected:
        print("M1C16_SCIENTIFIC_STATUS limitedLinear1=REJECTED upwind=ADMISSIBLE")
        print("PASS: M1C-16 diagnostic software properly rejected physically negative limitedLinear 1, six real v2606 3D runs audited; NOT a high-order scientific PASS")
    else:
        print("PASS: M1C-16 six real v2606 limited/upwind source-transport solves, positivity extrema, signed activity budgets and 3D fixed probes")

if __name__=="__main__":
    if len(sys.argv)==5 and sys.argv[1]=="--single":
        single(Path(sys.argv[2]),int(sys.argv[3]),sys.argv[4])
        print("PASS: M1C-16 single genuine 3D native field scheme/extrema check")
    elif len(sys.argv)==5 and sys.argv[1]=="--single-diagnostic":
        single(Path(sys.argv[2]),int(sys.argv[3]),sys.argv[4],diagnostic=True)
        print("PASS: M1C-16 native evidence parsed, candidate scientific admissibility recorded separately")
    elif len(sys.argv)==3 and sys.argv[1]=="--matrix":
        matrix(Path(sys.argv[2]))
    elif len(sys.argv)==3 and sys.argv[1]=="--matrix-diagnostic":
        matrix(Path(sys.argv[2]),diagnostic=True)
    else:
        raise SystemExit("usage: check_convection_scheme.py --single CASE N upwind|limited | --matrix ROOT")
