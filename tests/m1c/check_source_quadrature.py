#!/usr/bin/env python3
"""M1C-15 independent native Gaussian source distribution and matched receptor audit."""
import json
import math
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_fixed_receptors import audit
from check_gaussian_green import field,POS
CENTER=(.42,.55,.50)
SIG=.08

def fail(reason):
    raise SystemExit("FAIL: M1C-15 "+reason)

def exact_moments(c):
    m=SIG*math.sqrt(math.pi/2)*(math.erf((1-c)/(SIG*math.sqrt(2)))+math.erf(c/(SIG*math.sqrt(2))))
    e0=math.exp(-(c/SIG)**2/2)
    e1=math.exp(-((1-c)/SIG)**2/2)
    mean=c+SIG**2*(e0-e1)/m
    second=(c*c+SIG**2)+SIG**2*(c*e0-(1+c)*e1)/m
    return mean,second

def run_single(case,n,method):
    if n not in (8,16,24) or method not in ("cellCentre","cellTet4"):
        fail("invalid fixed physical configuration")
    p=(case/"constant/transportProperties").read_text()
    if "sourceShape gaussian;" not in p or "writeSourceProfile true;" not in p:
        fail("source diagnostic configuration not present")
    if method=="cellCentre":
        if "sourceSpatialQuadrature cellTet4;" in p:
            fail("legacy centre mode mistakenly enabled new option")
    elif "sourceSpatialQuadrature cellTet4;" not in p:
        fail("tetra-volume source integration missing from native setup")
    log=(case/"solve.log").read_text()
    if f"sourceSpatialQuadrature={method}" not in log:
        fail("native solver did not report the requested spatial integration")
    profile=field(case/"0/sourceProfile",n**3,"scalar","0 0 0 0 0 0 0")
    centre=field(case/"0/C",n**3,"vector","0 1 0 0 0 0 0")
    volume=field(case/"0/V",n**3,"scalar","0 3 0 0 0 0 0")
    if any(not math.isfinite(x) or x<0 for x in profile):
        fail("source profile negative/nonfinite")
    total=sum(a*v for a,v in zip(profile,volume))
    if abs(total-1)>1e-9:
        fail("source profile native volume-weighted mean differs from one")
    observed_mean=[sum(centre[i][axis]*volume[i]*profile[i] for i in range(n**3))/total for axis in range(3)]
    observed_second=[sum(centre[i][axis]**2*volume[i]*profile[i] for i in range(n**3))/total for axis in range(3)]
    exact=[exact_moments(v) for v in CENTER]
    error_mean=[observed_mean[i]-exact[i][0] for i in range(3)]
    error_second=[observed_second[i]-exact[i][1] for i in range(3)]
    baseline=audit(case,n)
    result={"n":n,"method":method,"native_source_profile_integral":total,
            "native_source_lumped_mean_xyz":observed_mean,
            "native_source_lumped_second_xyz":observed_second,
            "exact_truncated_gaussian_mean_xyz":[v[0] for v in exact],
            "exact_truncated_gaussian_second_xyz":[v[1] for v in exact],
            "lumped_mean_error_xyz":error_mean,"lumped_second_error_xyz":error_second,
            "source_profile_min":min(profile),"source_profile_max":max(profile),
            "receptors":baseline["receptors"],"max_signed_budget_Bq":baseline["max_budget_Bq"],
            "emitted_Bq":baseline["emitted_Bq"]}
    print(f"M1C15_CASE method={method} N={n} sumVp={total:.16g} "
          f"mean_x_err={error_mean[0]:.8g} second_x_err={error_second[0]:.8g} "
          f"R2_HT_rel={baseline['receptors']['R2']['species']['HT']['relative_error']:.9g} "
          f"R2_HTO_rel={baseline['receptors']['R2']['species']['HTO']['relative_error']:.9g}")
    return result

def matrix(root):
    results={}
    for n in (8,16,24):
        results[n]={method:run_single(root/f"{n}_{method}",n,method)
                    for method in ("cellCentre","cellTet4")}
        c=field(root/f"{n}_cellCentre/0/sourceProfile",n**3,"scalar","0 0 0 0 0 0 0")
        q=field(root/f"{n}_cellTet4/0/sourceProfile",n**3,"scalar","0 0 0 0 0 0 0")
        maxdiff=max(abs(x-y) for x,y in zip(c,q))
        if n==8 and maxdiff<=1e-4:
            fail("native tetra-source and center-source spatial profiles did not differ")
        print(f"M1C15_PROFILE_DIFF N={n} max_abs={maxdiff:.10g}")
    for n in (8,16,24):
        for sp in ("HT","HTO"):
            for r in POS:
                centre=results[n]["cellCentre"]["receptors"][r]["species"][sp]
                tet=results[n]["cellTet4"]["receptors"][r]["species"][sp]
                print(f"M1C15_FIXED_DIFF N={n} {r} {sp} "
                      f"centre_rel={centre['relative_error']:.10g} "
                      f"tet_rel={tet['relative_error']:.10g} "
                      f"tet_minus_centre_Bq_m3={tet['interpolated_native_Bq_m3']-centre['interpolated_native_Bq_m3']:.10g}")
    for sp in ("HT","HTO"):
        b=results[24]["cellCentre"]["receptors"]["R2"]["species"][sp]["relative_error"]
        t=results[24]["cellTet4"]["receptors"]["R2"]["species"][sp]["relative_error"]
        print(f"M1C15_R2_24_COMPARE {sp} baseline_rel={b:.12g} cellTet4_rel={t:.12g} "
              f"reduction_in_rel={b-t:.12g} improved={t<b}")
        if t>1.1*b:
            fail("cellTet4 R2 on 24 cubed worse than frozen 10% relative degradation limit "+sp)
    report={"scope":"synthetic warped periodic 3-D Gaussian source cell-volume tetra quadrature versus original centre sample; NOT independent field validation",
            "parent_sha":"b48a7a3d43c6925fdbd316fba2f2c1efb1d2c9f4",
            "cases":results,"conclusion":"SCOPED_OPTIONAL_GAUSSIAN_CELL_VOLUME_QUADRATURE_COMPARISON"}
    (root/"m1c15-report.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("M1C15_REPORT_JSON="+json.dumps(report,sort_keys=True,separators=(",",":")))
    print("PASS: M1C-15 six native Gaussian source-integration comparators, independent moments, real fixed receptors, strict HT/HTO conservation")

if __name__=="__main__":
    if len(sys.argv)==5 and sys.argv[1]=="--single":
        run_single(Path(sys.argv[2]),int(sys.argv[3]),sys.argv[4])
        print("PASS: M1C-15 genuine native Gaussian source-field audit")
    elif len(sys.argv)==3 and sys.argv[1]=="--matrix":
        matrix(Path(sys.argv[2]))
    else:
        raise SystemExit("usage: check_source_quadrature.py --single CASE N cellCentre|cellTet4 | --matrix ROOT")
