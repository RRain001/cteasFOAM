#!/usr/bin/env python3
"""M1C-9: independent native receptor-operator and directional-mesh diagnosis."""
from pathlib import Path
import json
import math
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_convergence import CASES, RECEPTORS, get
from check_gaussian_source import field

DIRECTIONS={
    "space_x":(12,20,6,.02,20),
    "space_y":(8,30,6,.02,20),
    "space_z":(8,20,8,.02,20),
}
ALL=tuple(CASES)+tuple(DIRECTIONS)
CASES.update(DIRECTIONS)
WIDTHS=(.05,.10,.15)
SHIFTS=(-.025,0.,.025)
SPECIES=("HT","HTO")
OLD=("space_lo","baseline","space_hi","time_lo","time_hi")
COMP={"space":("space_lo","baseline","space_hi"),
      "time":("time_lo","baseline","time_hi")}
EXPECTED_FALSE={("space","R2","HTO"),("space","R3","HT"),
                ("space","R3","HTO"),("time","R2","HT")}
REF=Path(__file__).with_name("m1c8_exact_observations.json")


def fail(msg):
    raise SystemExit("FAIL: M1C-9 "+msg)


def probe(case,name):
    # Reuse M1C-8's strict native mesh/field/log, full signed Bq budget
    # and every time-step source-integration audit, untouched from #36.
    previous=get(case,name)
    n=previous["cells"]
    centres=field(case/"0/C",n,"vector","0 1 0 0 0 0 0")
    vols=field(case/"0/V",n,"scalar","0 3 0 0 0 0 0")
    scalar={sp:field(case/"0.4"/sp,n,"scalar","0 -3 -1 0 0 0 0")
            for sp in SPECIES}
    total=sum(vols)
    if not .99999999<=total<=1.00000001:fail("native physical total volume not 1m3")
    source=(.38,.62,.55)
    source_w=[v*math.exp(-sum((a-b)**2 for a,b in zip(c,source))
                             /(2*.13**2))
              for v,c in zip(vols,centres)]
    denom=sum(source_w)
    if denom<=0 or not math.isfinite(denom):fail("source kernel support invalid")
    centroid=[sum(w*c[i] for w,c in zip(source_w,centres))/denom
              for i in range(3)]
    variance=[sum(w*(c[i]-source[i])**2 for w,c in zip(source_w,centres))/denom
              for i in range(3)]
    if any(not math.isfinite(x) or x<=0 for x in variance):
        fail("nonfinite source kernel spatial second moments")
    operators={}
    nearest={}
    for rec,xyz in RECEPTORS.items():
        nearest[rec]=min(math.dist(c,xyz) for c in centres)
        for shift in SHIFTS:
            position=(xyz[0]+shift,xyz[1],xyz[2])
            for width in WIDTHS if shift==0 else (.10,):
                weights=[v*math.exp(-sum((a-b)**2 for a,b in zip(c,position))
                                         /(2*width**2))
                         for v,c in zip(vols,centres)]
                normalization=sum(weights)
                if normalization<=0 or not math.isfinite(normalization):
                    fail(f"{name} {rec} receptor support 0")
                for sp in SPECIES:
                    conc=sum(w*s for w,s in zip(weights,scalar[sp]))/normalization
                    if not math.isfinite(conc) or conc<0:
                        fail(f"{name} {rec} {sp} invalid Bq/m3")
                    key=f"{rec}_{sp}_w{width:.2f}_dx{shift:+.3f}"
                    operators[key]=conc
    # Orthogonal check of the nominal M1C-8 averaging implementation.
    for rec in RECEPTORS:
        for sp in SPECIES:
            nominal=operators[f"{rec}_{sp}_w0.10_dx+0.000"]
            old=previous[f"{rec}_{sp}_Bq_m3"]
            if abs(old-nominal)>1e-7:
                fail(f"{name} independent nominal observation operator != M1C-8")
    return {
        "native_qa":previous,
        "centroid_m":centroid,
        "source_second_moment_m2":variance,
        "nearest_receiver_m":nearest,
        "cell_spacing_m":{k:1/v for k,v in
                          (("x",previous["nx"]),("y",previous["ny"]),
                           ("z",previous["nz"]))},
        "receptor_Bq_m3":operators,
    }


def assess(root):
    ref=json.loads(REF.read_text())
    if ref["source_sha"]!="68370fb4483608f96ceb7073f4fe8ea69ec0dc59":
        fail("wrong independently frozen prior experimental SHA")
    if set(ref["receptors"])!=set(OLD):
        fail("prior M1C-8 cases incomplete")
    native={name:probe(root/name,name) for name in ALL}
    # The published prior native test is a separate, SHA-verified source,
    # not self-generated values from these M1C-9 cases.
    max_old=0
    for name in OLD:
        for key,val in ref["receptors"][name].items():
            rec,sp=key.split("_")[:2]
            actual=native[name]["receptor_Bq_m3"][f"{rec}_{sp}_w0.10_dx+0.000"]
            max_old=max(max_old,abs(actual-val))
            if abs(actual-val)>1e-5:
                fail(f"prior exact-SHA M1C-8 result drift {name}:{key}: {actual-val}")
    print(f"PASS: M1C-9 exact independent M1C-8 prior reference 30 native concentrations max_abs_diff={max_old:.12g} Bq/m3")
    widths={}
    shifts={}
    def delta(axis,rec,sp,width,shift):
        low,mid,hi=COMP[axis]
        key=f"{rec}_{sp}_w{width:.2f}_dx{shift:+.3f}"
        a,b,c=[native[name]["receptor_Bq_m3"][key] for name in (low,mid,hi)]
        d1=abs(b-a);d2=abs(c-b)
        return {"low":a,"mid":b,"high":c,"delta_low_mid":d1,
                "delta_mid_high":d2,"ratio":d2/max(d1,1e-9),
                "improved":bool(d2<d1 and d1>1e-9)}
    for width in WIDTHS:
        label=f"{width:.2f}"
        widths[label]={axis:{
            f"{rec}_{sp}":delta(axis,rec,sp,width,0)
            for rec in RECEPTORS for sp in SPECIES}
            for axis in COMP}
    for shift in SHIFTS:
        label=f"{shift:+.3f}"
        shifts[label]={axis:{
            f"{rec}_{sp}":delta(axis,rec,sp,.10,shift)
            for rec in RECEPTORS for sp in SPECIES}
            for axis in COMP}
    original=widths["0.10"]
    previous_false=set()
    for axis,checks in original.items():
        for key,detail in checks.items():
            if not detail["improved"]:
                rec,sp=key.split("_")
                previous_false.add((axis,rec,sp))
            if ref["flags"][axis][key+"_Bq_m3"]!=detail["improved"]:
                fail(f"prior M1C-8 monotonicity flag drift {axis}/{key}")
    if previous_false!=EXPECTED_FALSE:
        fail(f"original four non-monotonic observables were suppressed: {previous_false}")
    directions={}
    baseline=native["baseline"]
    for direction in DIRECTIONS:
        case=native[direction]
        diffs={}
        for rec in ("R2","R3"):
            for sp in SPECIES:
                key=f"{rec}_{sp}_w0.10_dx+0.000"
                a=baseline["receptor_Bq_m3"][key]
                b=case["receptor_Bq_m3"][key]
                diffs[f"{rec}_{sp}"]=b-a
        directions[direction]={
            "delta_vs_baseline_Bq_m3":diffs,
            "kernel_centroid_delta_m":[a-b for a,b in
                                       zip(case["centroid_m"],
                                           baseline["centroid_m"])],
            "kernel_moment_delta_m2":[a-b for a,b in
                                     zip(case["source_second_moment_m2"],
                                         baseline["source_second_moment_m2"])],
            "receptor_nearest_cell_m":case["nearest_receiver_m"],
        }
    result={
        "parent_sha":"68370fb4483608f96ceb7073f4fe8ea69ec0dc59",
        "scope":"synthetic laminar PISO observation-operator and directional discretization sensitivity only",
        "original_four_failures":[{"axis":a,"receptor":r,"species":s} for a,r,s in sorted(previous_false)],
        "fixed_width_refinements":widths,
        "fixed_shift_refinements":shifts,
        "directional_refinements":directions,
        "per_native_case":native,
        "previous_native_field_consistency_max_Bq_m3":max_old,
        "conclusion":"PARTIAL_EXPLANATION; NO GLOBAL CONVERGENCE CLAIM",
    }
    (root/"m1c9-diagnostic-report.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    for axis in COMP:
        for width in WIDTHS:
            bad=[k for k,v in widths[f"{width:.2f}"][axis].items()
                 if not v["improved"]]
            print(f"M1C9_OPERATOR axis={axis} width={width:.2f} non_improving={bad}")
        for shift in SHIFTS:
            bad=[k for k,v in shifts[f"{shift:+.3f}"][axis].items()
                 if not v["improved"]]
            print(f"M1C9_POSITION axis={axis} dx={shift:+.3f} non_improving={bad}")
    for k,d in directions.items():
        print(f"M1C9_DIRECTION {k} deltaBq_m3={d['delta_vs_baseline_Bq_m3']}")
    print("M1C9_REPORT_JSON="+json.dumps(result,separators=(",",":"),sort_keys=True))
    print("PASS: M1C-9 eight native PISO cases, previous frozen failures retained, receptor kernel and directional mesh diagnosis")
    return result


if __name__=="__main__":
    if len(sys.argv)==3 and sys.argv[1]=="--matrix":
        assess(Path(sys.argv[2]))
    elif len(sys.argv)==4 and sys.argv[1]=="--single":
        name=sys.argv[3]
        if name not in DIRECTIONS:fail("single mode requires directional test")
        probe(Path(sys.argv[2]),name)
        print(f"PASS: M1C-9 directional {name} genuine native signed Bq and boundary checks")
    else:
        raise SystemExit("usage: check_receiver_diagnostics.py --matrix ROOT | --single CASE space_x|space_y|space_z")
