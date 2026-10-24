#!/usr/bin/env python3
"""Native M1C-21 bounded nonorthogonal correction: independent source/Bq/extrema."""
import math, re, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from check_robustness import records
from check_fixed_receptors import grid, weights
from check_gaussian_green import field, DIM, POS, green, GL32, NUM
from check_convection_scheme import EXTREMA

def fail(m):raise SystemExit("FAIL: M1C21 "+m)

def audit(case,n,count,dt):
    props=(case/"constant/transportProperties").read_text()
    if ("activityDiffusionCorrection faceMules;" not in props or
        "activityConvectionBlend faceMules;" not in props):
        fail("required explicit mode config not present")
    schemes=(case/"system/fvSchemes").read_text()
    for literal in (
        "laplacian(DHT,HTLowDiff) Gauss linear uncorrected;",
        "laplacian(DHTO,HTOLowDiff) Gauss linear uncorrected;",
        "snGrad(HTHighDiff) corrected;",
        "snGrad(HTLowDiff) uncorrected;",
        "snGrad(HTOHighDiff) corrected;",
        "snGrad(HTOLowDiff) uncorrected;",
    ):
        if literal not in schemes:fail("native declared gradient/positive low diffusion scheme missing "+literal)
    s=(case/"solve.log").read_text()
    inv,budget=records(case,count,dt)
    diffusion=re.findall(
       r"CTEAS_M1C21_DIFFUSION species=(HT|HTO) max_abs_high_minus_low_Bq_s=("+NUM+
       r") nonzero_internal_faces=(\d+)",s)
    if len(diffusion)!=2*count:fail("native corrected diffusive flux chronology absent")
    activated={"HT":0,"HTO":0}
    largest={"HT":0.,"HTO":0.}
    for i,(sp,flux,faces) in enumerate(diffusion):
        if sp!=("HT" if i%2==0 else "HTO"):fail("native species diffusion chronology")
        flux=float(flux);faces=int(faces)
        if not math.isfinite(flux) or flux<0 or faces<0:fail("invalid correction face amplitude/count")
        if faces>0 and flux>1e-12:activated[sp]+=1
        largest[sp]=max(largest[sp],flux)
    if min(activated.values())<count//4:
        fail("high-order physical diffusive correction never exercised")
    ext=[tuple(map(float,x)) for x in EXTREMA.findall(s)]
    if len(ext)!=count or any(not all(map(math.isfinite,row)) or row[0]<-1e-8 or row[2]<-1e-8 for row in ext):
        fail("true per-step native HT/HTO min/max negative or absent")
    centres,volumes,lookup=grid(case,n)
    for sp in ("HT","HTO"):
        data=field(case/"0.4"/sp,n**3,"scalar",DIM)
        if not all(map(math.isfinite,data)) or min(data)<-1e-8:
            fail("negative final true-cell species "+sp)
        idx=0 if sp=="HT" else 1
        if abs(sum(a*b for a,b in zip(volumes,data))-inv[-1][idx])>1e-5:
            fail("native real-cell weighted Bq inventory mismatches actual solver")
        for label,xyz in POS.items():
            _,interpolator=weights(xyz,centres,lookup,n)
            expected=green(xyz,GL32)[idx]
            err=abs(interpolator(data)-expected)/expected
            print(f"M1C21_PROBE N={n} dt={dt} {label} {sp} relative_green={err:.12g}")
            if not math.isfinite(err) or err>1.5:
                fail("physical fixed receptor error bound")
    for sp in ("HT","HTO"):
        print(f"M1C21_CORRECTION N={n} dt={dt} {sp} activated_steps={activated[sp]}/{count} max_face_difference_Bq_s={largest[sp]:.12g}")
    print(f"PASS: M1C21 genuine N={n}^3 dt={dt} steps={count} max_signed_budget_Bq={budget:.12g} nonnegative/physical source/Green/local face correction")
if __name__=="__main__":
    if len(sys.argv)!=5:raise SystemExit("usage CASE N STEP_COUNT DT")
    audit(Path(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4]))
