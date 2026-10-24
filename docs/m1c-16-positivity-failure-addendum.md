# M1C-16 failed positivity candidate — diagnostic protocol addendum

This addendum is recorded after the FIRST authentic v2606 8³ `Gauss limitedLinear 1` scientific FAIL and **before adapting the diagnostic runner**. Initial official v2606 optional-extrema build failed due to `Foam::messageStream` not exposing a `precision()` method; this was corrected using local `std::ostringstream`, with no transport/source equations or numerical thresholds changed. The next genuine OpenFOAM run produced `CTEAS_M1C extrema_HT_min_Bq_m3=-2.523170723` at an intermediate 8³ step; exact native last step had positive HT extrema. Thus **checking only final output would miss an unacceptable transient negative**. Existing `Gauss upwind` 8³ all-step minima are nonnegative.

The original acceptance stays binding: any `limitedLinear 1` native cell less than −1e−8 Bq/m³ is a **scientific REJECTED/FAIL**, not an alternative solver PASS. The new diagnostic runner should finish real 8³/16³/24³ upwind/limited pairs and report actual fixed-point R1/R2/R3 errors and signed inventories even for rejected candidates, but must:
- preserve a separate **strict** checker invocation that exits nonzero on real physical negativity and verify it rejects with expected specific reason;
- record `scientific_scheme_admissible=false` per offending case, show min HT/HTO and all timestep native extrema counts, never claim the limited solution is acceptable or safe;
- require upwind cases to satisfy strict original frozen positivity and conservation;
- verify even rejected limited runs have real native 80-step fields, exact source window and Bq conservation, but do not call their physical performance PASS;
- if any limited 24³ R2 relative error exceeds frozen 10%-relative degradation limit, retain it as a separate scientific failure, not a software crash;
- label the overall native script exit-0 as **diagnostic software PASS, limitedLinear scientific FAIL** *only if the unsupported candidate is explicitly classified and reported*; do not call the high-order scheme a successful physical validation.
- do not introduce clipping, normalize negative field values, relax thresholds, or switch upstream default `Gauss upwind`. Formal external scientific state remains NOT_RUN.

The initial failure log is `/home/lhy/dev/cteasFOAM/reports/m1c16/dev2.log`. Numerical comparisons on any rejected scheme describe a *diagnostic counterfactual*, not an admissible numerical prediction.
