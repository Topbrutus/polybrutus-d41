# Polybrutus D24 — Exact mathematical proof and independent counter-tests

**1. Orthogonal conservation.** The 28 immutable radius-one parent vertices are embedded with a zero appended storage coordinate. A fresh perpendicular unit vertex creates 28 squared distances equal to 2; all 378 previously labelled squared distances remain identical to their frozen D23 certificate.

**2. Euclidean validity.** The 29-point rational Gram matrix separates the 3-dimensional original core, 4-dimensional centered simplex, and 17 new orthogonal directions. Exact rank `24` and factorization `B=-(1/2)HDH=HGH` prove that the Euclidean metric is realizable. Because two core points are antipodal, the origin lies in their affine hull, and the minimum affine Euclidean realization dimension is **24**.

**3. Exact combinatorics.** `C(29,2)=406`; class histogram `{4:6,7:2,8:375,9:3,10:10,12:6,13:2,16:2}`; total squared-distance sum `3295/4` for unit radius. The 28 additional squared distances are 2.

**4. Weighted graph model.** Each edge is assigned the chosen hypothetical weight `w_ij=1/d_ij²`. Matrix rank computations in two programs establish exactly `nullity(L24-(29/2)I)=18` and `nullity(L24-14I)=4`. These are not observations of a physical oscillator.

**5. Two separately implemented verifiers.** `src/verify_d24.py` constructs Gram coefficients exactly and checks all labelled pairs; `tests/independent_d24.py` reconstructs them separately using `fractions.Fraction` and never imports the primary generator. Both compare the full 378-parent-pair labelled map, compute the Gram rank and graph multiplicities, and print the same SHA-256 certificate hash. Negative controls reject a corrupted pair, duplicate point, tilted point, and point-label permutation hidden by the unlabeled distance histogram.

**6. Limitations.** Orthogonal completion and uniform graph joins are established mathematical procedures. This is a specific synthetic case study, not a novel universal theorem. Nothing here demonstrates an extra physical spacetime dimension, a measured resonance, or human peer review.
