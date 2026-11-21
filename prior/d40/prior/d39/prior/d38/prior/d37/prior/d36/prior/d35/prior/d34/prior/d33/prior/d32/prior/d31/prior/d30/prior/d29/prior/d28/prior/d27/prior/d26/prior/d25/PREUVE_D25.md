# Polybrutus D25 — Exact construction and independent controls

**Frozen parent:** `Topbrutus/polybrutus-d24@v0.1.0` (commit `16d9f0ee705d15c220bc02e18563a5e5c5c1e5fc`, DOI `10.5281/zenodo.24410102`). The D24 exact certificate and source programs are archived unchanged under `prior/d24/` and the parent certificate is checked by SHA-256.

**1. Parent preservation.** The new perpendicular unit point changes no scalar product or squared distance between the 29 old points; exactly `406/406` labelled predecessor pairs remain invariant. All 29 new squared distances equal `2R²`.

**2. Metric and affine dimension.** An exact rational positive-semidefinite Gram matrix `G` has rank 25 and unit diagonal, with `D_ij=2-2G_ij`. The centered squared-distance Gram matrix `B=-(1/2)HDH=HGH` is positive semidefinite of rank **25**. Two antipodal core points place the origin in the affine span, proving the minimum affine Euclidean dimension is exactly 25.

**3. Arithmetic checks.** For `R=1` there are 435 distinct vertex pairs. The histogram is `{4:6,7:2,8:404,9:3,10:10,12:6,13:2,16:2}` for `4d²/R²`, and total squared-distance sum is `3527/4`. Every label and distance is represented in `results/certificate_d25.json`.

**4. Ideal weighted graph.** With explicitly chosen `w_ij=1/d_ij²`, the weighted Laplacian is the uniform one-half join of the D24 model to a new vertex. Its exact eigenvalue nullities are `nullity(L25-15I)=19` and `nullity(L25-(29/2)I)=4`, checked by both verifiers. These are synthetic, dimensionless, model-dependent eigenvalues.

**5. Independent reproduction.** `src/verify_d25.py` derives squared distances from an exact Gram matrix; `tests/independent_d25.py` uses independently written `fractions.Fraction` distance rules and rechecks the parent certificate. Both compare 435 labelled pairs, preserve 406 predecessor pairs and calculate the Laplacian rank. The primary checker rejects four deliberately malformed configurations, while the secondary independently rejects corrupted values and a label permutation invisible to histograms.

**Limits:** Orthogonal extension and uniform graph joins are established techniques. This reproducible metric is not a novel universal theorem, an experiment, evidence for a physical spacetime dimension, a measured wave frequency or a peer-reviewed result.
