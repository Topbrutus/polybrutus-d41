# Polybrutus D34 — Exact constructive proof and independent tests

**Frozen source:** `Topbrutus/polybrutus-d33@v0.1.0`, scientific commit `bc78e6c6fe1755a069d8860a0144cbebbb1168a6`, DOI `10.5281/zenodo.24411569`, original certificate SHA-256 `86516b8a6393bf9fdf4e698d6e439d4db783a71956730ac35bbc89e2c8e7b734`.

**1. Labelled conservation.** Appending a zero coordinate to each of the 38 parent norm-one points leaves every one of their 703 labelled squared pair distances unchanged. A new perpendicular unit point contributes 38 additional squared distances `2R²`. Exactly `C(39,2)=741` unordered labelled pairs result.

**2. Positive semidefiniteness and minimal dimension.** The explicit exact rational Gram `G` has three canonical independent directions, four centered simplex directions and 27 additional orthogonal axes. Its diagonal equals one and `rank(G)=34`. The unique centered Euclidean Gram `B=-(1/2)HDH=HGH` is PSD, with exact rank 34. Because two original canonical vertices are antipodal, the origin is in the affine hull. No Euclidean realization of the same complete labelled squared-distance metric can have affine dimension below 34.

**3. Exhaustive invariants.** The class histogram `k=4d²/R²` is `{4:6,7:2,8:710,9:3,10:10,12:6,13:2,16:2}`; it sums to 741. For `R=1`, squared pair distances sum to `5975/4`. All 703 frozen predecessor labels are independently checked, not only the histogram.

**4. Chosen dimensionless ideal graph.** With scalar masses one and stipulated weights `w_ij=1/d_ij²`, the graph Laplacian gains a uniform join of coupling 1/2. Two independent exact-rational matrix-rank calculations verify eigenvalue `39/2` of multiplicity 28 and `19` of multiplicity 4. No physical frequency is measured or inferred.

**5. Primary and independent verification.** `src/verify_d34.py` constructs a rational Gram and all 741 pair labels, confirms all 703 conserved parent labels, checks Euclidean ranks and the chosen graph eigenspaces. Four adversarial controls reject a corrupted squared distance, a duplicate point, a non-orthogonal tilt and a relabelling with the same histogram. Independently, `tests/independent_d34.py` reconstructs every distance with `fractions.Fraction` rules without importing the primary Gram generator; it checks the SHA-256 of the frozen parent certificate, parent conservation, exact ranks and selected model eigenmodes. Both scripts print the same deterministic SHA-256 of `results/certificate_d34.json`.

**Limitations:** Orthogonal Euclidean extensions and uniform graph joins are established mathematical techniques. This specific exact example is not a novel universal theorem, a measured resonance, a particle discovery, evidence of 34 physical spacetime dimensions, or independently human peer-reviewed work.
