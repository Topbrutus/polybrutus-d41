# Polybrutus D32 — Exact metric proof and independent negative controls

**Frozen predecessor:** `Topbrutus/polybrutus-d31@v0.1.0`; scientific commit `0c1f546a6c5442bc83a5bdde506132cce20b6e05`; DOI `10.5281/zenodo.24411160`; certificate SHA-256 `7c373c50d61ca2964eb832fc6fb8178590715e603f5fdada122575f624f60e8c`.

**Conservation.** The zero-coordinate embedding of 36 frozen norm-one vectors preserves all 630 labelled predecessor squared-distance values. Appending one orthogonal unit point adds 36 squared distances equal to `2R²`. Therefore the exact 37-point squared-distance metric has 666 labelled pairs.

**Euclidean realizability and minimal dimension.** The explicit rational Gram `G` is positive semidefinite with unit diagonal, comprising the canonical rank-three block, rank-four centered simplex and 25 independent unit-axis directions. Its rank is 32. The centered squared-distance Gram `B=-(1/2)HDH=HGH` has exact rank 32; the antipodal canonical pair places the origin in the affine hull. Because the complete labelled squared-distance matrix uniquely determines `B`, every Euclidean realization of this metric has affine rank 32.

**Arithmetic checks.** The normalized distance-family histogram is `{4:6,7:2,8:635,9:3,10:10,12:6,13:2,16:2}`, summing to 666 pairs. Their squared distances sum to `5375/4` at `R=1`.

**Specified ideal graph.** Under the freely selected mathematical coupling `w_ij=1/d_ij²` with scalar masses one, independently computed exact matrix nullities prove eigenvalue `37/2` of multiplicity 26 and eigenvalue `18` of multiplicity 4. No physical frequency or new physical dimension follows.

**Independent reproducibility.** Primary `src/verify_d32.py` builds the rational Gram and all 666 labels; independent `tests/independent_d32.py` uses Python `fractions.Fraction` metric rules without importing the Gram generator. Both check each of 630 frozen predecessor labels and exact Euclidean/spectral ranks, and both output the same SHA-256 for `results/certificate_d32.json`. Primary negative tests reject a corrupted distance, a duplicate vertex, a nonperpendicular tilt and histogram-preserving relabelling. Secondary adversarial tests reject corruption and relabelling.

**Scientific limits.** Orthogonal Euclidean extension and uniformly joined weighted graph Laplacians are known techniques. This is one exact reproducible synthetic case study, not a new universal theorem, observed physical spacetime dimension, measured resonance or independently peer-reviewed human research.
