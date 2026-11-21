# Polybrutus D41 — Constructive exact proof and controls

**Frozen parent:** `Topbrutus/polybrutus-d40@v0.1.0`, scientific commit `2a9fa8265b709d6d971afb6c63b7b024011e67a1`, DOI `10.5281/zenodo.24414238`, parent certificate SHA-256 `af361ae471397b6d91beed87fe65c268a7139c02fbc0187400152fb8b202ffca`.

**Exact label conservation.** The zero-coordinate extension of 45 prior unit vectors preserves each of 990 predecessor squared-distance labels. One new perpendicular unit point contributes 45 new squared distances `2R²`, yielding **1035=C(46,2)** unordered labelled pairs.

**Constructive Euclidean rank proof.** Exact rational Gram `G` has canonical rank-three, centered regular-simplex rank-four, and 34 independent perpendicular unit-axis blocks. It is PSD by explicit Gram construction with diagonal one and rank 41. The uniquely determined centered Gram `B=-(1/2)HDH=HGH` is PSD with exact rank 41. The inherited antipodal canonical pair ensures the origin belongs to the affine hull. Therefore no Euclidean embedding of this complete labelled metric into fewer than 41 affine dimensions exists.

**Exhaustive arithmetic.** The normalized histogram `k=4d²/R²` is `{4:6,7:2,8:1004,9:3,10:10,12:6,13:2,16:2}`, totaling 1035 pairs. The squared distances sum to `8327/4` for R=1. All 990 predecessor values are checked with their original pair labels, beyond just matching an unlabeled histogram.

**Separately selected model.** Scalar unit masses and ideal complete graph edge weights `w_ij=1/d_ij²` yield exact Laplacian eigenvalue 23 of multiplicity 35, and 45/2 of multiplicity 4. These are mathematical dimensionless modes, not measured hertz.

**Two different verification implementations.** The primary `src/verify_d41.py` generates exact rational Gram entries and 1035 squared-distance pairs, checks 990/990 frozen predecessor labels against the checksum-verified parent, and verifies centered Euclidean and ideal graph ranks. Four negative controls reject corrupted pair distance, duplicated point, nonorthogonal tilt, and histogram-preserving relabelling. The separate `tests/independent_d41.py` derives all pair distances with Python `fractions.Fraction` structural rules, *without importing the primary Gram generator*, and confirms the parent SHA, labels and ranks. Both programs output the same deterministic certificate SHA-256.

**Scientific limits.** Standard orthogonal Euclidean extension and uniform graph joining give a particular synthetic construction. This does not establish observed physical extra dimensions, frequencies, particles, a new universal theorem, or independent human peer review.
