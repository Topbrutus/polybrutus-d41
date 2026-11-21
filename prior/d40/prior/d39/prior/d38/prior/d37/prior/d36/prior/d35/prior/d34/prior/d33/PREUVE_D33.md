# Polybrutus D33 — Exact proof and independent adversarial tests

**Frozen scientific parent:** `Topbrutus/polybrutus-d32@v0.1.0`; scientific commit `cb63ed7ab447b13e251d92571fe0f1062db31cca`; DOI `10.5281/zenodo.24411375`; original certificate SHA-256 `66b7605ab168d99ee33014ed1d61b5fbb3915c32041bdf44dd529c0bb56bb016`.

**Metric conservation.** Zero-coordinate inclusion of 37 labelled unit points preserves the 666 old squared pair distances exactly. The new perpendicular norm-one point adds 37 squared distances all equal to `2R²`, totaling **703 = C(38,2)** unordered labelled pairs.

**Euclidean realizability and minimum dimension.** The exact rational PSD Gram `G` has a rank-three canonical block, rank-four centered simplex block and 26 mutually orthogonal unit axes, and unit diagonal. Thus `rank(G)=33`. Its complete squared-distance matrix uniquely determines centered Gram `B=-(1/2)HDH=HGH`. An antipodal pair places the origin in the affine hull and exact rational rank gives `rank(B)=33`; the minimum possible affine Euclidean realization dimension of this specific complete labelled metric is therefore 33.

**Arithmetic invariants.** The normalized histogram `k=4d²/R²` is `{4:6,7:2,8:672,9:3,10:10,12:6,13:2,16:2}`, for a total of 703 pairs. The squared distances sum to `5671/4` at unit radius.

**Ideal chosen model graph.** Assign unit scalar masses and stipulated complete weights `w_ij=1/d_ij²`. The new vertex joins the previous 37 at weight 1/2, giving two exact-rational eigenmode multiplicities: eigenvalue **19 multiplicity 27** and eigenvalue **37/2 multiplicity 4**. These are dimensionless mathematical graph properties, not measured physical resonances or evidence of physical extra dimensions.

**Two independent programs.** `src/verify_d33.py` builds the rational Gram and 703 exact labels and checks 666 frozen parent values, centered Gram rank and Laplacian eigenspaces. `tests/independent_d33.py` uses Python `fractions.Fraction` structural distance rules without importing the first Gram generator; it independently checks every old and new label, parent certificate integrity, centered rank and graph nullities. Both generate identical deterministic certificate SHA-256. The primary checker rejects four intentionally invalid configurations: a corrupted distance, a duplicated vertex, a nonorthogonal tilt and an unlabeled-histogram-preserving label swap; the independent checker rejects corruption and relabelling.

**Scientific limits.** Orthogonal Euclidean extensions and uniform graph joins are known methods. The reproducible case study does not claim a new universal theorem, observed particles, measured frequencies, 33 physical spacetime dimensions, or independent human peer review.
