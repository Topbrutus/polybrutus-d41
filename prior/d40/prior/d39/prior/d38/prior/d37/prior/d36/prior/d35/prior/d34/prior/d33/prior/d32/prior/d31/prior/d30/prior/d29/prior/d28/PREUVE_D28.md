# Polybrutus 28D - Exact proof and independently checked negative controls

**Frozen parent:** `Topbrutus/polybrutus-d27@v0.1.0`, source commit `a311d63291d9254664340b434c5a13162d0658c1`, DOI `10.5281/zenodo.24410593`, certificate SHA-256 `e464dab0dc5d6586a94a714f1fee06b93710b3ba2763d710501970e6a3e399c4`.

**Conservation.** Embedding 32 frozen radius-one vectors by appending a zero coordinate does not change the existing 496 labelled squared distances. A new perpendicular unit vector adds 32 squared distances of exactly `2R²`. Thus `C(33,2)=528`.

**Metric realization.** An explicit rational PSD Gram matrix has three independent canonical directions, four regular-simplex directions, and 21 further orthogonal unit-axis directions, for exact rank 28. Its diagonal entries equal one. `D_ij=2-2G_ij` and the centered Euclidean squared-distance Gram `B=-(1/2)HDH=HGH` has exact rank 28. An antipodal canonical pair places the origin in the affine hull, proving minimum Euclidean affine dimension 28.

**Arithmetic invariants.** Normalized index `k=4d²/R²` yields `{4:6,7:2,8:497,9:3,10:10,12:6,13:2,16:2}`, total 528 pairs. Sum of squared distances at unit radius is `4271/4`.

**Synthetic spectral model.** In the explicitly chosen ideal complete weighted graph `w_ij=1/D_ij`, the new point uniformly joins predecessors at weight 1/2. Exact nullities `nullity(L28-(33/2)I)=22` and `nullity(L28-16I)=4` are independently verified. Neither eigenvalue represents a measured wave frequency.

**Independent checks.** The primary `src/verify_d28.py` constructs exact rational Gram entries. The separately implemented `tests/independent_d28.py` calculates all labelled distances using Python `fractions.Fraction` structural rules without importing the Gram generator. Both check every one of the 496 old labels, all 528 current labels, Gram ranks, and graph nullities; both print the same SHA-256. Primary adversarial tests cover a corrupted distance, a duplicate vertex, a nonperpendicular tilt, and a histogram-preserving label swap.

**Limitations.** This is a particular reproducible mathematical configuration using known techniques, not a new universal theorem. It does not demonstrate additional physical spacetime dimensions, observed resonance, a particle, or independent human peer review.
