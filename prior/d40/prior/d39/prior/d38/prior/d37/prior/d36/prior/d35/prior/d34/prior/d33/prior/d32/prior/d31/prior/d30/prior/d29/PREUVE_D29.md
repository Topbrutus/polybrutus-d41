# Polybrutus D29 — Exact proof and independent adversarial checks

**Source of truth:** Frozen Polybrutus D28 science tag `v0.1.0`, source commit `496e6373abbabfad5c2e99adfcb3994ce8e39531`, DOI `10.5281/zenodo.24410652`, parent certificate SHA-256 `d21c21645bc96aa480c840fa24e12e755f6a21d8df8a1d51dac989c7caed1788`.

**1. Labelled distance conservation.** Zero-coordinate inclusion of 33 old unit vectors preserves all **528/528** prior labelled squared distances exactly. The new perpendicular norm-one point contributes 33 squared distances equal to `2R²` and there are `C(34,2)=561` total pairs.

**2. Positive semidefinite Gram construction.** The exact rational `G` consists of 3 canonical directions, 4 simplex directions, and 22 independent orthogonal unit-axis directions. It has unit diagonal, is positive semidefinite, and has rank 29. The exact centered squared-distance Gram `B=-(1/2)HDH=HGH` is PSD of rank 29. The origin is in the affine hull of an antipodal canonical pair. Since complete labelled squared distances uniquely determine `B`, no Euclidean realization of those distances can have affine dimension less than 29.

**3. Exact invariants.** The 561 labelled squared distances comprise the normalized histogram `{4:6,7:2,8:530,9:3,10:10,12:6,13:2,16:2}`, with sum `4535/4` for unit radius.

**4. Chosen ideal graph.** Deliberately assign complete weights `w_ij=1/d_ij²` and unit masses. The new vertex is a uniform join with weight 1/2, and exact rank calculations independently verify `nullity(L29-17I)=23`, `nullity(L29-(33/2)I)=4`. These are abstract, dimensionless model eigenvalues, not experimental frequencies.

**5. Two separately implemented verifiers.** `src/verify_d29.py` constructs exact Gram scalar products and the complete pair map; `tests/independent_d29.py` uses separate `fractions.Fraction` distance rules and never imports the primary Gram generator. Both check all 561 pair labels, all 528 original parent labels, centered Gram rank and Laplacian nullities, and output the identical SHA-256. The primary checker rejects a damaged distance, a duplicated point, nonorthogonal tilt, and a histogram-preserving label permutation; the secondary independently rejects corruption and relabelling.

**Limits.** Standard orthogonal extension and uniform weighted graph joins are applied to one reproducible mathematical case study, not a new universal theorem, measured physical phenomenon, physical spacetime dimensionality or human peer-reviewed proof.
