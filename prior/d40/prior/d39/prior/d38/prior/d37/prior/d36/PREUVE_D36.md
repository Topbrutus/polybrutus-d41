# Polybrutus D36 — Constructive exact proof and independent controls

**Immutable parent:** `Topbrutus/polybrutus-d35@v0.1.0` at scientific commit `72989a1fd2a7fb6e28522a0f2dcbe168924ca78d`, DOI `10.5281/zenodo.24412435`, original certificate SHA-256 `7390c10f46132faeca0ccd8b52050d5ed899fbd4a078b9dd02407d922a4e4f22`.

**Conservation:** Extend the 40 frozen norm-one points by a zero coordinate; all 780 pairwise squared distances remain exactly unchanged. Append one perpendicular unit point giving 40 new distances squared `2R²`, for `C(41,2)=820` labelled pairs.

**Positive semidefiniteness and minimal dimension:** Construct exact rational Gram `G` with rank-three canonical seven-point block, rank-four centered simplex five-point block, and 29 mutually independent orthogonal unit axes. The Gram is PSD by construction with unit diagonal and rank `36`. The centered squared-distance Gram `B=-(1/2)HDH=HGH` has exact rank **36**, proving minimum Euclidean affine realization dimension 36 for this *specific* full labelled metric.

**Distance arithmetic:** Normalized classes `k=4d²/R²` are `{4:6,7:2,8:789,9:3,10:10,12:6,13:2,16:2}`, totaling 820. Their squared distances sum to `6607/4` at R=1.

**Selected ideal graph:** For the deliberately chosen complete graph `w_ij=1/d_ij²` and unit scalar masses, exact matrix ranks give eigenvalue `41/2` of multiplicity **30** and `20` of multiplicity **4**. These are dimensionless model values rather than experimentally measured resonances.

**Two distinct implementations:** Primary `src/verify_d36.py` constructs all exact rational Gram entries, 820 labelled squared distances, checks all 780 parent labels and both Euclidean and selected-model ranks. Its four adversarial negative controls reject a corrupted distance, duplicated point, nonorthogonal tilt, and histogram-preserving label permutation. Independently, `tests/independent_d36.py` derives each labelled distance using Python `fractions.Fraction` rules, **without importing the primary Gram generator**, and checks the parent checksum, all 820 labels, rank, and model eigenspaces. Both output the identical SHA-256 for `results/certificate_d36.json`.

**Limitations:** This uses established orthogonal Euclidean embedding and uniformly joined graph techniques to describe one exact synthetic metric, not a novel universal theorem, additional physical spacetime dimensions, measured frequencies, particles, or independently human-peer-reviewed research.
