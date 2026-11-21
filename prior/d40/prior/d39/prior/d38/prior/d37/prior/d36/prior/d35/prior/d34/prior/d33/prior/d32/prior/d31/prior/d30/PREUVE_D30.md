# Polybrutus 30D — Exact constructive proof and independent checks

**Immutable source:** `Topbrutus/polybrutus-d29@v0.1.0`, commit `cbe179697cc1a287fc0637912ea6a6fd890eeca4`, DOI `10.5281/zenodo.24410725`, original certificate SHA-256 `d718f784bf20dccc1b8c08821466b352316054900fba20705bc28d79763a4731`.

**1 — Labelled conservation.** Embed the original 34 radius-one vertices in a new coordinate space by appending zero. Adding a perpendicular norm-one vector leaves all **561 old labelled squared distances unchanged** and contributes 34 squared distances equal to `2R²`, for **C(35,2)=595**.

**2 — Euclidean PSD and minimum affine dimension.** An explicit rational positive-semidefinite Gram matrix `G` has a rank-three canonical block, a rank-four simplex block, and 23 independent orthogonal basis axes; `rank(G)=30`, with unit diagonal. Its complete labelled squared-distance matrix `D` determines `B=-(1/2)HDH=HGH`. Two core points are antipodal, placing the origin in their affine hull. Exact rank `rank(B)=30` establishes the minimum affine Euclidean realization dimension **30** for this one metric.

**3 — Arithmetic counts.** The normalized squared-distance histogram `k=4d²/R²` is `{4:6,7:2,8:564,9:3,10:10,12:6,13:2,16:2}`, totaling 595. The sum of squared pair distances at radius one is `4807/4`.

**4 — Chosen ideal graph model.** For stipulated unit masses and complete edge weights `w_ij=1/D_ij`, the new vertex joins all 34 old vertices at weight 1/2. Exact rational ranks give eigenvalue `35/2` of multiplicity **24** and `17` of multiplicity **4**. They are abstract dimensionless graph modes, not measured oscillations.

**5 — Two separately implemented verifiers.** `src/verify_d30.py` constructs Gram scalar products and all 595 labelled pair values, checks all 561 immutable parent labels, centered rank and Laplacian nullities, and detects four intentionally false configurations: corrupted distance, duplicate point, non-orthogonal tilt, and histogram-preserving relabelling. `tests/independent_d30.py` independently reconstructs the labelled distance metric using `fractions.Fraction` rules, without importing the primary Gram generator, and checks all pairs and ranks. Both print the same deterministic SHA-256 of `results/certificate_d30.json`.

**Limits:** Orthogonal Euclidean extension and uniform graph joins are established mathematical techniques. This is a reproducible specific case study, not a novel universal theorem or experimental evidence of additional physical spacetime dimensions, particles or measured resonance. No independent human peer review has occurred.
