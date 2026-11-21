# Polybrutus 35D — Constructive proof and independent controls

**Immutable parent:** `Topbrutus/polybrutus-d34@v0.1.0`, commit `21385a5875d472af8f035372096a46f5501cee1a`, DOI `10.5281/zenodo.24412288`, exact certificate SHA-256 `95274561481b735515a43dd185282cbe6b2424611f06efa5c6e5d100d433cfb6`.

**1. Labelled conservation.** Append a zero coordinate to each predecessor point, preserving all **741** old labelled squared distances. Add a new perpendicular norm-one point, contributing **39** squared distances equal to `2R²`. The total is **780=C(40,2)**. The verifiers compare every old labelled pair, not merely the distance histogram.

**2. Euclidean realizability and minimal dimension.** Construct an explicit positive-semidefinite rational Gram `G` with a rank-three canonical seven-point block, a rank-four centered five-point simplex block and 28 independent unit-axis directions. It has unit diagonal and rank 35. The unique centered distance Gram `B=-(1/2)HDH=HGH` is PSD and has exact rank **35**, which is the minimum Euclidean affine realization dimension of this complete labelled metric. An antipodal canonical pair places the origin in the affine hull.

**3. Exact pair counts.** The normalized squared-distance classes `k=4d²/R²` are `{4:6,7:2,8:749,9:3,10:10,12:6,13:2,16:2}`, totaling 780. The sum of squared pair distances at unit radius is **6287/4**.

**4. Separate graph assumption and eigenspaces.** The **stipulated** equal-mass complete graph with edge weights `w_ij=1/d_ij²` has Laplacian eigenvalue 20 of exact multiplicity **29** and 39/2 of multiplicity **4**, as checked by exact-rational matrix nullities in both programs. No physical frequency is measured.

**5. Falsification controls and reproducibility.** The primary `src/verify_d35.py` builds exact rational Gram products, derives all 780 labelled distances, validates every one of 741 frozen predecessor pairs, checks centered Gram and Laplacian ranks, and tests four deliberately invalid configurations: damaged distance, duplicate point, non-orthogonal tilt, and histogram-preserving label swap. The independent `tests/independent_d35.py` derives the labelled metric using Python `fractions.Fraction` structural rules **without importing the main Gram generator**. It checks all present pairs, parent SHA, conserved labels and exact ranks. Both output the same deterministic SHA-256 of `results/certificate_d35.json`.

**Scientific boundaries.** Orthogonal Euclidean extensions and uniform graph joins are established techniques. This one reproducible example is not a new universal theorem, evidence of additional physical spacetime dimensions, observed particles or measured resonances. There has been no independent human peer review.
