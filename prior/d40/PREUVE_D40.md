# Polybrutus D40 — Constructive exact proof and negative controls

**Frozen parent:** `Topbrutus/polybrutus-d39@v0.1.0`, scientific commit `a43cd3c67017fb6162005fd521c156cb78795beb`, DOI `10.5281/zenodo.24413712`, original certificate SHA-256 `8cd6731c294537d74c04a3501fe0da76a2cee6d496a2d2d5dfcfb0a26bb9304f`.

**Complete labelled conservation.** Embed each predecessor unit vector with one extra coordinate equal to zero. Every one of the **946** previous squared-distance pair labels is conserved. Append one perpendicular unit vector contributing **44** new squared distances of `2R²`. Exactly **990=C(45,2)** unordered labelled squared-distance pairs result.

**Euclidean PSD and minimal realization rank.** The full rational Gram `G` contains independent canonical rank-three, centered simplex rank-four, and 33 orthogonal unit-axis blocks. It is PSD by explicit Gram construction, has unit diagonal and exact rank **40**. The uniquely reconstructed centered Gram `B=-(1/2)HDH=HGH` has exact rank **40**. The original antipodal canonical pair contains the origin in its affine hull; hence this *specific complete labelled metric* cannot be embedded in fewer than 40 affine Euclidean dimensions.

**Exhaustive arithmetic.** The normalized squared-distance histogram `k=4d²/R²` is `{4:6,7:2,8:959,9:3,10:10,12:6,13:2,16:2}` and sums to 990. The exact sum of all squared distances at R=1 is **7967/4**. All 946 original values are compared by label, not merely by their unlabelled histogram.

**Separately stipulated mathematical graph.** Select scalar unit masses and complete symmetric edge weights `w_ij=1/d_ij²`. Exact rational matrix nullities establish weighted Laplacian eigenvalue **45/2** of multiplicity **34** and **22** of multiplicity **4**. These dimensionless modes are properties of one chosen mathematical graph and not experimental hertz.

**Two distinct programs and adversarial controls.** Primary `src/verify_d40.py` builds the rational Gram, derives all 990 labelled squared distances, validates every one of 946 predecessor labels against the SHA-256-checked frozen certificate, checks centered rank and ideal graph eigenspaces, and rejects four deliberately false inputs: corrupted pair distance, duplicate point, nonorthogonal tilt, and histogram-preserving relabelling. Independent `tests/independent_d40.py` reconstructs all 990 exact pair distances with Python `fractions.Fraction`, **without importing the primary Gram generator**, independently validates the parent's SHA-256 and 946 preserved labels, and verifies Euclidean and graph ranks. Both output the same deterministic certificate SHA-256.

**Scientific boundaries:** Orthogonal Euclidean extensions and uniform graph joins are established techniques. This result does not establish a new universal theorem, additional physical spacetime dimensions, observed particles, measured resonances or independent human peer review.
