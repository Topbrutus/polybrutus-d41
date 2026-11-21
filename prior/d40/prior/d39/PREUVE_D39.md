# Polybrutus D39 — Constructive exact proof and negative checks

**Frozen parent:** `Topbrutus/polybrutus-d38@v0.1.0`, science commit `ac8fd112fcd314092f106b210162ba1cc9879e9b`, DOI `10.5281/zenodo.24413606`, predecessor certificate SHA-256 `3be0e4aaa1487e82f35a97a3521abe70b05992888c550c8fdd95b5d50f306dbc`.

**Exact labelled conservation.** Embed all 43 predecessor unit points with a zero final coordinate. Their 903 labelled squared pair distances remain exactly the same. Add a perpendicular unit point giving 43 new squared distances equal to `2R²`. Exactly `C(44,2)=946` unordered labelled pairs result.

**Minimum affine Euclidean dimension.** Construct a rational Gram `G` comprising canonical seven-point rank three, centered five-point regular-simplex rank four, and 32 independent perpendicular unit axes. This Gram is PSD by explicit construction with diagonal one and rank 39. The unique centered Euclidean Gram `B=-(1/2)HDH=HGH` has exact rank **39**. Original antipodal points put the origin in the affine hull. No Euclidean realization of this complete labelled metric can have smaller affine dimension.

**Full exact arithmetic.** Histogram of `k=4d²/R²` is `{4:6,7:2,8:915,9:3,10:10,12:6,13:2,16:2}` and totals 946. All squared distances sum to `7615/4` at R=1. Every one of the 903 predecessor labels is checked against the exact frozen parent, not just the histogram.

**Separately chosen model.** For unit scalar masses and ideal weighted complete graph edges `w_ij=1/d_ij²`, exact Laplacian nullity proves eigenvalue **22** multiplicity **33**, and **43/2** multiplicity **4**. These are theoretical dimensionless modes, not measured frequencies.

**Two different checking implementations.** Primary `src/verify_d39.py` constructs rational Gram and all 946 labels, verifies 903 predecessor labels against the SHA-256 checked parent certificate, and checks exact Gram and Laplacian ranks. Four deliberately false alternatives are detected: corrupt pair distance, duplicate point, nonorthogonal tilt, and histogram-preserving relabelling. Independent `tests/independent_d39.py` reconstructs 946 pair distances using Python `fractions.Fraction`, **without importing the primary Gram generator**, verifying all 903 parent pairs and exact model nullities. Both generate matching SHA-256 of `results/certificate_d39.json`.

**Scientific limits.** Standard perpendicular Euclidean extensions and uniform graph joins establish one specific synthetic mathematical example, not new physical dimensions, measured resonances, observed particles, a novel universal theorem, or independent human peer review.
