# Polybrutus D23 — Exact proof and independent controls

**1. Conservation.** Append an orthogonal unit vector to the 27 frozen D22 unit points. Every previous 351 labelled squared distances is preserved; exactly 27 new distances are 2. Two checkers compare all 351 old labels and values to `prior/d22/results/certificate_d22.json`.

**2. Euclidean metric.** A rational scalar-product Gram matrix `G` encodes three original orthogonal directions, four centered simplex directions, and 16 extra orthogonal directions. Each of the 28 unit vectors has squared norm 1, all points are distinct, and `D_ij=2-2G_ij`. The centered squared-distance Gram matrix satisfies `B=-(1/2)HDH=HGH`, hence is PSD. Its exact rank is 23, so any Euclidean realization of the labelled squared-distance metric has minimum affine dimension 23.

**3. Arithmetic invariants.** The primary checker independently verifies `C(28,2)=378`, the eight normalized distance classes `{4:6,7:2,8:347,9:3,10:10,12:6,13:2,16:2}` and `sum D_ij=3071/4` for `R=1`. The second checker derives every distance solely from rational construction rules using `fractions.Fraction`, never importing the primary Gram generator.

**4. Chosen graph model.** Defining complete edge weights `w_ij=1/D_ij` yields `L23` with old top-left block `L22+I/2`. The exact nullities `nullity(L23-14I)=17` and `nullity(L23-(27/2)I)=4` are established by independent rational rank computations in both programs, and cannot be mistaken for measured resonance. The ideal graph characteristic recurrence follows from the uniform join.

**5. Adversarial tests.** The primary checker rejects a corrupted labelled pair, a duplicated new vertex, a tilted new vertex, and a label-swap despite unchanged unlabeled distance histogram. The independent checker also detects corrupted pairs and an index permutation.

**6. Reproducibility.** Run `python src/verify_d23.py` then `python tests/independent_d23.py`. Assert that both print the identical `CERTIFICATE_SHA256` and that the output is unchanged in Git. GitHub Actions repeats these checks.

**Scope.** This is a reproducible particular synthetic metric. Orthogonal Euclidean extensions and uniform graph joins are established techniques; no novel universal theorem, measured wave, or additional physical spacetime dimension is asserted. Human peer review has not taken place.
