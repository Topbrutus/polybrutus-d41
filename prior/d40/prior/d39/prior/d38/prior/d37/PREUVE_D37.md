# Polybrutus D37 — Exact constructive proof and independent verification

**Frozen parent:** `Topbrutus/polybrutus-d36@v0.1.0`, scientific commit `b1f8be2cdcdcae0d9911d90381a7b1304c84abf0`, DOI `10.5281/zenodo.24413083`, predecessor certificate SHA-256 `b189f74f7b38b378fa8f34011660817019c02601e5889378103fd9424207b331`.

**Exact conservation.** The zero-coordinate embedding of 41 frozen unit points preserves each of 820 predecessor labelled squared distances; the newly appended perpendicular unit point contributes 41 additional squared distances of `2R²`, yielding **861 = C(42,2)** unordered labelled pairs.

**Constructive positive semidefiniteness and minimal rank.** The exact rational Gram `G` has independent rank-three canonical, rank-four centered regular-simplex, and 30 unit-coordinate blocks. Thus `G` is positive semidefinite, has unit diagonal and rank **37**. The complete pair matrix determines centered Euclidean Gram `B=-(1/2)HDH=HGH` uniquely, and exact calculation gives rank **37**. The original antipodal pair ensures the origin belongs to the affine hull; this proves minimum Euclidean affine dimension 37 for this particular complete labelled metric.

**Full arithmetic invariants.** The normalized squared-distance histogram `k=4d²/R²` is `{4:6,7:2,8:830,9:3,10:10,12:6,13:2,16:2}`, totaling 861 pairs. At R=1 all squared distances sum to `6935/4`. Every old pair is checked by label, a stronger condition than matching just an unlabeled histogram.

**Selected ideal graph.** Choose unit scalar masses and complete graph weights `w_ij=1/d_ij²`. This adds uniform coupling 1/2 from the new point to each older point. Independent rational Laplacian rank calculations verify eigenvalue **21** of multiplicity **31** and **41/2** of multiplicity **4**. They are theoretical dimensionless modes of the stipulated graph, not observed physical frequencies.

**Independent programs and falsification controls.** Primary `src/verify_d37.py` builds exact rational Gram entries and all 861 labelled squared distances, checksum-verifies the frozen predecessor, checks all 820 parent labels, centered rank and ideal graph eigenspaces. Four deliberately invalid configurations are rejected: modified pair distance, duplicate point, nonorthogonal tilt and histogram-preserving relabelling. Independent `tests/independent_d37.py` reconstructs all labelled pair distances using Python `fractions.Fraction` rules **without importing the primary Gram generator**, validates the exact frozen parent checksum, all parent labels, Euclidean and graph ranks. Both print the same deterministic SHA-256 for `results/certificate_d37.json`.

**Scope.** Standard orthogonal Euclidean extensions and uniform graph joins construct one reproducible synthetic mathematical case. No physical additional spacetime dimension, experimental resonance, new particle, universal new theorem or independent human peer review is claimed.
