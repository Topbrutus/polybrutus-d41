# Polybrutus 31D — Exact proof and independently reproduced checks

**Scientific predecessor:** `Topbrutus/polybrutus-d30@v0.1.0`, commit `3c7ce4019a49ff4c54adb9271367aff5c62c1326`, DOI `10.5281/zenodo.24411018`, certificate SHA-256 `96c2f0693d849b2bab83c7ef4b25f647ecd3ee3a3ae460d3a32c5e9a73f1fb1c`.

**Conservation.** A zero-coordinate embedding of 35 labelled norm-one vectors preserves all 595 pairwise squared distances exactly. One new perpendicular norm-one vertex adds 35 squared distances `2R²`, for `C(36,2)=630` labelled distances.

**Gram rank and PSD.** The exact rational Gram G is PSD by construction with three original core dimensions, four centered simplex dimensions, and 24 independent perpendicular unit axes. Its exact rank is 31. The full labelled distance matrix D uniquely determines `B=-(1/2)HDH=HGH`, centered Euclidean rank 31. The original canonical first and fifth points are antipodal and their affine span contains the origin. Consequently no Euclidean realization of the complete labelled metric has affine dimension less than 31.

**Arithmetic invariants.** The normalized distance classes `k=4d²/R²` form the histogram `{4:6,7:2,8:599,9:3,10:10,12:6,13:2,16:2}`, totaling 630. The exact sum of squared distances at R=1 is `5087/4`.

**Chosen ideal graph spectrum.** For deliberately stipulated `w_ij=1/d_ij²` and unit masses, exact rational ranks give eigenvalue 18 with multiplicity 25 and eigenvalue 35/2 with multiplicity 4. The numbers are not experimentally measured frequencies or physical spacetime dimensions.

**Independent computational reproducibility.** `src/verify_d31.py` constructs all rational Gram entries, derives all 630 squared-distance labels, compares all 595 predecessor values and checks centered rank and graph eigenspaces. `tests/independent_d31.py` reconstructs the same complete labelled metric using Python `fractions.Fraction` structural rules without importing the primary Gram construction, then checks predecessor integrity, both ranks and eigenmodes. Both print the same deterministic certificate SHA-256. Four negative controls in the main verifier reject a corrupted distance, a repeated vertex, non-orthogonal tilt and histogram-preserving mislabelling; the independent checker also rejects corruption and label permutation.

**Limitations.** Orthogonal Euclidean extension and uniform graph join are established techniques; this is one exact, falsifiable mathematical case study, not a novel universal theorem, a measured physical signal, or independently human-peer-reviewed research.
