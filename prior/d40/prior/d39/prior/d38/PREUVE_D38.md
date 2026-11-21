# Polybrutus 38D — Constructive exact proof and adversarial checks

**Immutable parent:** `Topbrutus/polybrutus-d37@v0.1.0`, scientific commit `294be7f814f2c8770af806833e98194fb192afff`, DOI `10.5281/zenodo.24413346`, certificate SHA-256 `f649594103c43927126174dc5d807777bb976bf4feae2e9a1c0ee73aa67f3b2b`.

**Exact label conservation:** Zero-coordinate embedding of 42 old unit points preserves each of 861 predecessor labelled squared distances. Adding one perpendicular unit point gives 42 new squared distances `2R²`, for a total of exactly `C(43,2)=903` pairs.

**Euclidean PSD and minimal dimension:** An explicit exact rational Gram `G` is formed from three independent canonical directions, four centered-simplex directions, and 31 additional orthogonal unit axes. It has unit diagonal, is PSD by construction, and has exact rank 38. The unique centered distance Gram `B=-(1/2)HDH=HGH` is exact PSD of rank 38; the original antipodal pair places the origin in the affine hull. Thus no Euclidean realization of this complete labelled squared-distance metric has affine dimension below 38.

**Full arithmetic invariants:** The normalized histogram `k=4d²/R²` is `{4:6,7:2,8:872,9:3,10:10,12:6,13:2,16:2}` and sums to 903. At R=1 the squared distances sum to `7271/4`. The predecessor conservation check uses pair labels rather than only the histogram.

**Chosen ideal graph:** With independently stipulated scalar unit masses and `w_ij=1/d_ij²`, the appended point joins uniformly at coupling 1/2. Exact rational matrix ranks give eigenvalue `43/2` of multiplicity 32 and eigenvalue `21` of multiplicity 4. No physical frequency is inferred.

**Two separate implementations and negative tests:** The primary `src/verify_d38.py` reconstructs all Gram entries and 903 labelled distances, checks all 861 parent labels against the checksum-verified frozen predecessor, centered ranks and selected model eigenspaces. Four adversarial primary controls reject a corrupted pair, duplicated vertex, non-orthogonal tilt and histogram-preserving label swap. Independently, `tests/independent_d38.py` derives all pair distances using Python `fractions.Fraction` structural rules **without importing the primary Gram generator**, verifies the predecessor certificate SHA-256, all old and new labels, Euclidean ranks and graph eigenspaces. Both print identical deterministic SHA-256 for `results/certificate_d38.json`.

**Scientific limits:** Orthogonal Euclidean extensions and uniform graph joins are established techniques. This proves facts about one particular synthetic mathematical construction, not a new universal theorem, physical spacetime dimensions, measured resonances, particles or independent human peer review.
