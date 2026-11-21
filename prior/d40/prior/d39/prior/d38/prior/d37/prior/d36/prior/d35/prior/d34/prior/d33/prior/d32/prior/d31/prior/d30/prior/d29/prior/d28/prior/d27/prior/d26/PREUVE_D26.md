# Polybrutus D26 — Formal exact proof and negative controls

**Provenance:** Polybrutus D25 scientific tag `v0.1.0` at `0dda9ed43622385d3e5fdddef9336458a32ed5ed`, public DOI `10.5281/zenodo.24410435`. Exact D25 certificate SHA-256 `1cad349945e85a0700906762b21c3644ee9188308098b340e99acc525cac4ffd`. Unmodified scientific tag files archived under `prior/d25/`.

**1. Conservation.** Embed the 30 old norm-one vertices in one extra storage coordinate equal to zero, append a new perpendicular norm-one point. Every one of the 435 old labelled squared distances is identical; 30 additional squared distances are each `2R²`.

**2. Euclidean realization and minimal dimension.** The full exact rational Gram `G` is PSD, rank `3+4+19=26`, and has diagonal ones. The complete labelled Euclidean squared-distance matrix `D` satisfies `B=-(1/2)HDH=HGH`. Two antipodal canonical vertices place the origin in the affine hull. Therefore the exact rank of `B` is 26; every Euclidean realization of these complete labelled distances needs affine dimension at least 26.

**3. Exact combinatorics.** `C(31,2)=465`; normalized histogram `{4:6,7:2,8:434,9:3,10:10,12:6,13:2,16:2}`; `sum_{i<j} d_ij²=3767/4` for radius one.

**4. Chosen ideal graph spectrum.** Define `w_ij=1/D_ij` for every pair, with graph Laplacian uniform join to the parent. Exact-rational ranks yield nullities `20` for eigenvalue `31/2` and `4` for eigenvalue `15`. These are model-specific dimensionless values, not observations or real vibration frequencies.

**5. Independent proof checks.** The main script `src/verify_d26.py` derives all 465 squared distances from an exact Gram construction and compares all 435 predecessor-labelled distances. The second `tests/independent_d26.py` reconstructs each distance from independent `fractions.Fraction` rules and does not import the Gram generator. Both compute exact Euclidean centered rank and Laplacian eigenspace dimensions and print the identical deterministic certificate hash. The main checker rejects a wrong squared distance, duplicate vertex, tilted nonperpendicular point and histogram-preserving label permutation; the second also rejects corruption and a label permutation.

**Limits:** This is a particular synthetic mathematical metric and one stipulated graph model, not a novel universal theorem, a measured physical effect, a proof of spacetime dimensions, or human peer-reviewed research.
