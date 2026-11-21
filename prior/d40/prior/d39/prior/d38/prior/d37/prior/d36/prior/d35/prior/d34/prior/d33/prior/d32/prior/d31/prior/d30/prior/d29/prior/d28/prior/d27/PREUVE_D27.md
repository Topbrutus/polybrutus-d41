# Polybrutus D27 — Constructive exact proof and adversarial verification

**Frozen predecessor:** Polybrutus D26 scientific tag `v0.1.0`, source commit `5334c0fae387338314ea6e0fb8450bb3c1a5fcfa`, DOI `10.5281/zenodo.24410491`; SHA-256 of certificate `1b109ad0e5e7df765a00ca3fc5566d9d3381c96be1dd2e83bfe76c049159707a`. The parent archive is under `prior/d26/`.

**1. Conservation.** A perpendicular unit basis point added to 31 parent norm-one points yields 31 new squared distances each `2R²`, while every one of the 465 frozen previous labelled squared distances is unchanged.

**2. Minimal Euclidean affine dimension.** Construct the rational PSD Gram matrix `G` with three independent core dimensions, four simplex dimensions and twenty mutually orthogonal axis dimensions. Hence `rank(G)=27`, with diagonal ones and `D_ij=2-2G_ij`. The centered squared-distance Gram `B=-(1/2)HDH=HGH` has exact rank 27. Two core points are antipodal; the origin lies in their affine hull. Any realization of the same complete labelled distances must have affine rank 27.

**3. Arithmetic invariants.** The 496 pair labels partition into 465 preserved old distances and 31 new distances; normalized class histogram `{4:6,7:2,8:465,9:3,10:10,12:6,13:2,16:2}`. The exact squared-distance sum for unit radius is `4015/4`.

**4. Specified ideal graph.** Under deliberately chosen complete weights `w_ij=1/d_ij²`, the exact-rational Laplacian rank verifications show `nullity(L27-16I)=21` and `nullity(L27-(31/2)I)=4`, both reproduced independently. Neither is a measured vibration frequency.

**5. Independent negative controls.** `src/verify_d27.py` constructs squared distances via rational Gram scalar products; `tests/independent_d27.py` separately computes labelled distances from `fractions.Fraction` structure, without importing the primary Gram generator. Both check the complete 465-pair frozen parent map, all 496 present pairs, exact ranks and graph modes. The main program rejects a corrupted squared distance, duplicate vertex, non-orthogonal tilt, and point relabelling that preserves the histogram; the secondary also rejects corrupted pairs and a histogram-preserving relabelling.

**Limits:** Orthogonal extension and uniform weighted graph joins are standard mathematics. This is one reproducible synthetic metric, not a universal new theorem, demonstration of new physical spacetime dimensions, observation of resonances, or human peer-reviewed proof.
