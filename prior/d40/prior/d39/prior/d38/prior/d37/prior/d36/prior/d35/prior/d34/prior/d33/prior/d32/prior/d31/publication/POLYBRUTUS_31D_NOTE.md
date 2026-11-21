# Polybrutus 31D: Exact 36-point Euclidean distance geometry and a twenty-fivefold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**20 November 2026**  
**Mathematical case study / scientific preprint**  
**AI assistance disclosed:** OpenAI ChatGPT (Astra) assisted exact symbolic mathematics, independent verification code, scientific writing and publication operations. AI is not a scientific author.

## Abstract

We construct an exact Euclidean labelled squared-distance metric on 36 unit vertices by conservatively appending one perpendicular unit point to the frozen Polybrutus 30D scientific configuration. All 595 old labelled squared distances remain unchanged; the 35 new squared distances are equal to 2 at unit radius, giving 630 unordered labelled pairs. Its rational centered Euclidean Gram matrix is positive semidefinite with exact rank 31, establishing the minimum affine dimension of a realization of this complete labelled metric. For a separately stipulated ideal equal-mass complete weighted graph with edge weights w_ij=1/d_ij², its Laplacian has eigenvalue 18 with exact multiplicity 25 and eigenvalue 35/2 with multiplicity 4. Two independently written exact-rational verifiers confirm pair labels, predecessor conservation, Euclidean ranks and these chosen graph modes. No physical spacetime dimension or measured resonance is inferred.

## 1. Immutable predecessor and perpendicular extension

Our frozen predecessor is Topbrutus/polybrutus-d30 scientific source tag v0.1.0, exact commit 3c7ce4019a49ff4c54adb9271367aff5c62c1326, DOI 10.5281/zenodo.24411018. Its certificate contains 35 unit vertices and 595 labelled squared distances. The original parent certificate SHA-256 is 96c2f0693d849b2bab83c7ef4b25f647ecd3ee3a3ae460d3a32c5e9a73f1fb1c. All frozen scientific tag files and their ancestral references are reproduced unmodified as ordinary files in prior/d30/.

Embed each prior unit point p_i into one more coordinate slot by q_i=(p_i,0), i=1,...,35. Append a unit vector q36=e32 perpendicular to every predecessor. The new squared distances are ||q36-q_i||²=||q36||²+||q_i||²-2<q36,q_i>=2. Zero-coordinate embedding leaves every predecessor pair unchanged. There are C(36,2)=630 pairs, partitioned into 595 preserved old pairs and 35 new pairs of squared distance 2.

## 2. Rational positive-semidefinite Gram construction

The original seven canonical unit vertices produce a rank-three rational Gram block, with scalar products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vertices produce an independent rank-four block, with diagonal products 1 and off-diagonal products -1/4. Twenty-four additional mutually orthogonal unit-axis points contribute 24 more independent directions with zero cross-block products. Thus the Gram matrix G is exact rational PSD with unit diagonal and rank 3+4+24=31.

For the complete labelled squared-distance matrix D, put D_ii=0 and D_ij=2-2G_ij for distinct points. Let H=I-(1/36)11^T. The exact Euclidean centered Gram B=-(1/2)HDH=HGH has rank 31. The first and fifth canonical points are antipodal, so their affine hull contains the origin. Because the complete labelled distance matrix uniquely determines B, every Euclidean realization of this same metric has minimum affine dimension exactly 31.

## 3. Labelled pair arithmetic invariants

For radius R=1, normalized squared-distance index k=4d²/R² has eight classes: 4 occurs 6 times; 7 occurs 2 times; 8 occurs 599 times; 9 occurs 3 times; 10 occurs 10 times; 12 occurs 6 times; 13 occurs 2 times; 16 occurs 2 times. The counts sum to 630. The exact sum of squared pair distances is 5087/4.

Conservation of all predecessor distances is verified for each labelled pair, not merely for the histogram. A label permutation may preserve the unlabeled histogram while violating the original labelled metric; this is explicitly detected by our adversarial controls.

## 4. A chosen ideal complete weighted graph

We independently select a mathematical graph model with unit scalar mass per vertex and complete edge weights w_ij=1/d_ij². The new perpendicular point is joined to all 35 previous vertices with coupling 1/2. Therefore its weighted graph Laplacian satisfies the uniform-join block identity:

L31 = [ L30+(1/2)I35    -(1/2)1_35 ]
      [ -(1/2)1_35^T        35/2  ].

The characteristic polynomial recurrence reads chi36(z)=z(z-18)chi35(z-1/2)/(z-1/2). In the inherited canonical seven-vertex factorization, this is chi36(z)=z(z-18)^25(z-35/2)^4*f6(z-29/2), where f6(z)=chi_P7(z)/z represents the nonzero canonical model factor.

Independent exact-rational calculations give nullity(L31-18I)=25 and nullity(L31-(35/2)I)=4. These are dimensionless modes of a selected graph and are not physical measurements of hertz, sound, particles or spacetime.

## 5. Two separately implemented verification programs

The primary program src/verify_d31.py constructs every exact rational Gram entry, derives all 630 labelled squared pair values, checks every one of 595 predecessor pairs from the frozen checksum-verified parent, and calculates exact centered Euclidean Gram and weighted Laplacian ranks. It includes four deliberately incorrect alternatives: a damaged distance, duplicated vertex, nonperpendicular tilt and histogram-preserving label swap.

A second independently implemented checker tests/independent_d31.py reconstructs the complete labelled distance map from Python fractions.Fraction rules, without importing the primary Gram implementation, and checks every old label, all present labels, centered rank and graph eigenspace dimensions. Both checkers generate the same SHA-256 for results/certificate_d31.json.

Reproduce with Python and sympy==1.14.0: run python src/verify_d31.py followed by python tests/independent_d31.py. GitHub Actions reruns the programs and checks that the committed certificate is unchanged.

## 6. Scientific limitations and references

Orthogonal embedding and uniformly joined weighted Laplacians are established mathematical constructions. This note documents one exact reproducible labelled-metric case study and one chosen ideal graph. It is not a novel universal theorem, direct evidence of 31-dimensional physical spacetime, a measured resonance or a particle; no independent human peer review has taken place.

References: Gabriel St-Pierre, Polybrutus 30D (2026), DOI 10.5281/zenodo.24411018; frozen source github.com/Topbrutus/polybrutus-d30, tag v0.1.0; canonical source github.com/Topbrutus/polybrutus.
