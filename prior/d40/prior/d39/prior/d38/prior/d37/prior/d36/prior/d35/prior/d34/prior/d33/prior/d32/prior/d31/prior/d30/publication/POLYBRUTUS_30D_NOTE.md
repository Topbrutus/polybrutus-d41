# Polybrutus 30D: Exact 35-point Euclidean metric and twenty-fourfold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**20 November 2026**  
**Mathematical case study / scientific preprint**  
**AI disclosure:** OpenAI ChatGPT (Astra) substantially assisted exact calculation, independent verification programs, drafting, and publication; AI is not an author.

## Abstract

This note describes an exact 35-point Euclidean labelled squared-distance metric produced by conservatively extending the frozen Polybrutus 29D scientific source with one perpendicular unit vertex. All 561 predecessor labelled squared distances are identical, while 34 new squared distances equal 2 at unit radius, giving 595 unordered labelled pairs. Its exact centered Euclidean Gram matrix is positive semidefinite of rank 30, proving minimum affine realization dimension 30. Separately, a stipulated dimensionless complete graph with scalar unit masses and weights w_ij=1/d_ij² has Laplacian eigenvalue 35/2 of exact multiplicity 24 and eigenvalue 17 of multiplicity 4. Two independently implemented exact-rational verifiers reproduce the metric, predecessor conservation, affine rank, and model graph eigenspaces. These are mathematical properties, not measured physical dimensions or resonances.

## 1. Immutable predecessor and new point

The frozen predecessor is Topbrutus/polybrutus-d29 at scientific tag v0.1.0, exact source commit cbe179697cc1a287fc0637912ea6a6fd890eeca4, DOI 10.5281/zenodo.24410725. Its immutable certificate has SHA-256 d718f784bf20dccc1b8c08821466b352316054900fba20705bc28d79763a4731, recording 34 norm-one points and 561 exact labelled squared distances. All source tag files, scripts and ancestor records are reproduced unchanged under prior/d29/.

Embed each old p_i as q_i=(p_i,0) for i=1,...,34 in a larger Euclidean coordinate space. Append a new perpendicular norm-one basis vector q35=e31. Then ||q35-q_i||²=||q35||²+||q_i||²-2<q35,q_i>=2 for every previous label i. Appending zero to the former vectors does not change any old squared distance. The 595 pairs therefore partition into 561 conserved old pairs and 34 new pairs with squared distance 2.

## 2. Constructive rational Gram and minimum affine rank

Seven original canonical points form an exact rank-three Gram block with inner products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vertices have rank-four Gram block, diagonal 1 and off-diagonal -1/4. Twenty-three additional mutually perpendicular unit axes provide 23 independent directions, with zero cross-block products. Thus the rational Gram matrix G has unit diagonal, is positive semidefinite, and has rank 3+4+23=30.

Set D_ii=0 and D_ij=2-2G_ij for distinct labels. With H=I-(1/35)11^T the centered Euclidean distance Gram matrix satisfies B=-(1/2)HDH=HGH. Two canonical points are antipodal, placing the origin in their affine hull. Exact rank calculations give rank(B)=30. Because a complete labelled squared-distance matrix determines B uniquely, every Euclidean realization of this metric requires affine dimension at least 30.

## 3. Exact labelled-distance arithmetic

For unit radius R=1, normalized squared-distance class k=4d²/R² has eight frequencies: 4 appears 6 times, 7 appears 2 times, 8 appears 564 times, 9 appears 3 times, 10 appears 10 times, 12 appears 6 times, 13 appears 2 times, and 16 appears 2 times. These sum to 595. The exact sum of unordered squared distances is 4807/4.

The main and independent programs compare all 561 predecessor labels individually to the frozen original certificate, not merely an unlabeled histogram. Deliberate relabelling can preserve a histogram while corrupting labelled correspondence; the verification controls reject it.

## 4. A separately stipulated ideal weighted graph

In one chosen mathematical model, each vertex has scalar mass one, and every edge has weight w_ij=1/d_ij². The new point joins all 34 predecessor vertices with weight 1/2. The model graph Laplacian is the uniform block join:

L30 = [ L29+(1/2)I34    -(1/2)1_34 ]
      [ -(1/2)1_34^T         17  ].

Its characteristic polynomial obeys chi35(z)=z(z-35/2)chi34(z-1/2)/(z-1/2). In the inherited canonical factor form, chi35(z)=z(z-35/2)^24(z-17)^4*f6(z-14), with f6(z)=chi_P7(z)/z representing the nonzero seven-vertex canonical model factor.

Two exact-rational matrix-rank computations independently give nullity(L30-(35/2)I)=24 and nullity(L30-17I)=4. These eigenvalues are dimensionless modes of a stipulated weighted graph, not experimentally measured sound, particles, or spacetime geometry.

## 5. Reproducibility and adversarial tests

The primary verifier src/verify_d30.py builds every exact Gram entry, derives all 595 labelled squared pair values, checks the 561 frozen predecessor values, and computes centered Gram and ideal graph eigenspace dimensions. It also rejects four invalid configurations: a damaged distance, a duplicated point, a nonperpendicular tilt and an index permutation preserving histogram counts.

The independent verifier tests/independent_d30.py reconstructs the labelled metric separately using Python fractions.Fraction and rules for the three coordinate blocks, without importing the primary Gram generator. It independently verifies the predecessor SHA-256, all pair labels, centered Euclidean rank, and both graph nullities.

For reproduction, install Python and sympy==1.14.0, run python src/verify_d30.py and python tests/independent_d30.py. Both emit the same deterministic SHA-256 for results/certificate_d30.json. GitHub Actions also runs both scripts, checks frozen predecessor integrity and rejects unexpected certificate edits.

## 6. Limitations and references

Orthogonal Euclidean extension and uniform joins of weighted graph Laplacians are established mathematical techniques. This work records one reproducible specific metric and one stipulated graph, not a new universal theorem. It does not demonstrate 30-dimensional physical spacetime, measure a real resonance or identify a particle. There has not been independent human peer review.

References: Gabriel St-Pierre, Polybrutus 29D (2026), DOI 10.5281/zenodo.24410725; frozen source github.com/Topbrutus/polybrutus-d29 at v0.1.0; canonical geometry github.com/Topbrutus/polybrutus.
