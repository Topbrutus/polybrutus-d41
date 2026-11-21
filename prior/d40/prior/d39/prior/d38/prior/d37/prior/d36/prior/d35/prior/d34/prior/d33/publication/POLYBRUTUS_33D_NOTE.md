# Polybrutus 33D: Exact 38-point Euclidean distance geometry and a twenty-sevenfold ideal graph mode

**Gabriel St-Pierre** - Independent researcher  
**20 November 2026**  
**Exact mathematical case study / scientific preprint**  
**AI assistance:** OpenAI ChatGPT (Astra) materially assisted symbolic calculations, two separately implemented checkers, manuscript preparation and controlled release. AI is not a scientific author.

## Abstract

We provide an exact Euclidean labelled squared-distance metric on 38 unit vertices by conservatively adding one perpendicular unit point to the frozen Polybrutus 32D configuration. All 666 predecessor labelled squared distances remain identical, while 37 new squared distances equal 2 at unit radius, for 703 unordered labelled pairs. An exact rational centered Euclidean Gram matrix is positive semidefinite of rank 33, establishing minimum affine realization dimension 33 for this complete labelled metric. Separately, in one stipulated unit-mass complete weighted graph with w_ij=1/d_ij², the Laplacian has eigenvalue 19 of exact multiplicity 27 and eigenvalue 37/2 of multiplicity 4. Two independently written exact-rational programs reproduce all labels, conservation, Gram ranks and ideal graph modes. These statements are mathematical, not evidence of 33-dimensional physical spacetime or measured resonances.

## 1. Frozen predecessor and conservative construction

The frozen predecessor is Topbrutus/polybrutus-d32 at scientific tag v0.1.0, commit cb63ed7ab447b13e251d92571fe0f1062db31cca, DOI 10.5281/zenodo.24411375. Its certificate records 37 unit vertices and 666 labelled squared distances, SHA-256 66b7605ab168d99ee33014ed1d61b5fbb3915c32041bdf44dd529c0bb56bb016. Its complete scientific tag source and all ancestor documents are preserved under prior/d32/ without nested Git metadata.

For each old vector p_i, i=1,...,37, define q_i=(p_i,0). Append new perpendicular norm-one point q38=e34. Every former dot product is unchanged. By orthogonality, ||q38-q_i||²=||q38||²+||q_i||²-2<q38,q_i>=2. All old labelled squared-distance values are exactly conserved, and the number of unordered pairs equals C(38,2)=703=666+37.

## 2. Exact rational Euclidean Gram and affine dimension

Seven original canonical unit vertices form a rational rank-three Gram block whose scalar products are 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vectors form an independent rank-four Gram block, with diagonal product 1 and distinct-pair products -1/4. Twenty-six mutually perpendicular norm-one unit axes are orthogonal to both other blocks. The full 38-by-38 rational Gram G is positive semidefinite, has unit diagonal and rank 3+4+26=33.

For every distinct pair define D_ij=2-2G_ij, with D_ii=0. Let H=I-(1/38)11^T be the centering matrix. The exact centered Euclidean distance Gram is B=-(1/2)HDH=HGH. Exact rank computation yields rank(B)=33. The first and fifth original canonical points are antipodal, so their affine span contains the origin. Because a complete labelled squared-distance matrix uniquely determines its centered Gram B, every Euclidean realization of this metric requires at least 33 affine dimensions.

## 3. Exhaustive labelled pair invariants

For R=1, normalize squared distances by k=4d²/R². The exact eight-class histogram is k=4 six times, k=7 twice, k=8 672 times, k=9 three times, k=10 ten times, k=12 six times, k=13 twice, and k=16 twice. The counts total 703. The exact sum of all unordered squared distances is 5671/4.

Both checkers compare each of the 666 predecessor labels against the checksum-verified frozen certificate; this is stronger than comparing histograms. A permutation of two labels can leave the multiset of distances unchanged while violating the exact predecessor correspondences, and is included among the negative controls.

## 4. Independently chosen ideal weighted graph

We separately stipulate one complete unit-mass weighted graph with edges w_ij=1/d_ij². The new perpendicular unit point joins every previous vertex at weight 1/2. The associated weighted graph Laplacian satisfies the exact uniform-join block form:

L33 = [ L32 + (1/2)I37      -(1/2)1_37 ]
      [ -(1/2)1_37^T             37/2 ].

For the characteristic polynomial chi_n(z), the uniform-join recurrence is chi38(z)=z(z-19)chi37(z-1/2)/(z-1/2). With the original canonical seven-point nonzero polynomial factor f6(z)=chi_P7(z)/z, the inherited expression is chi38(z)=z(z-19)^27(z-37/2)^4*f6(z-31/2).

Two independently implemented rational rank computations give nullity(L33-19I)=27 and nullity(L33-(37/2)I)=4. These are dimensionless ideal-model graph eigenvalues, not measured sound, real particles, or additional physical spacetime dimensions.

## 5. Reproducibility and adversarial checks

The primary program src/verify_d33.py constructs exact rational Gram entries, all 703 labelled squared pair distances, compares each of 666 old pairs with the frozen source, and computes exact Euclidean centered Gram and graph nullities. Four deliberate counterexamples are rejected: a corrupted distance, a duplicated vertex, a nonperpendicular tilt and a relabelling that preserves the unlabeled histogram.

The independently implemented tests/independent_d33.py reconstructs all pair labels using Python fractions.Fraction rules without importing the main Gram construction. It independently checks the frozen parent certificate SHA-256, every predecessor label, exact centered Gram rank, and graph eigenspace dimensions. Both programs print identical SHA-256 for results/certificate_d33.json. GitHub Actions verifies the same steps on the published scientific source and tag, including deterministic certificate output.

## 6. Scientific scope and references

Orthogonal Euclidean extension and uniform weighted-graph joining are established mathematics; this is one reproducible exact metric and a chosen ideal graph, not a novel universal theorem. No experiments establish new physical spacetime dimensions, measured resonances or particles. Independent human peer review has not occurred.

References: Gabriel St-Pierre, Polybrutus 32D (2026), DOI 10.5281/zenodo.24411375; source github.com/Topbrutus/polybrutus-d32, tag v0.1.0; original canonical geometry github.com/Topbrutus/polybrutus.
