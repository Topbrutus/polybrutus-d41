# Polybrutus 27D: Exact 32-point Euclidean distance geometry and twenty-onefold ideal graph modes

**Gabriel St-Pierre** — Independent researcher  
**20 November 2026**  
**Exact mathematical case study / scientific preprint**  
**AI contribution disclosure:** OpenAI ChatGPT (Astra) materially assisted symbolic mathematics, independent verification code, writing and release operations; AI is not a scientific author.

## Abstract

We describe a concrete exact Euclidean metric on 32 labelled unit vertices, obtained by appending one perpendicular unit point to the immutable Polybrutus 26D configuration. The extension preserves all 465 predecessor labelled squared distances and creates 31 new squared distances, all equal to 2 at unit radius. There are 496 unordered point pairs in total. The exact centered Euclidean squared-distance Gram matrix is positive semidefinite of rank 27. Consequently no Euclidean realization of this full labelled metric can have affine dimension less than 27. A separately stipulated unit-mass, complete weighted graph with weights w(i,j)=1/d(i,j)^2 has Laplacian eigenvalue 16 of exact multiplicity 21, and eigenvalue 31/2 of multiplicity 4. Two independently coded exact-rational verifiers reproduce all labelled distances, affine rank and these model-dependent spectral eigenspace dimensions. No physical spacetime dimension or measured resonance is inferred.

## 1. Immutable source and conservative orthogonal completion

The predecessor is Topbrutus/polybrutus-d26 at frozen scientific source tag v0.1.0, commit 5334c0fae387338314ea6e0fb8450bb3c1a5fcfa, DOI 10.5281/zenodo.24410491. The predecessor JSON certificate has 31 labelled unit points, 465 squared distances, and SHA-256 1b109ad0e5e7df765a00ca3fc5566d9d3381c96be1dd2e83bfe76c049159707a. The original frozen predecessor scientific files, verifiers and ancestor references are included unchanged under prior/d26/.

Define q_i=(p_i,0) for all 31 old norm-one vertices and append q32=e28, a norm-one vector orthogonal to all previous coordinates. For each old i, the squared distance is ||q32-q_i||²=||q32||²+||q_i||²-2<q32,q_i>=2. All old pairwise squared distances are conserved by the zero-coordinate embedding. Thus C(32,2)=496=465+31 labelled pairs; the new pairs each have squared distance 2.

## 2. Exact rational Gram representation and minimum affine dimension

The original seven canonical points occupy a rank-three block with pair products computed exactly as 3*x_i*x_j+y_i*y_j+z_i*z_j. The next five centered simplex unit points occupy a rank-four block; their self-products are 1 and the distinct-pair products are -1/4. Twenty further mutually orthogonal unit points contribute independent coordinate axes, with no cross-products to the preceding blocks. The exact 32-by-32 Gram matrix G is rational positive semidefinite with rank 3+4+20=27. Every point has squared norm one.

Define D as the full labelled squared-distance matrix, D_ii=0 and D_ij=2-2G_ij for distinct indices. The centering projection H=I-(1/32)11^T gives the Euclidean centered Gram matrix B=-(1/2)HDH=HGH. The first and fifth original core points are antipodal; this places the origin in the affine span. The exact rank of B is 27. Since a full labelled Euclidean squared-distance matrix determines the centered Gram matrix uniquely, every Euclidean realization of the same labelled metric has minimal affine dimension exactly 27.

## 3. Complete labelled-distance invariants

For unit radius R=1 the normalized squared-distance index k=4d²/R² has eight exact classes: 4 occurs 6 times; 7 occurs 2 times; 8 occurs 465 times; 9 occurs 3 times; 10 occurs 10 times; 12 occurs 6 times; 13 occurs 2 times; 16 occurs 2 times. These counts sum to 496. The sum of all unordered squared distances is exactly 4015/4.

Both verification programs compare all predecessor labels and values to the immutable source certificate rather than merely checking an unlabeled histogram. A permutation of vertex names can leave the histogram identical while changing a labelled metric; a deliberately false renaming is rejected by the adversarial tests.

## 4. Chosen complete weighted graph and ideal eigenspaces

In a deliberately stipulated mathematical graph model, every vertex has scalar mass one, and each distinct pair has edge weight w_ij=1/d_ij². The new perpendicular point joins all 31 predecessor points with uniform weight 1/2. Accordingly its weighted graph Laplacian has the exact uniform-join block form:

L27 = [ L26 + (1/2)I31    -(1/2)1_31 ]
      [ -(1/2)1_31^T          31/2  ].

The characteristic polynomial recurrence is chi32(z)=z(z-16)chi31(z-1/2)/(z-1/2). In the inherited canonical seven-vertex factorization the same spectrum takes the form chi32(z)=z(z-16)^21(z-31/2)^4*f6(z-25/2), where f6(z)=chi_P7(z)/z is the nonzero canonical graph factor for the chosen model.

Two exact-rational Laplacian rank calculations yield nullity(L27-16I)=21 and nullity(L27-(31/2)I)=4. They are **dimensionless eigenvalues of a mathematical model**, not measured vibrations, sound frequencies, particle states, or proof of physical dimensions.

## 5. Independent exact validation and adversarial tests

The primary verifier src/verify_d27.py constructs all exact rational Gram entries, derives all 496 labelled squared distances, checks preservation of all 465 frozen predecessor pairs, and computes centered Euclidean Gram and weighted Laplacian ranks. It rejects four deliberately invalid configurations: one corrupted distance entry, one duplicated vertex, one non-orthogonal tilted vertex, and a label swap preserving the unlabeled histogram.

The independent verifier tests/independent_d27.py separately reconstructs every labelled squared distance using Python fractions.Fraction structural rules without importing the primary Gram generator. It verifies the original predecessor SHA-256 and every old labelled pair, independently computes centered Euclidean and graph ranks, and compares the complete JSON certificate.

To reproduce the results, install Python with sympy==1.14.0, then run python src/verify_d27.py and python tests/independent_d27.py. Both print the same deterministic SHA-256 for results/certificate_d27.json. The accompanying GitHub Actions workflow reruns both and checks the committed certificate has no changes.

## 6. Scientific limitations

Orthogonal Euclidean completion and uniform graph joins are established mathematical techniques. This work gives one specific reproducible labelled metric and one selected ideal graph, not a novel universal theorem. It neither demonstrates 27-dimensional physical spacetime nor measures any oscillation, particle or resonant material. No independent human peer review has been completed.

## References

1. Gabriel St-Pierre, Polybrutus 26D: Exact 31-Point Euclidean Metric and Twentyfold Ideal Graph Mode, Zenodo, 2026, https://doi.org/10.5281/zenodo.24410491.
2. Topbrutus/polybrutus-d26, scientific source frozen at v0.1.0, https://github.com/Topbrutus/polybrutus-d26/tree/v0.1.0.
3. Topbrutus/polybrutus canonical specification, https://github.com/Topbrutus/polybrutus/blob/main/SPECIFICATION.md.
