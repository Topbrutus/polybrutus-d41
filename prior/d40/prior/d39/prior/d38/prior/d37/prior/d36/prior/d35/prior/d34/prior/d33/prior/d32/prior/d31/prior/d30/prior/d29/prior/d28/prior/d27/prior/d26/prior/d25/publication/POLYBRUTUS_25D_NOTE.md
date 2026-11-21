# Polybrutus 25D: A reproducible 30-vertex Euclidean metric and a nineteenfold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**20 November 2026**  
**Mathematical case study / preprint**  
**AI disclosure:** OpenAI ChatGPT (Astra) substantially assisted exact verification, independent algorithms, manuscript writing and release. AI is not a scientific author.

## Abstract

This paper gives an explicitly reproducible Euclidean labelled squared-distance metric on 30 unit vertices, conservatively extended from the immutable Polybrutus 24D scientific source. The new point is orthogonal to the preceding 29 unit vertices, so it preserves all 406 previously labelled squared distances and adds exactly 29 distances of square 2. There are 435 unordered pairs in total. Its exact centered Euclidean Gram matrix is positive semidefinite with rank 25, proving that no Euclidean realization of these complete labelled squared distances has affine dimension less than 25. A separately stipulated complete ideal weighted graph with weights w(i,j)=1/d(i,j)^2 has Laplacian eigenvalues 15 (multiplicity 19) and 29/2 (multiplicity 4), independently certified using exact arithmetic. Neither the mathematical dimension nor the synthetic graph modes represent empirically observed physics.

## 1. Frozen predecessor and conservative construction

The frozen predecessor is Topbrutus/polybrutus-d24 at scientific tag v0.1.0, commit 16d9f0ee705d15c220bc02e18563a5e5c5c1e5fc, public Zenodo DOI 10.5281/zenodo.24410102. Its certificate provides 29 labelled unit points with 406 exact labelled squared distances. The D24 scientific files and certificate are preserved unaltered inside prior/d24/. Both new programs require an exact SHA-256 match of the predecessor's scientific certificate.

Embed every preceding p_i as q_i=(p_i,0) in a new 26-coordinate storage space, and append q_30=e_26 orthogonal to all preceding q_i. Then for each i=1,...,29, ||q_30-q_i||^2=||q_30||^2+||q_i||^2-2<q_30,q_i>=2. All earlier pairwise distances remain identical because the old points are merely embedded with a zero in the extra coordinate. Thus there are C(30,2)=435 total unordered pairs: 406 old, 29 new.

## 2. Exact Gram construction and minimum affine dimension

The original seven canonical unit vectors contribute a three-dimensional core with scalar products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vectors contribute a rank-four block, with inner products 1 on the diagonal and -1/4 for every pair of distinct simplex vertices. Eighteen mutually orthogonal unit basis vectors provide another 18 independent directions. The three blocks are mutually orthogonal and every Gram diagonal entry equals 1.

Therefore the explicit rational Gram matrix G is positive semidefinite and has rank 3+4+18=25. Let D denote the squared-distance matrix, D_ii=0 and D_ij=2-2G_ij for i not equal j. With H=I-(1/30)11^T, the centered Euclidean squared-distance Gram matrix obeys B=-(1/2)HDH=HGH. Two canonical core vertices are antipodal, so the origin belongs to the affine hull. Two independent exact rank computations show rank(B)=25. Since a complete labelled squared-distance matrix uniquely determines B, any Euclidean realization of this labelled distance metric has minimum affine dimension 25.

## 3. Invariants and independently verified distances

For radius R=1 the eight normalized squared-distance class counts, indexed by k=4d^2/R^2, are 4:6, 7:2, 8:404, 9:3, 10:10, 12:6, 13:2, and 16:2. They total 435 unordered vertex pairs; their sum of squared distances is exactly 3527/4.

The 406 ancestor labelled squared distances are checked against the source D24 JSON by full pair labels rather than only a histogram. A label permutation can leave the eight-class histogram unchanged while changing the correct pair-to-distance mapping. Four negative controls are evaluated by the main verifier; the separate Fraction-based program also checks deliberately corrupted labels and measurements.

## 4. Chosen ideal complete weighted graph

One possible mathematical graph model gives each point a scalar mass of one and each distinct pair an edge weight w_ij=1/d_ij^2. These couplings are stipulated rather than measured. The new point joins each preceding vertex with weight 1/2. Therefore the 30-by-30 Laplacian is the uniform join:

L25 = [L24+(1/2)I29   -(1/2)1_29]
      [-(1/2)1_29^T        29/2 ].

With chi_n(z) the characteristic polynomial, the uniform-join recurrence is chi30(z)=z(z-15)chi29(z-1/2)/(z-1/2). The preceding canonical factorization yields chi30(z)=z(z-15)^19(z-29/2)^4*f6(z-23/2), with f6(z)=chi_P7(z)/z for the model's canonical seven-vertex graph factor.

The exact rank tests independently return nullity(L25-15I)=19 and nullity(L25-(29/2)I)=4. These values are dimensionless eigenvalues of this prescribed graph. They are not frequencies in hertz, evidence of particles, or physical new dimensions.

## 5. Verification, adversarial controls and reproducibility

The first verifier, src/verify_d25.py, uses an exact rational Gram matrix and SymPy to check all 435 labelled squared distances and the 406 predecessor pairs, the centered Gram rank and the two spectral eigenspace dimensions. It rejects four deliberately false proposals, including one corrupted pair, one duplicated point, a non-orthogonal tilt, and a label swap that retains the overall distance histogram.

The second separately authored verifier, tests/independent_d25.py, uses the built-in Fraction rational datatype to reconstruct every labelled squared distance from structural rules, without importing the first program's Gram generator. It also independently checks the frozen D24 certificate and the Gram and Laplacian rank conditions.

To reproduce, install Python and sympy==1.14.0, then run python src/verify_d25.py followed by python tests/independent_d25.py. Both emit the identical SHA-256 of results/certificate_d25.json. A GitHub Actions workflow repeats these tests and rejects any nondeterministic change to the checked-in certificate.

## 6. Scientific limits

Orthogonal extension and the uniform complete-graph join are established mathematical techniques. The present manuscript documents a particular exact labelled metric and a deliberately chosen ideal graph; it does not establish a novel universal mathematical theorem. No physical extra spacetime dimension, measured oscillation, particle, or independent human peer review is claimed.

## References

1. Gabriel St-Pierre, Polybrutus 24D: Exact 29-Vertex Euclidean Metric and Eighteenfold Ideal Graph Mode (2026), https://doi.org/10.5281/zenodo.24410102.
2. Topbrutus/polybrutus-d24, frozen scientific source v0.1.0, https://github.com/Topbrutus/polybrutus-d24/tree/v0.1.0.
3. Topbrutus/polybrutus canonical seven-point specification, https://github.com/Topbrutus/polybrutus/blob/main/SPECIFICATION.md.
