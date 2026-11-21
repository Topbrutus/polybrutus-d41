# Polybrutus 29D: Exact 34-vertex Euclidean distance geometry and a twenty-threefold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**20 November 2026**  
**Exact mathematical case study / preprint**  
**AI assistance:** OpenAI ChatGPT (Astra) materially assisted exact symbolic computation, separate verification programs, manuscript preparation and release operations. AI is not a scientific author.

## Abstract

We construct a particular exact Euclidean labelled squared-distance metric on 34 unit points by appending one perpendicular unit vector to the frozen Polybrutus 28D scientific configuration. All 528 predecessor labelled squared distances are conserved exactly; 33 new squared distances equal 2 at unit radius. There are 561 unordered labelled pairs. An explicit rational positive-semidefinite centered Euclidean Gram matrix has rank 29, so the minimum Euclidean affine realization dimension of this full labelled metric is exactly 29. In a separately stipulated complete ideal weighted graph with unit scalar masses and weights w(i,j)=1/d(i,j)^2, the weighted graph Laplacian has eigenvalue 17 of exact multiplicity 23 and eigenvalue 33/2 of multiplicity 4. Two separately implemented rational verifiers check all labels, centered ranks, and selected graph mode multiplicities. No physical spacetime extension or measured resonance is inferred.

## 1. Frozen source and conservative orthogonal extension

The predecessor is Topbrutus/polybrutus-d28 at frozen scientific source tag v0.1.0, commit 496e6373abbabfad5c2e99adfcb3994ce8e39531, DOI 10.5281/zenodo.24410652. The immutable parent certificate contains 33 norm-one unit points and 528 labelled squared distances. Its SHA-256 is d21c21645bc96aa480c840fa24e12e755f6a21d8df8a1d51dac989c7caed1788. The complete source tagged v0.1.0 is preserved without Git metadata under prior/d28/.

Embed each old point p_i as q_i=(p_i,0) for i=1,...,33. Append one new perpendicular norm-one point q34=e30. By Euclidean orthogonality, ||q34-q_i||^2=||q34||^2+||q_i||^2-2<q34,q_i>=2 for every predecessor. The zero-coordinate embedding leaves every old pair squared distance unchanged, so the 561 pairs consist of 528 frozen old values and 33 new values equal to 2.

## 2. Exact rational Gram matrix and minimal dimension

The original seven canonical points define a rank-three Gram block with scalar products computed as 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vectors yield a rank-four block with diagonal product one and all distinct-pair inner products -1/4. Twenty-two additional perpendicular norm-one basis vectors provide 22 further independent directions, with all cross-block products zero.

Thus the complete rational 34-by-34 Gram matrix G has diagonal one, is positive semidefinite, and has rank 3+4+22=29. Define D as the full labelled squared-distance matrix with D_ii=0 and D_ij=2-2G_ij for distinct points. Write H=I-(1/34)11^T for the centering projection. Then B=-(1/2)HDH=HGH is positive semidefinite with exact rank 29. The first and fifth original core vertices are antipodal, so the origin lies in the affine hull. Consequently affine rank 29 is the minimum dimension of every Euclidean realization of this same complete labelled metric, rather than merely the number of chosen storage coordinates.

## 3. Exhaustive labelled pair invariants

For radius R=1, use normalized squared-distance index k=4d^2/R^2. The exact eight-class histogram reads: 4 occurs 6 times; 7 occurs 2 times; 8 occurs 530 times; 9 occurs 3 times; 10 occurs 10 times; 12 occurs 6 times; 13 occurs 2 times; 16 occurs 2 times. The class counts sum to 561. The sum of all 561 squared pair distances equals 4535/4.

Conservation is checked on the full labelled pair map rather than solely through histograms. A permutation of labels can leave the distance multiset unchanged while damaging the labelled correspondence, and this test is explicitly included in the negative controls.

## 4. Chosen ideal weighted graph and eigenmode degeneracy

Specify a complete graph on the 34 vertices, scalar mass one at each vertex, and weights w_ij=1/d_ij^2 for every distinct pair. This is a mathematical modeling choice, not a physical coupling measured in a laboratory. The new perpendicular unit point connects uniformly to all 33 predecessor vertices at weight 1/2. Therefore the Laplacian has the exact uniform-join form:

L29 = [ L28+(1/2)I33     -(1/2)1_33 ]
      [ -(1/2)1_33^T          33/2  ].

Its characteristic polynomial obeys chi34(z)=z(z-17)chi33(z-1/2)/(z-1/2). With f6(z)=chi_P7(z)/z denoting the nonzero canonical seven-point ideal-graph factor, the inherited decomposition is chi34(z)=z(z-17)^23(z-33/2)^4*f6(z-27/2).

Two independent exact-rational Laplacian rank computations yield nullity(L29-17I)=23 and nullity(L29-(33/2)I)=4. These are dimensionless spectral modes of the explicitly stipulated graph model. They are not observed sound frequencies in hertz, particles, or proof of extra physical spacetime dimensions.

## 5. Independent exact reproducibility and adversarial controls

The main script src/verify_d29.py builds the exact rational Gram matrix and computes all 561 labelled squared distances; it verifies the unmodified parent certificate and conservation of all 528 old pair values, positive-semidefinite centered Gram rank, and the two spectral nullities. Four deliberately invalid cases are rejected, covering corruption of one distance, duplication of a vertex, a nonorthogonal tilt and a point-label permutation which keeps histogram counts unchanged.

The secondary tests/independent_d29.py independently reconstructs the full 561-pair labelled metric from Python fractions.Fraction rules without importing the primary Gram generator. It also checks the parent certificate SHA-256 and every frozen pair label, centered Gram rank and model graph multiplicities. Both scripts emit identical SHA-256 values for results/certificate_d29.json.

To reproduce: install Python and sympy==1.14.0, execute python src/verify_d29.py and python tests/independent_d29.py. GitHub Actions executes the same pair of checks on the scientific source and tag and verifies that the deterministic certificate remains unchanged after recomputation.

## 6. Scientific limitations and references

Orthogonal Euclidean extensions and uniformly joined graph Laplacians are known constructions. This particular case study shows exact results for a selected metric and model, not a novel universal mathematical theorem. It does not demonstrate physical 29-dimensional spacetime, measured oscillatory modes, particles, or independent human peer review.

References: G. St-Pierre, Polybrutus 28D (2026), DOI 10.5281/zenodo.24410652; source: github.com/Topbrutus/polybrutus-d28, v0.1.0; canonical source: github.com/Topbrutus/polybrutus.
