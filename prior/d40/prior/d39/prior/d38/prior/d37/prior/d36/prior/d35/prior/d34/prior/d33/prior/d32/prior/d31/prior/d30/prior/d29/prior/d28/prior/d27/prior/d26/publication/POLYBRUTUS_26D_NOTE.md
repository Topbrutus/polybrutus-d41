# Polybrutus 26D: Exact 31-point Euclidean distance geometry and a twentyfold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**20 November 2026**  
**Mathematical case study / scientific preprint**  
**AI assistance disclosure:** OpenAI ChatGPT (Astra) substantively assisted exact rational mathematics, independent verification, drafting, and controlled scientific release. AI is not an author.

## Abstract

We construct a specific exact 31-point unit-sphere Euclidean labelled squared-distance metric via one conservative perpendicular extension of the frozen Polybrutus 25D release. The extension preserves all 435 parent labelled squared distances and adds 30 new squared distances equal to 2, for 465 total labelled pairs. An exact positive-semidefinite Euclidean centered Gram matrix has rank 26, hence every Euclidean realization of this complete labelled metric has affine dimension at least 26. In a separately stipulated dimensionless complete graph model with unit scalar masses and edge weights w_ij=1/d_ij², the weighted Laplacian has eigenvalue 31/2 of exact multiplicity 20 and eigenvalue 15 of multiplicity 4. Two separately implemented rational verifiers check the distance labels, Gram ranks and eigenspace dimensions. This is a mathematical construction, not physical evidence of extra spacetime dimensions or measured oscillations.

## 1. Frozen predecessor and extension method

The frozen predecessor is Topbrutus/polybrutus-d25 scientific tag v0.1.0, source commit 0dda9ed43622385d3e5fdddef9336458a32ed5ed (Zenodo DOI 10.5281/zenodo.24410435). The immutable D25 certificate contains 30 distinct unit vertices and 435 exact labelled squared distances. Its SHA-256 is 1cad349945e85a0700906762b21c3644ee9188308098b340e99acc525cac4ffd. The full source and original verifiers are archived unmodified under prior/d25/.

Embed the old vertices p_i as q_i=(p_i,0) for i=1,...,30 in one additional coordinate slot. Append the new orthogonal unit basis vertex q31=e27. The norm of every existing and new point is one, and each new cross product equals zero, so ||q31-q_i||²=1+1-0=2 for all parent vertices. All former pairwise squared distances remain unchanged under the embedding. Thus C(31,2)=465 labelled pairs partition into 435 conserved pairs and 30 new pairs of squared distance exactly 2.

## 2. Exact rational Gram construction and minimal affine rank

Seven canonical core vertices have pairwise scalar products derived exactly as 3*x_i*x_j+y_i*y_j+z_i*z_j in a rank-three block. Five centered regular-simplex unit vertices form an independent rank-four block: diagonal products 1 and off-diagonal products -1/4. Nineteen mutually orthogonal unit basis vectors occupy remaining independent axes. All three blocks have zero mutual cross products. Thus the Gram matrix G is positive semidefinite, rational, and has exact rank 3+4+19=26 with every diagonal entry one.

For the complete squared-distance matrix D, D_ii=0 and D_ij=2-2G_ij for distinct points. Define the centering matrix H=I-(1/31)11^T. The exact Euclidean centered Gram matrix equals B=-(1/2)HDH=HGH and has rank 26. The first and fifth core vertices are antipodal, so the origin is contained in the affine hull. Every Euclidean realization of a complete labelled squared-distance metric has the same unique centered Gram matrix B; consequently affine rank 26 is a minimal and exact fact of this case, rather than a coordinate artefact.

## 3. Exact combinatorial and arithmetic invariants

At unit radius, the normalized squared-distance index k=4d²/R² has the eight-class histogram: 4:6, 7:2, 8:434, 9:3, 10:10, 12:6, 13:2, 16:2. These add to 465. The exact sum of all unordered pair squared distances is 3767/4.

Crucially, both verification programs check pair labels individually against the frozen D25 JSON certificate, rather than only comparing class counts. A relabelling of the points can preserve the overall distance histogram while violating labelled pair conservation; this deliberate counter-test is detected.

## 4. Deliberately chosen weighted complete graph and model spectrum

A particular theoretical complete graph assigns each vertex a scalar mass one and each distinct pair edge weight w_ij=1/d_ij². These are stipulated mathematical graph couplings, not physical measurements. The appended perpendicular unit vertex has edge weight 1/2 to every old vertex. Thus its Laplacian has the uniform-join block form:

L26 = [ L25 + (1/2)I30    -(1/2)1_30 ]
      [ -(1/2)1_30^T             15  ].

The characteristic-polynomial recurrence for the uniform join is chi31(z)=z(z-31/2)chi30(z-1/2)/(z-1/2). The inherited canonical factorization yields chi31(z)=z(z-31/2)^20(z-15)^4*f6(z-12), where f6(z)=chi_P7(z)/z is the nonzero canonical seven-vertex factor of the chosen mathematical model.

Independent exact-rational rank calculations show nullity(L26-(31/2)I)=20 and nullity(L26-15I)=4. The values are dimensionless eigenvalues of the chosen graph. They are not frequencies in hertz, an observed resonant material, or new particles.

## 5. Independent computational evidence and negative controls

The primary verifier src/verify_d26.py constructs every rational Gram entry directly, checks 465 exact labelled squared distances, asserts all 435 parent pairs remain unchanged, and computes both Euclidean and weighted Laplacian ranks using exact symbolic arithmetic. The checker includes four deliberately invalid candidates: a corrupted labelled pair, duplicated vertex, nonperpendicular tilt, and index permutation undetectable by an unlabeled histogram.

The independent verifier tests/independent_d26.py reconstructs all 465 labelled distances separately using the Python fractions.Fraction rational datatype and structural construction rules without importing the primary Gram implementation. It checks the frozen predecessor checksum, reproduces all distance labels, and independently recomputes Gram and Laplacian ranks. The two scripts print the same SHA-256 of results/certificate_d26.json.

Reproduction instructions: install Python and sympy==1.14.0, run python src/verify_d26.py followed by python tests/independent_d26.py. Continuous integration repeats both checks, requires their success and rejects any unexpected edit to the deterministic JSON certificate.

## 6. Limitations and prior work

Orthogonal Euclidean extensions and uniformly joined graph Laplacians are established mathematics. This note presents one particular exact, reproducible synthetic metric and one specified ideal graph model. It does not establish a universal new theorem, demonstrate a physical 26-dimensional spacetime, measure a vibration, identify a physical particle, or constitute independent human peer review.

## References

1. Gabriel St-Pierre, Polybrutus 25D: Exact 30-Vertex Euclidean Metric and Nineteenfold Ideal Graph Mode (2026), https://doi.org/10.5281/zenodo.24410435.
2. Topbrutus/polybrutus-d25, frozen scientific source v0.1.0, https://github.com/Topbrutus/polybrutus-d25/tree/v0.1.0.
3. Topbrutus/polybrutus canonical seven-point specification, https://github.com/Topbrutus/polybrutus/blob/main/SPECIFICATION.md.
