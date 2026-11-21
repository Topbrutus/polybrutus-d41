# Polybrutus 24D: An exact 29-vertex Euclidean metric and an eighteenfold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**20 November 2026**  
**Mathematical case study / scientific preprint**  
**AI contribution disclosure:** OpenAI ChatGPT (Astra) assisted symbolic computation, independent verification code, manuscript preparation and publication. AI is not an author.

## Abstract

This note describes a particular 29-point unit-sphere Euclidean metric constructed by one exact perpendicular extension of the frozen Polybrutus 23D configuration. The extension conserves all 378 predecessor labelled squared distances, and introduces 28 new squared distances equal to 2, for a total of 406 unordered vertex pairs. Its centered Euclidean Gram matrix is positive semidefinite with exact rank 24, implying that any Euclidean realization of this labelled metric has minimum affine dimension 24. A deliberately stipulated equal-mass complete weighted graph with weights w_ij = 1/d_ij² has a Laplacian eigenvalue 29/2 of exact multiplicity 18 and a value 14 of multiplicity 4. The results are independently reproduced in exact rational arithmetic. This is not evidence of a new physical spacetime dimension or measured resonance.

## 1. Immutable predecessor

The scientific predecessor is Topbrutus/polybrutus-d23 tag v0.1.0, source commit 23f71c9b9932894c502fb67c1326271daa0692e7, Zenodo DOI 10.5281/zenodo.24410035. The labelled D23 certificate contains 28 radius-one unit vertices and 378 squared pair distances. The original full sources and original parent JSON are included as an unchanged reference under prior/d23/.

The old 28 vertices p1 through p28 are embedded as q_i=(p_i,0) in one extra storage slot, and a unit vector q29 is appended orthogonally. Every previous pair keeps exactly the same squared distance, while ||q29-q_i||²=||q29||²+||q_i||²-2<q29,q_i>=2. This creates 28 new labelled pair distances, making 378+28=406=C(29,2) pairs.

## 2. Exact Euclidean Gram construction and minimal affine dimension

The seven canonical points form the original three-dimensional core with inner products computed as 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex vectors contribute a rank-four block, with inner products 1 on the diagonal and -1/4 for any two distinct simplex vertices. Seventeen mutually perpendicular unit axes occupy additional orthogonal slots. The Gram matrix G is block-diagonal in this decomposition, is positive semidefinite, and has rank 3+4+17=24. Every diagonal Gram entry is 1.

The full exact squared-distance matrix D has D_ii=0 and D_ij=2-2G_ij. The centered Euclidean squared-distance Gram matrix satisfies B=-(1/2) H D H = H G H for H=I-11^T/29. The first and fifth core points are antipodal, so the origin lies in the affine hull; rational symbolic rank gives rank(B)=24. Since complete labelled squared distances uniquely determine the centered Gram matrix, no lower-dimensional Euclidean realization of this same full labelled metric exists.

## 3. Distance invariants

For R=1 the normalized class counts of k=4d²/R² are: 4:6, 7:2, 8:375, 9:3, 10:10, 12:6, 13:2, 16:2. These eight classes total exactly 406 unordered labelled pairs. Their total squared-distance sum is 3295/4.

Two separately written programs check all 406 labelled distances against the exact certificate and separately compare all 378 old labelled pairs against the D23 immutable predecessor certificate. Checking the labelled metric is stronger than checking its histogram: swapping labels leaves the multiset of distances unchanged but fails the complete labelled comparison.

## 4. Deliberately chosen weighted graph spectral model

A complete weighted graph is defined with scalar masses set to 1 and edge weights w_ij=1/d_ij²; those couplings are mathematical choices, not laboratory measurements. The perpendicular new vertex attaches to each predecessor with weight 1/2. Thus its weighted Laplacian is:

L24 = [ L23 + (1/2)I28     -(1/2)1_28 ]
      [ -(1/2)1_28^T          14      ].

The characteristic polynomial obeys the uniform-join recurrence chi29(z)=z(z-29/2)chi28(z-1/2)/(z-1/2). The parent factorization yields chi29(z)=z(z-29/2)^18(z-14)^4*f6(z-11), where f6(z)=chi_P7(z)/z is the nonzero canonical seven-point graph factor.

Independent exact-rational Laplacian rank computations establish nullity(L24-(29/2)I)=18 and nullity(L24-14I)=4. They are dimensionless eigenvalues of the stipulated mathematical model; they are not measured oscillations or discovered particles.

## 5. Independent reproduction and adversarial controls

The first script src/verify_d24.py constructs the full Gram matrix rationally, establishes all 406 squared distances and the exact centered Gram rank, verifies the 378 old labelled distances, and checks the two spectral multiplicities by exact ranks. Four false configurations are rejected: a corrupted labelled pair, a duplicated vertex, a tilted non-perpendicular vertex, and an index relabelling that preserves the histogram.

The second program tests/independent_d24.py computes all 406 distances independently from Fraction-based rational structural rules rather than importing the primary Gram generator. It verifies the old certificate separately, centered Gram rank and Laplacian ranks, and prints the same SHA-256 certificate value.

Reproduce using Python and sympy==1.14.0: python src/verify_d24.py, then python tests/independent_d24.py. Verify that both output the same CERTIFICATE_SHA256 and that the generated certificate has no git diff. GitHub Actions repeats these checks for source commits and the frozen scientific tag.

## 6. Limits and prior work

Orthogonal Euclidean extension and the uniform-join recurrence are known mathematical techniques. This note documents one exact reproducible labelled metric and chosen ideal graph. It does not present a new universal mathematical theorem, measurement of physical frequency, evidence for an extra physical spacetime dimension or independent human peer review.

## References

1. Gabriel St-Pierre, Polybrutus 23D: Exact 28-Point Euclidean Metric and Seventeenfold Ideal Graph Mode. Zenodo, 2026. DOI: https://doi.org/10.5281/zenodo.24410035.
2. Topbrutus/polybrutus-d23, immutable scientific source v0.1.0, https://github.com/Topbrutus/polybrutus-d23/tree/v0.1.0.
3. Topbrutus/polybrutus canonical seven-point specification, https://github.com/Topbrutus/polybrutus/blob/main/SPECIFICATION.md.
