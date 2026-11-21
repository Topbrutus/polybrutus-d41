# Polybrutus 23D: exact 28-point Euclidean metric and ideal weighted-graph spectrum

**Gabriel St-Pierre** - Independent researcher  
**20 November 2026**  
**Type:** Mathematical case study / preprint  
**AI disclosure:** ChatGPT (Astra) assisted computation, verification, manuscript preparation, and release operations. The research author and scientific responsibility are human.

## Abstract

We describe a labelled Euclidean metric on 28 unit vertices obtained as one exact perpendicular extension of the previously released Polybrutus 22D configuration. All 351 predecessor squared distances are preserved label by label. The appended vertex contributes exactly 27 new squared distances equal to 2 (at unit radius), yielding 378 labelled pairs. The exact centered Gram matrix is positive semidefinite of rank 23; this proves that the minimum Euclidean affine realization dimension of this particular metric is 23. For a separately chosen complete weighted graph with weights w(i,j)=1/d(i,j)^2, two rational calculations independently verify eigenvalue 14 of multiplicity 17 and eigenvalue 27/2 of multiplicity 4. These are synthetic mathematical constructions rather than measurements or claims of new physical dimensions.

## 1. Immutable predecessor and objective

The predecessor is Topbrutus/polybrutus-d22 at the frozen scientific tag v0.1.0, source commit 4d4e80a7156a118f8f2bebf63aca76ac5f0f06eb (DOI 10.5281/zenodo.24409636). Its certificate lists 27 unit vertices and 351 labelled squared distances. The predecessor certificate and its verification sources are archived unaltered in the present repository under prior/d22/. We do not modify any previous source or publication.

Let p1,...,p27 be the original unit points in 23 storage-coordinate slots. Embed them as q_i=(p_i,0) in 24 storage slots. Append the new perpendicular unit basis point q28=e24. Therefore, for every i in 1,...,27:

  ||q28 - q_i||^2 = ||q28||^2 + ||q_i||^2 - 2<q28,q_i> = 2.

Every old labelled squared distance is identical to its predecessor value, and 351+27=378=C(28,2).

## 2. Exact rational Gram construction

The seven canonical points have scalar products 3*x_i*x_j+y_i*y_j+z_i*z_j in their three-dimensional block. Five centered regular-simplex points have pairwise products 1 on the diagonal and -1/4 off diagonal; this block has rank 4. Sixteen further mutually perpendicular unit axes contribute rank 16. All cross-block products vanish.

Consequently, the explicit 28 by 28 Gram matrix G is rational, positive semidefinite and of rank 3+4+16=23. Its diagonal is identically 1, so the squared-distance matrix is D_ij=2-2G_ij for distinct i,j, with zeros on the diagonal. The canonical first and fifth vertices are antipodal, placing the origin in the affine hull and showing affine rank 23.

For H=I-(1/28)11^T, the exact centered Euclidean squared-distance Gram matrix is

  B=-(1/2) H D H = H G H.

The two independently implemented exact-rational verifiers establish rank(B)=23. Because a complete labelled Euclidean squared-distance matrix uniquely determines B, no Euclidean embedding of the same exact labelled metric has smaller affine dimension.

## 3. Distance invariants and traceability

For radius R=1, the normalized classes k=4d^2/R^2 and their exact labelled counts are:
k=4:6, k=7:2, k=8:347, k=9:3, k=10:10, k=12:6, k=13:2, k=16:2.
The counts sum to 378, and the sum over all unordered pairs of squared distances equals 3071/4.

The predecessor's 351 labelled distances are verified exactly against its immutable certificate, not merely compared through a histogram. A label permutation preserves the histogram but fails full labelled-distance comparison.

## 4. Chosen graph spectral model

In a deliberately stipulated mathematical model, give each vertex a unit scalar mass and each distinct pair an edge weight w_ij=1/D_ij. This defines a complete weighted graph, not a physically measured oscillator.

The newly added vertex is coupled uniformly with weight 1/2 to each of the old 27 vertices. If L22 is the old model Laplacian, the new Laplacian is the block matrix

  L23 = [ L22 + (1/2)I27   -(1/2)1_27 ]
        [ -(1/2)1_27^T       27/2    ].

With chi_n(z) the characteristic polynomial, the uniform-join recursion is

  chi28(z) = z(z-14) chi27(z-1/2) / (z-1/2).

Equivalently the form inherited from the canonical 7-vertex graph factor f6 is

  chi28(z)=z(z-14)^17(z-27/2)^4 f6(z-21/2).

Independent rational ranks show nullity(L23-14I)=17 and nullity(L23-(27/2)I)=4. The two values are dimensionless eigenvalues of this explicitly chosen ideal graph, not frequencies in hertz or empirically observed particles.

## 5. Independent computational verification

The primary program src/verify_d23.py uses rational Gram entries and exact linear algebra. The independent program tests/independent_d23.py reconstructs every labelled squared distance separately from fractions.Fraction rules and never imports the primary coordinate or Gram generator.

Both verify all 378 labelled distances, preservation of the 351 parent pairs, centered Gram rank 23, the two stated graph eigenvalue multiplicities, and the total squared-distance sum. The primary verification additionally rejects four purposely incorrect variants: a corrupted pair, a duplicated new vertex, a tilted new vertex, and a label swap hidden from an unlabeled histogram.

The scripts emit and independently check one deterministic JSON certificate, results/certificate_d23.json, identified by its SHA-256. Continuous integration reruns both and checks the certificate for a clean reproducible diff.

## 6. Scientific limitations and reproducibility

Orthogonal extension and the uniform-join graph recurrence are standard mathematical techniques. The present work is a specific exactly reproducible synthetic case study, not a newly established universal theorem. Euclidean affine dimension 23 is not evidence for a 23-dimensional physical spacetime. No measured resonance, experimental apparatus, or independent human peer review is claimed.

To reproduce: install Python and sympy==1.14.0, run python src/verify_d23.py, then python tests/independent_d23.py. Check that the printed SHA-256 hashes of the resulting certificate agree. The parent files in prior/d22/ are frozen copies of the D22 scientific source.

## References

1. Gabriel St-Pierre, Polybrutus 22D - Exact 27-point Euclidean metric, Zenodo (2026), https://doi.org/10.5281/zenodo.24409636.
2. Topbrutus/polybrutus-d22 scientific source tag v0.1.0, https://github.com/Topbrutus/polybrutus-d22/tree/v0.1.0.
3. Topbrutus/polybrutus canonical specification, https://github.com/Topbrutus/polybrutus/blob/main/SPECIFICATION.md.
