# Polybrutus 41D: Exact 46-point Euclidean distance geometry and a thirty-fivefold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**21 November 2026**  
**Exact mathematical case study / scientific preprint**  
**AI disclosure:** OpenAI ChatGPT (Astra) materially assisted exact rational computation, two separate verifier implementations, manuscript writing and controlled publication. AI is not a scientific author.

## Abstract

We construct one exact Euclidean labelled squared-distance metric on 46 norm-one vertices by conservatively adjoining a perpendicular unit vector to the frozen Polybrutus 40D scientific source. All 990 predecessor labelled squared distances are preserved, while the appended point adds 45 new squared distances equal to 2 at unit radius. There are thus 1035 unordered labelled pairs. The complete exact rational centered Euclidean Gram is positive semidefinite of rank 41, proving minimum affine realization dimension 41 for this particular metric. For a separately stipulated ideal complete unit-mass weighted graph with edges w_ij=1/d_ij², the graph Laplacian has eigenvalue 23 of multiplicity 35 and 45/2 of multiplicity 4. Two distinct exact-rational verifier implementations check all pair labels, parent conservation, Euclidean rank and ideal graph eigenspaces. No physical extra spacetime dimensions or experimentally measured resonances are inferred.

## 1. Immutable parent and conservative perpendicular extension

The frozen scientific parent is Topbrutus/polybrutus-d40 at tag v0.1.0, scientific commit 2a9fa8265b709d6d971afb6c63b7b024011e67a1, DOI 10.5281/zenodo.24414238. The parent exact certificate lists 45 norm-one points and 990 unordered labelled squared pair distances, SHA-256 af361ae471397b6d91beed87fe65c268a7139c02fbc0187400152fb8b202ffca. Its complete frozen scientific source and older ancestors remain unmodified inside prior/d40/.

Map every old unit point p_i, i=1,...,45, to q_i=(p_i,0) in an additional coordinate. Append a new perpendicular norm-one vertex q46=e42. Every old dot product and pair squared distance is unchanged. By orthogonality, each new squared distance is ||q46-q_i||²=||q46||²+||q_i||²-2<q46,q_i>=2. Consequently C(46,2)=1035=990+45 unordered labelled pair squared distances are present.

## 2. Exact rational Gram and minimal affine Euclidean dimension

Seven canonical norm-one points have a rational rank-three Gram block with inner products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex norm-one points form an independent rank-four block, diagonal one and off-diagonal products -1/4. Thirty-four further mutually orthogonal unit coordinate vectors are orthogonal to these two blocks. The full 46-by-46 exact rational Gram G is positive semidefinite with unit diagonal and exact rank 3+4+34=41.

Let D_ii=0 and D_ij=2-2G_ij for every distinct pair, and let the centering projector be H=I-(1/46)11^T. The unique centered Euclidean Gram B=-(1/2)HDH=HGH has exact rational rank 41. Two original canonical unit points are antipodal, so their affine span contains the origin. Since every realization of the same complete labelled distances determines the same centered Gram rank, the minimum affine Euclidean dimension is 41.

## 3. Complete exact distance-label arithmetic

At unit radius define k=4d²/R². The exact distance classes are k=4 on six pairs; k=7 on two; k=8 on 1004; k=9 on three; k=10 on ten; k=12 on six; k=13 on two; and k=16 on two. Their counts sum to 1035. The exact sum of all squared distances is 8327/4.

Both verifiers compare every one of the 990 frozen predecessor pair distances with its original label, not just the distribution of values. A deliberate relabelling can leave the histogram invariant while changing which ordered label receives a value; the negative controls detect it.

## 4. Separately stipulated ideal weighted complete graph

Independently stipulate a complete weighted graph with scalar unit masses and symmetric edge weights w_ij=1/d_ij². The new perpendicular point joins every one of the 45 old points with uniform edge weight 1/2. The exact Laplacian block identity is

L41 = [ L40+(1/2)I45        -(1/2)1_45 ]
      [ -(1/2)1_45^T            45/2 ].

The characteristic polynomial satisfies chi46(z)=z(z-23)chi45(z-1/2)/(z-1/2). With f6(z)=chi_P7(z)/z the original canonical seven-point nonzero graph factor, the exact model-specific expression is chi46(z)=z(z-23)^35(z-45/2)^4*f6(z-39/2).

Two separate exact rational matrix rank computations give nullity(L41-23I)=35 and nullity(L41-(45/2)I)=4. These are dimensionless eigenvalues of an ideal chosen mathematical graph, not experimental measurements of hertz, particle frequencies, resonances or physical spacetime dimensions.

## 5. Distinct verification implementations and adversarial checks

The primary src/verify_d41.py program constructs exact rational Gram entries and all 1035 labelled squared distances, verifies the parent certificate's original SHA-256, checks the 990 inherited labels, centered Euclidean rank and stipulated graph eigenspaces. Four negative cases are rejected: one corrupted squared pair distance, duplicated point, non-orthogonal tilt, and a histogram-preserving label permutation.

The independently implemented tests/independent_d41.py program reconstructs all labelled squared distances using Python fractions.Fraction structural rules, without importing the primary Gram generator. It independently checks the frozen parent SHA-256, all 990 parent labels, 1035 present pairs, exact centered Euclidean rank and selected graph eigenspace multiplicities. Both scripts emit the same deterministic SHA-256 for results/certificate_d41.json, and a public GitHub Actions workflow reruns them and checks the committed bytes remain unchanged.

Reproduce by installing Python and sympy==1.14.0, then running python src/verify_d41.py and python tests/independent_d41.py. The entire parent and ancestral science source is preserved in prior/d40/.

## 6. Scientific boundaries and references

Perpendicular Euclidean extensions and uniform joins of weighted complete graphs are established mathematical methods. This preprint reports one synthetic, reproducible labelled Euclidean metric and one separately stipulated ideal graph, not a new universal theorem or experimental confirmation of additional physical dimensions. No observed particles, measured resonances, or independent human peer review are claimed.

References: Gabriel St-Pierre, Polybrutus 40D (2026), DOI 10.5281/zenodo.24414238; frozen scientific source github.com/Topbrutus/polybrutus-d40 v0.1.0; original seven-point geometry github.com/Topbrutus/polybrutus.
