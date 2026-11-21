# Polybrutus 40D: Exact 45-point Euclidean distance geometry and a thirty-fourfold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**21 November 2026**  
**Exact mathematical case study / scientific preprint**  
**Material AI assistance:** OpenAI ChatGPT (Astra) assisted exact rational mathematics, two distinctly implemented verifiers, scientific writing and controlled publication. AI is not a scientific author.

## Abstract

We construct an exact Euclidean labelled squared-distance metric on 45 norm-one points by conservatively adding one perpendicular unit point to frozen Polybrutus 39D. The complete 946 predecessor labelled squared distances remain unchanged. A new perpendicular unit point contributes 44 additional squared distances equal to 2 at unit radius, giving 990 unordered labelled pairs. The exact rational centered Euclidean Gram is positive semidefinite of rank 40, proving minimum affine realization dimension 40 for this particular metric. For a separately stipulated ideal complete unit-mass weighted graph with w_ij=1/d_ij², the weighted Laplacian has eigenvalue 45/2 of multiplicity 34 and eigenvalue 22 of multiplicity 4. Two separately implemented exact-rational verification programs reproduce all 990 labels, predecessor conservation, Euclidean rank and graph eigenspaces. These exact mathematical results do not establish physical extra spacetime dimensions, observed particles or measured frequencies.

## 1. Frozen scientific parent and orthogonal extension

The sole frozen scientific parent is Topbrutus/polybrutus-d39 at tag v0.1.0, source commit a43cd3c67017fb6162005fd521c156cb78795beb, DOI 10.5281/zenodo.24413712. The parent certificate covers 44 unit points and 946 labelled unordered squared pair distances, SHA-256 8cd6731c294537d74c04a3501fe0da76a2cee6d496a2d2d5dfcfb0a26bb9304f. Its complete source and all inherited ancestors remain unmodified in prior/d39/.

Map every previous unit point p_i, i=1,...,44, to q_i=(p_i,0) in a new ambient coordinate, and adjoin a new perpendicular norm-one vector q45=e41. All original dot products and squared distances remain exactly equal by embedding. Each new squared pair distance is ||q45-q_i||²=||q45||²+||q_i||²-2<q45,q_i>=2. We obtain C(45,2)=990=946+44 exact unordered labelled squared distances.

## 2. Exact rational Euclidean Gram and minimal affine rank

The seven original canonical unit vertices form a rational rank-three Gram block with products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vertices form an independent rank-four Gram block, with diagonal products one and off-diagonal products -1/4. Thirty-three more mutually perpendicular unit coordinate vectors are orthogonal to both preceding blocks. The full 45-by-45 rational Gram matrix G is positive semidefinite, has unit diagonal, and rank 3+4+33=40.

For the complete labelled squared-distance matrix D, define D_ii=0 and D_ij=2-2G_ij for each distinct pair. Let H=I-(1/45)11^T. Its uniquely determined centered Euclidean Gram B=-(1/2)HDH=HGH has exact rank 40. Two original canonical points are antipodal, placing the origin in their affine span. Since every Euclidean realization of these same complete labelled distances has the same centered Gram rank, the minimum affine Euclidean dimension is 40.

## 3. Exhaustive exact labelled-distance arithmetic

For R=1 define normalized squared-distance index k=4d²/R². The eight exact classes are: k=4 six pairs; k=7 two pairs; k=8 959 pairs; k=9 three pairs; k=10 ten pairs; k=12 six pairs; k=13 two pairs; and k=16 two pairs. These total 990 labelled unordered pairs. The sum of all squared pair distances is 7967/4.

The two programs check each of 946 frozen predecessor values with its original pair label and verified parent SHA-256, not just the histogram. A deliberately permuted labelling leaves the histogram unchanged but breaks the pair-label conservation and is rejected.

## 4. Separately stipulated ideal unit-mass weighted graph

Independently stipulate a complete weighted graph with scalar unit masses and symmetric edge weights w_ij=1/d_ij². The new perpendicular point joins all 44 predecessor vertices uniformly with weight 1/2. The weighted Laplacian satisfies the exact block identity

L40 = [ L39 + I44/2       -1_44/2 ]
      [ -1_44^T               22 ].

Its characteristic polynomial satisfies chi45(z)=z(z-45/2)chi44(z-1/2)/(z-1/2). Using f6(z)=chi_P7(z)/z for the inherited nonzero seven-point canonical factor, its model-specific closed expression is chi45(z)=z(z-45/2)^34(z-22)^4*f6(z-19).

Two separately implemented exact rational rank checks establish nullity(L40-(45/2)I)=34 and nullity(L40-22I)=4. These are dimensionless eigenvalues of a freely stipulated mathematical graph, not sound frequencies, measured resonances, particles or additional physical dimensions.

## 5. Independent reproducibility and negative controls

The primary program src/verify_d40.py constructs all exact rational Gram entries and 990 labelled squared distances, checksum-verifies the frozen predecessor, checks all 946 inherited labels, centered Euclidean rank and selected graph eigenspaces. Four adversarial negative cases are rejected: corrupt distance, duplicate point, nonorthogonal tilt and label permutation retaining the histogram.

The separately implemented tests/independent_d40.py program calculates each labelled squared distance using Python fractions.Fraction structural rules and never imports the primary Gram generator. It validates the original predecessor SHA-256, all 946 parent labels, all 990 present labels and exact Euclidean and graph ranks. Both programs generate the same deterministic SHA-256 for results/certificate_d40.json. GitHub Actions runs both on the published source and checks the parent certificate hash and output determinism.

Reproduce with Python and sympy==1.14.0, running python src/verify_d40.py and python tests/independent_d40.py. The exact scientific parent and inherited ancestors are preserved in prior/d39/.

## 6. Scientific scope and references

Perpendicular extensions of Euclidean point configurations and uniform graph joins are established mathematical techniques. This is one reproducible synthetic metric and one independently stipulated ideal graph, not a new universal theorem, experimental observation of forty physical spacetime dimensions, measured frequencies, particles or independent human peer review.

References: Gabriel St-Pierre, Polybrutus 39D (2026), DOI 10.5281/zenodo.24413712; frozen source github.com/Topbrutus/polybrutus-d39 v0.1.0; original seven-point geometry github.com/Topbrutus/polybrutus.
