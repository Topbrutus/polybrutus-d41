# Polybrutus 37D: Exact 42-point Euclidean distance geometry and a thirty-onefold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**21 November 2026**  
**Exact mathematical case study / scientific preprint**  
**Material AI assistance:** OpenAI ChatGPT (Astra) helped exact rational calculations, two independently implemented checking scripts, scientific writing and publication. AI is not a scientific author.

## Abstract

We construct one exact Euclidean labelled squared-distance metric on 42 unit vertices by conservatively adding one perpendicular unit point to frozen Polybrutus 36D. All 820 predecessor labelled squared distances remain unchanged. The newly appended point adds 41 squared distances equal to 2 at unit radius, for 861 unordered labelled pairs. The centered Euclidean Gram is positive semidefinite of exact rank 37, establishing minimum affine realization dimension 37. Separately, a specifically chosen ideal complete unit-mass weighted graph with edge weights w_ij=1/d_ij² has Laplacian eigenvalue 21 of multiplicity 31 and eigenvalue 41/2 of multiplicity 4. Two separately implemented exact-rational verifiers reproduce all labels, original pair conservation, Gram ranks and chosen graph nullities. These claims do not establish physical extra spacetime dimensions or measured resonances.

## 1. Frozen scientific source and exact perpendicular extension

The immutable predecessor is Topbrutus/polybrutus-d36 at science tag v0.1.0, scientific commit b1f8be2cdcdcae0d9911d90381a7b1304c84abf0, DOI 10.5281/zenodo.24413083. Its original certificate contains 41 unit points and 820 exact labelled squared pair distances, with SHA-256 b189f74f7b38b378fa8f34011660817019c02601e5889378103fd9424207b331. Its scientific source and all earlier ancestors are archived without alteration under prior/d36/.

For each old unit point p_i, i=1,...,41, embed q_i=(p_i,0) and append the norm-one perpendicular basis point q42=e38. Every old scalar product and labelled pair squared distance is preserved. Orthogonality yields ||q42-q_i||²=||q42||²+||q_i||²-2<q42,q_i>=2. Therefore C(42,2)=861=820+41 exact unordered labelled squared-distance pairs.

## 2. Exact Euclidean Gram and minimum affine dimension

Seven canonical unit vertices form a rational Gram block of rank three with scalar products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vectors form an independent rank-four block with diagonal entry 1 and off-diagonal entry -1/4. Thirty additional mutually perpendicular unit axes are orthogonal to these two blocks. Hence the 42-by-42 rational Gram matrix G is positive semidefinite with unit diagonal and exact rank 3+4+30=37.

For the complete squared-distance matrix, define D_ii=0 and D_ij=2-2G_ij for distinct pairs. Let H=I-(1/42)11^T be the centering projection. The unique centered Euclidean Gram matrix is B=-(1/2)HDH=HGH, with exact rank 37. Two original canonical points are antipodal and place the origin in their affine hull. Since the complete labelled distances uniquely determine B, every Euclidean realization of this metric has affine dimension at least 37.

## 3. Exhaustive labelled distance invariants

At unit radius R=1, set the normalized squared-distance index k=4d²/R². There are eight exact classes: k=4 six times; k=7 twice; k=8 830 times; k=9 three times; k=10 ten times; k=12 six times; k=13 twice; k=16 twice. Their counts total 861. The exact sum of the 861 unordered squared pair distances is 6935/4.

Both verifiers compare all 820 predecessor squared-distance values by their original labels, not just a histogram. A deliberate label permutation can preserve the unlabeled histogram yet change which pair holds a value; the adversarial checks reject that case.

## 4. Separately stipulated ideal complete weighted graph

For one selected mathematical model, stipulate a complete weighted graph with unit scalar masses and symmetric weights w_ij=1/d_ij². The new perpendicular unit point joins every older point uniformly at edge weight 1/2. The weighted Laplacian has exact block form

L37 = [ L36+(1/2)I41     -(1/2)1_41 ]
      [ -(1/2)1_41^T         41/2 ].

Its characteristic polynomial satisfies the exact uniform-join recurrence chi42(z)=z(z-21)chi41(z-1/2)/(z-1/2). With original nonzero canonical seven-point graph factor f6(z)=chi_P7(z)/z, its inherited closed expression is chi42(z)=z(z-21)^31(z-41/2)^4*f6(z-35/2).

Two independent exact rational matrix rank calculations establish nullity(L37-21I)=31 and nullity(L37-(41/2)I)=4. These are dimensionless mathematical properties of the specified graph, not measured hertz, sound resonances, particle observations or additional physical spacetime dimensions.

## 5. Independent reproducibility and negative controls

The primary src/verify_d37.py derives an exact rational Gram and all 861 pairwise labelled squared distances, compares all 820 parent labels against the checksum-verified frozen predecessor, and checks centered Gram rank and selected graph eigenspaces. Four adverse inputs are rejected: one corrupted distance, duplicated point, nonperpendicular tilt and histogram-preserving label swapping.

Independently, tests/independent_d37.py computes every labelled squared distance using Python fractions.Fraction structural rules without importing the primary Gram generator. It checks the original parent SHA-256, exact 820 label conservation, 861 present pairs, centered Euclidean rank and graph multiplicities. Both programs print the same SHA-256 for results/certificate_d37.json. GitHub Actions reruns both checkers, rechecks the frozen parent integrity and verifies deterministic certificate output.

Reproduce with Python, sympy==1.14.0, python src/verify_d37.py and python tests/independent_d37.py. The frozen scientific parent and all ancestor records reside inside prior/d36/.

## 6. Scientific boundaries and references

Orthogonal Euclidean extensions and uniform graph joining are established mathematical techniques, here applied to one exactly reproducible synthetic metric and a separately stipulated graph. This is not a new universal theorem, evidence of 37 physical spacetime dimensions, experimental resonance, particle discovery or independent human peer-reviewed research.

References: Gabriel St-Pierre, Polybrutus 36D (2026), DOI 10.5281/zenodo.24413083; frozen source github.com/Topbrutus/polybrutus-d36 tag v0.1.0; original seven-point geometry github.com/Topbrutus/polybrutus.
