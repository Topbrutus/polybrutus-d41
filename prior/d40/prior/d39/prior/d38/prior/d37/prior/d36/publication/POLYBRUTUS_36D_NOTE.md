# Polybrutus 36D: Exact 41-point Euclidean distance geometry and a thirtyfold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**21 November 2026**  
**Exact mathematical case study / scientific preprint**  
**AI disclosure:** OpenAI ChatGPT (Astra) materially assisted exact-rational calculations, two independently implemented verifiers, scientific text and controlled release. AI is not a scientific author.

## Abstract

We construct a specific exact Euclidean labelled squared-distance metric on 41 unit vertices by conservatively adding one perpendicular unit point to frozen Polybrutus 35D. All 780 frozen predecessor labelled squared distances remain unchanged, while 40 new squared distances equal 2 at unit radius, yielding 820 unordered labelled pairs. The exact rational centered Euclidean Gram is positive semidefinite of rank 36, establishing minimum affine realization dimension 36 for this complete labelled metric. Separately, a stipulated ideal equal-mass complete weighted graph with w_ij=1/d_ij² has Laplacian eigenvalue 41/2 of multiplicity 30 and eigenvalue 20 of multiplicity 4. Two separately implemented exact-rational checkers confirm all pair labels, parent conservation, centered Gram rank and ideal graph nullities. No physical extra dimension, measured resonance, particle discovery or universal new theorem is inferred.

## 1. Immutable predecessor and orthogonal extension

The only predecessor is Topbrutus/polybrutus-d35 at scientific tag v0.1.0, frozen source commit 72989a1fd2a7fb6e28522a0f2dcbe168924ca78d, DOI 10.5281/zenodo.24412435. Its original certificate lists 40 unit points and 780 exact labelled squared pair distances, SHA-256 7390c10f46132faeca0ccd8b52050d5ed899fbd4a078b9dd02407d922a4e4f22. Its complete source tag and inherited ancestors are preserved under prior/d35/.

For the original points p_i, i=1,...,40, use q_i=(p_i,0) in an extra coordinate. Append one perpendicular norm-one point q41=e37. All old dot products and labelled squared distances are preserved. Orthogonality gives ||q41-q_i||²=||q41||²+||q_i||²-2<q41,q_i>=2. Therefore C(41,2)=820=780+40 exact unordered labelled squared-distance pairs.

## 2. Exact rational Gram construction and minimal realization

Seven original canonical unit vertices form a rational Gram block of rank three, with products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vectors form a rank-four Gram block with diagonal products 1 and off-diagonal products -1/4. Twenty-nine more mutually perpendicular unit axes are independent of these blocks. Hence the exact rational 41-by-41 Gram matrix G is positive semidefinite with unit diagonal and rank 3+4+29=36.

Set the complete labelled squared-distance matrix D_ii=0 and D_ij=2-2G_ij for distinct points. Let H=I-(1/41)11^T. Its unique centered Euclidean Gram is B=-(1/2)HDH=HGH, with exact rank 36. The original first and fifth canonical points are antipodal, putting the origin in their affine hull. The complete labelled metric therefore requires minimum Euclidean affine dimension 36 in every realization.

## 3. Exhaustive label-based invariants

At unit radius, normalize squared distances by k=4d²/R². Exactly eight classes occur: k=4 six times; k=7 twice; k=8 789 times; k=9 three times; k=10 ten times; k=12 six times; k=13 twice; and k=16 twice. The counts sum to 820. The exact sum of all unordered squared pair distances for R=1 is 6607/4.

Both verifying programs compare the 780 predecessor pair values by their exact labels, not only by the histogram. A deliberate label permutation preserves the unlabeled histogram but breaks the frozen metric correspondence and is rejected by adversarial checks.

## 4. Stipulated ideal weighted complete graph

Independently choose a mathematical model with unit scalar masses and complete edge weights w_ij=1/d_ij². The new perpendicular unit vertex joins all 40 original vertices at weight 1/2, yielding the exact block Laplacian

L36 = [ L35+(1/2)I40      -(1/2)1_40 ]
      [ -(1/2)1_40^T             20 ].

The characteristic polynomial obeys chi41(z)=z(z-41/2)chi40(z-1/2)/(z-1/2). Using f6(z)=chi_P7(z)/z for the inherited nonzero canonical seven-point polynomial factor, the exact model expression is chi41(z)=z(z-41/2)^30(z-20)^4*f6(z-17).

Two independent exact-rational rank checks establish nullity(L36-(41/2)I)=30 and nullity(L36-20I)=4. These are dimensionless values of the freely chosen complete-graph Laplacian, not measured frequencies, sound, particles or physical spacetime dimensions.

## 5. Two exact checkers and adversarial controls

The primary src/verify_d36.py constructs rational Gram entries, derives all 820 squared pair labels, confirms each of 780 old parent values using the SHA-256-checked frozen certificate, and verifies Euclidean and ideal graph ranks. Four negative controls reject one corrupted pair, repeated point, non-orthogonal tilt and a label permutation with unchanged histogram.

The separately implemented tests/independent_d36.py reconstructs each exact pair distance with Python fractions.Fraction rules, without importing the primary Gram generator. It independently checks frozen predecessor SHA-256 and all 780 old labels, all 820 present labels, centered rank and selected-model eigenvalue multiplicities. Both print the identical deterministic SHA-256 of results/certificate_d36.json.

Reproduction: install Python and sympy==1.14.0, execute python src/verify_d36.py and python tests/independent_d36.py. GitHub Actions reruns both verifiers and asserts the committed certificate remains unchanged.

## 6. Scientific scope and references

Orthogonal embedding of Euclidean metrics and uniform weighted graph joins are established techniques. The result is one reproducible exact synthetic metric and one stipulated graph model; it is not a new universal theorem, a physical 36-dimensional spacetime observation, a measured resonance or a particle discovery. Independent human peer review has not occurred.

References: Gabriel St-Pierre, Polybrutus 35D (2026), DOI 10.5281/zenodo.24412435; frozen scientific source github.com/Topbrutus/polybrutus-d35, tag v0.1.0; original canonical seven-point geometry github.com/Topbrutus/polybrutus.
