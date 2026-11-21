# Polybrutus 32D: Exact 37-point Euclidean distance geometry and a twenty-sixfold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**20 November 2026**  
**Mathematical case study / scientific preprint**  
**AI assistance:** OpenAI ChatGPT (Astra) materially assisted exact symbolic calculations, independent verification programs, writing and release. AI is not a scientific author.

## Abstract

We give one exact Euclidean labelled squared-distance metric on 37 unit vertices, conservatively extending the frozen Polybrutus 31D scientific configuration by one perpendicular unit point. All 630 predecessor labelled squared distances are unchanged and 36 new squared distances equal 2 at unit radius, giving 666 unordered labelled pairs. The exact centered Euclidean Gram matrix is positive semidefinite of rank 32, determining minimum affine Euclidean realization dimension 32. Separately, one stipulated equal-mass ideal complete weighted graph with w_ij=1/d_ij² has Laplacian eigenvalue 37/2 of exact multiplicity 26 and eigenvalue 18 of multiplicity 4. Two distinct programs verify the complete rational labelled metric, frozen parent conservation, centered Gram rank and these model eigenspaces. This construction is not evidence of physical spacetime dimensions or measured resonances.

## 1. Frozen source and orthogonal completion

The sole scientific predecessor is Topbrutus/polybrutus-d31 at scientific tag v0.1.0, source commit 0c1f546a6c5442bc83a5bdde506132cce20b6e05, DOI 10.5281/zenodo.24411160. Its certificate contains 36 unit points and 630 exact labelled squared distances; SHA-256 7c373c50d61ca2964eb832fc6fb8178590715e603f5fdada122575f624f60e8c. Its source and all prior ancestor records are preserved under prior/d31/.

Embed all old unit vertices p_i as q_i=(p_i,0), for i=1,...,36. Append the perpendicular norm-one basis point q37=e33. By orthogonality, ||q37-q_i||²=||q37||²+||q_i||²-2<q37,q_i>=2. Every old squared pair value is exactly conserved. The 666 unordered pairs comprise 630 prior labels and 36 new pairs with squared distance 2.

## 2. Exact Euclidean Gram and minimum affine dimension

The original seven canonical unit vertices form a rational rank-three Gram block with inner products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered unit simplex vertices form a rank-four block with diagonal inner products 1 and distinct-pair products -1/4. Twenty-five additional mutually perpendicular norm-one vectors contribute independent axes orthogonal to these two original blocks. Thus the 37-by-37 rational Gram matrix G is positive semidefinite, has unit diagonal and rank 3+4+25=32.

Define D_ii=0 and D_ij=2-2G_ij for every distinct labelled pair. For the centering projection H=I-(1/37)11^T, the centered squared-distance Gram matrix satisfies B=-(1/2)HDH=HGH. Its exact rank is 32. An original antipodal canonical pair shows that the origin belongs to the affine hull. Since a complete labelled Euclidean distance matrix determines B uniquely, the minimum affine Euclidean realization dimension of this metric is exactly 32.

## 3. Exhaustive labelled-distance invariants

For unit radius R=1, set k=4d²/R². There are eight exact distance families: k=4 occurs 6 times; k=7 occurs 2 times; k=8 occurs 635 times; k=9 occurs 3 times; k=10 occurs 10 times; k=12 occurs 6 times; k=13 occurs 2 times; and k=16 occurs 2 times. The sum is 666. The total sum of all unordered squared distances equals 5375/4.

Both verifying programs check all 630 frozen predecessor pair labels individually, not just a histogram. A deliberate label permutation preserves an unlabeled histogram but changes the original metric correspondence; the negative controls catch that discrepancy.

## 4. Selected ideal weighted graph and dimensionless modes

Independently of the Euclidean distance geometry, we stipulate a complete ideal graph with unit scalar masses and symmetric edge weights w_ij=1/d_ij². The new unit point connects uniformly to all 36 predecessor points with weight 1/2. Its exact weighted Laplacian has the block form:

L32 = [ L31 + (1/2)I36     -(1/2)1_36 ]
      [ -(1/2)1_36^T             18  ].

The characteristic-polynomial recurrence is chi37(z)=z(z-37/2)chi36(z-1/2)/(z-1/2). Using the inherited canonical seven-point factor f6(z)=chi_P7(z)/z, the decomposition is chi37(z)=z(z-37/2)^26(z-18)^4*f6(z-15).

The two independently implemented exact-rational Laplacian calculations establish nullity(L32-(37/2)I)=26 and nullity(L32-18I)=4. These modes are dimensionless mathematical values of this chosen graph, not sound frequencies, measured physical oscillations, particles, or evidence of 32 physical spacetime dimensions.

## 5. Exact independent verification and adversarial controls

The primary verifier src/verify_d32.py constructs all rational Gram matrix entries, all 666 labelled squared pairs, compares 630 original predecessor distances against the checksum-verified frozen D31 certificate, and checks centered Gram and Laplacian ranks. Four negative controls reject one corrupted pair, duplicate unit point, nonperpendicular tilt and a label swap that conserves the unlabeled histogram.

The independently implemented verifier tests/independent_d32.py reconstructs every labelled distance using Python fractions.Fraction structural rules, without importing the primary Gram generator. It also independently verifies the predecessor SHA-256 and all 630 labels, the centered Euclidean rank and the graph eigenspace dimensions. Both programs print the same deterministic SHA-256 for results/certificate_d32.json.

Install Python and sympy==1.14.0, then execute python src/verify_d32.py and python tests/independent_d32.py. GitHub Actions runs both scripts on the frozen source, checks the archived predecessor checksum and rejects changes to the generated exact certificate.

## 6. Scientific limitations and sources

Orthogonal Euclidean extension and uniform weighted graph joining are already established mathematics. The presented results describe one concrete exact synthetic metric and one deliberately selected graph. There is no new universal theorem, observed physical dimension, measured resonance, particle discovery or independent human peer review.

References: Gabriel St-Pierre, Polybrutus 31D (2026), DOI 10.5281/zenodo.24411160; frozen source github.com/Topbrutus/polybrutus-d31, v0.1.0; original seven-point geometry github.com/Topbrutus/polybrutus.
