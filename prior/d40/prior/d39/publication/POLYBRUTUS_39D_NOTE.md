# Polybrutus 39D: Exact 44-point Euclidean distance geometry and a thirty-threefold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**21 November 2026**  
**Exact mathematical case study / scientific preprint**  
**AI disclosure:** OpenAI ChatGPT (Astra) materially assisted exact rational calculations, two separately implemented verifiers, manuscript writing and controlled publication. AI is not a scientific author.

## Abstract

We construct a specific exact Euclidean labelled squared-distance metric on 44 norm-one points by conservatively adjoining one perpendicular unit point to the frozen Polybrutus D38 source. All 903 predecessor labelled squared distances are unchanged, and 43 new distances squared equal 2 at unit radius, giving 946 unordered labelled pairs. The exact rational centered Euclidean Gram is positive semidefinite of rank 39, establishing minimum affine realization dimension 39 for this metric. Separately, a deliberately stipulated unit-mass complete weighted graph with edges w_ij=1/d_ij² has Laplacian eigenvalue 22 of multiplicity 33 and eigenvalue 43/2 of multiplicity 4. Two separately implemented exact-rational programs reproduce every labelled pair, all frozen predecessor values, Euclidean ranks and selected graph eigenspaces. These are mathematical properties of a synthetic construction, not extra physical spacetime dimensions or measured resonances.

## 1. Frozen source and conservative perpendicular extension

The sole scientific predecessor is Topbrutus/polybrutus-d38 at tag v0.1.0, frozen scientific commit ac8fd112fcd314092f106b210162ba1cc9879e9b, DOI 10.5281/zenodo.24413606. Its exact original certificate records 43 unit points and 903 labelled unordered squared distances, SHA-256 3be0e4aaa1487e82f35a97a3521abe70b05992888c550c8fdd95b5d50f306dbc. All frozen parent source files and older ancestors are archived unchanged under prior/d38/.

Embed each predecessor unit point p_i, i=1,...,43, as q_i=(p_i,0) in a new ambient coordinate. Adjoin an orthogonal norm-one point q44=e40. Every old dot product and every one of 903 parent labelled distances are unchanged. The new point has squared distance ||q44-q_i||²=||q44||²+||q_i||²-2<q44,q_i>=2 for all old indices. There are therefore 946=C(44,2) exact unordered labelled squared pairs: 903 inherited labels and 43 new pairs of squared distance 2.

## 2. Exact rational Euclidean Gram and dimension minimality

Seven canonical norm-one points form a rational rank-three Gram block whose inner products are 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vectors form an independent rank-four block with diagonal inner products one and off-diagonal products -1/4. Thirty-two more pairwise perpendicular unit coordinate vectors are orthogonal to both blocks. Thus the 44-by-44 exact rational Gram matrix G is positive semidefinite with unit diagonal and rank 3+4+32=39.

Define full squared-distance matrix D_ii=0 and D_ij=2-2G_ij for each distinct pair. Let H=I-(1/44)11^T. The uniquely determined centered Euclidean Gram is B=-(1/2)HDH=HGH, with exact rank 39. Two original canonical points are antipodal, so the origin lies in the affine hull. As B is uniquely determined by the complete labelled distances, every Euclidean realization of this particular metric requires affine dimension at least 39.

## 3. Complete labelled distance arithmetic

At unit radius R=1 define normalized squared-distance index k=4d²/R². Exactly eight classes occur: k=4 six times; k=7 twice; k=8 915 times; k=9 three times; k=10 ten times; k=12 six times; k=13 twice; and k=16 twice. Counts sum to 946. The exact sum of all unordered squared pair distances is 7615/4.

Every one of the 903 predecessor labelled squared distances is separately checked against the SHA-256-verified frozen certificate. Matching a histogram alone would not demonstrate conservation, because labels can be permuted without changing an unlabeled distribution; an adversarial relabelling is rejected.

## 4. Separately stipulated complete weighted graph and spectral modes

Independently choose an ideal complete graph with unit scalar masses and symmetric edge weights w_ij=1/d_ij². The appended perpendicular point joins the 43 older vertices at edge weight 1/2. Consequently its Laplacian has the exact uniform-join block identity

L39 = [ L38+(1/2)I43      -(1/2)1_43 ]
      [ -(1/2)1_43^T          43/2 ].

The characteristic polynomial recurrence is chi44(z)=z(z-22)chi43(z-1/2)/(z-1/2). With the inherited original seven-point nonzero model factor f6(z)=chi_P7(z)/z, its closed form is chi44(z)=z(z-22)^33(z-43/2)^4*f6(z-37/2).

Two independently implemented exact-rational matrix-rank checks establish nullity(L39-22I)=33 and nullity(L39-(43/2)I)=4. These are dimensionless eigenmodes of a stipulated mathematical graph and do not describe experimental hertz, observed resonance, particle discoveries or extra physical spacetime dimensions.

## 5. Independent reproducibility and negative controls

The primary program src/verify_d39.py constructs exact Gram entries and all 946 labelled distances, checksum-verifies the original frozen D38 certificate, checks all 903 inherited labels and exact Euclidean and ideal graph ranks. Four false configurations are rejected: a damaged pair distance, duplicated vertex, non-orthogonal tilt and histogram-preserving label swap.

The independently implemented tests/independent_d39.py program reconstructs every labelled squared distance with Python fractions.Fraction structural rules without importing the primary Gram generator. It verifies the frozen parent SHA-256, all 903 predecessor labels, all 946 present labels, the exact centered Euclidean rank and stipulated graph nullities. Both programs print the same deterministic SHA-256 of results/certificate_d39.json. Public GitHub Actions runs both on the frozen source, verifies predecessor SHA-256 and fails if certificate bytes differ.

To reproduce: install Python and sympy==1.14.0; run python src/verify_d39.py and python tests/independent_d39.py. The frozen predecessor and its older source history are included under prior/d38/.

## 6. Scientific limits and references

Orthogonal extension of Euclidean metrics and uniform joining of complete weighted graphs are established mathematical techniques. This example does not constitute a novel universal theorem, independent human peer review, observed additional physical spacetime dimensions, particles or experimental resonances.

References: Gabriel St-Pierre, Polybrutus D38 (2026), DOI 10.5281/zenodo.24413606; frozen source github.com/Topbrutus/polybrutus-d38 tag v0.1.0; original seven-point geometry github.com/Topbrutus/polybrutus.
