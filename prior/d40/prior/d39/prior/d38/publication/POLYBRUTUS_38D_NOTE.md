# Polybrutus 38D: Exact 43-point Euclidean distance geometry and a thirty-twofold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**21 November 2026**  
**Exact mathematical case study / scientific preprint**  
**Material AI assistance:** OpenAI ChatGPT (Astra) assisted exact rational computations, two separately implemented verification programs, manuscript writing and controlled publication. AI is not a scientific author.

## Abstract

We construct one exact Euclidean labelled squared-distance metric on 43 norm-one vertices by conservatively adjoining one perpendicular unit point to the frozen Polybrutus 37D scientific source. All 861 predecessor labelled squared distances remain unchanged; 42 new squared distances equal 2 at unit radius, giving 903 unordered labelled pairs. The exact rational centered Euclidean Gram is PSD with rank 38, proving minimum affine realization dimension 38 for this complete labelled metric. For a separately stipulated ideal unit-mass complete weighted graph with edges w_ij=1/d_ij², the Laplacian has eigenvalue 43/2 of multiplicity 32 and eigenvalue 21 of multiplicity 4. Two independently implemented exact-rational programs verify all pairs, parent conservation, Gram rank and selected graph eigenspaces. None of these mathematical statements implies physical extra spacetime dimensions or measured resonances.

## 1. Frozen source and conservative orthogonal extension

The immutable predecessor is Topbrutus/polybrutus-d37 at scientific tag v0.1.0 and source commit 294be7f814f2c8770af806833e98194fb192afff, DOI 10.5281/zenodo.24413346. Its original certificate records 42 unit points and 861 labelled squared distances, SHA-256 f649594103c43927126174dc5d807777bb976bf4feae2e9a1c0ee73aa67f3b2b. The complete frozen science tag and all older ancestors are preserved unmodified under prior/d37/.

For each predecessor unit vector p_i, i=1,...,42, embed q_i=(p_i,0) into one additional coordinate. Append one perpendicular norm-one vector q43=e39. Every old pairwise dot product is unchanged; thus every one of 861 labelled parent squared distances remains identical. By perpendicularity ||q43-q_i||²=||q43||²+||q_i||²-2<q43,q_i>=2. Hence exactly C(43,2)=903 unordered labelled pair squared distances arise, consisting of 861 inherited and 42 new pairs.

## 2. Exact rational Euclidean Gram and minimal affine rank

Seven canonical unit points form a rational rank-three Gram block with products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit points span an independent rank-four block with diagonal inner product one and distinct-pair products -1/4. Thirty-one further mutually perpendicular unit axes are orthogonal to the first two blocks. The 43-by-43 exact rational Gram G is positive semidefinite by construction, has unit diagonal and rank 3+4+31=38.

For complete squared-distance matrix D, set D_ii=0 and D_ij=2-2G_ij for distinct labels. With H=I-(1/43)11^T, the uniquely determined centered Gram is B=-(1/2)HDH=HGH; exact rank(B)=38. Two canonical vertices are antipodal, so their affine span contains the origin. Since B is uniquely determined by the full labelled squared-distance metric, every Euclidean realization requires affine dimension at least 38.

## 3. Exhaustive labelled-distance arithmetic

For unit radius R=1 let k=4d²/R² be the normalized squared-distance index. Exactly eight classes occur: k=4 six times; k=7 twice; k=8 872 times; k=9 three times; k=10 ten times; k=12 six times; k=13 twice; and k=16 twice. These counts total 903. The exact sum of all squared unordered pair distances is 7271/4.

Both checkers verify all 861 predecessor labelled values against the frozen certificate, rather than comparing histograms alone. A deliberate label permutation retains the unlabeled histogram but violates the parent correspondence; it is rejected by the controls.

## 4. Separately selected ideal complete graph

Independently choose a mathematical complete weighted graph with scalar unit masses and symmetric weights w_ij=1/d_ij². The new perpendicular vertex joins each of 42 predecessor vertices uniformly with edge weight 1/2. The weighted Laplacian therefore has exact block form:

L38 = [ L37+(1/2)I42        -(1/2)1_42 ]
      [ -(1/2)1_42^T               21 ].

Its exact uniform-join characteristic polynomial recurrence is chi43(z)=z(z-43/2)chi42(z-1/2)/(z-1/2). With f6(z)=chi_P7(z)/z the inherited nonzero canonical seven-point factor, the closed expression is chi43(z)=z(z-43/2)^32(z-21)^4*f6(z-18).

Two separate exact-rational matrix-rank calculations confirm nullity(L38-(43/2)I)=32 and nullity(L38-21I)=4. These are dimensionless theoretical eigenmodes of the stipulated graph, not measured hertz or physical observables.

## 5. Reproducibility, separate implementations and adversarial checks

The primary src/verify_d38.py constructs exact Gram entries, derives all 903 labelled squared distances, confirms all 861 predecessor values against the SHA-256 checked frozen parent, and checks centered Gram and chosen graph eigenspaces. Four deliberately incorrect cases are rejected: a corrupted distance, repeated vertex, non-orthogonal tilt and histogram-preserving label swap.

The independent tests/independent_d38.py program directly constructs each squared distance using Python fractions.Fraction structural rules without importing the main Gram generator. It checks the original parent SHA-256, all 861 conserved labels, all 903 present pair labels, exact centered Euclidean rank and selected graph multiplicities. Both scripts output the same SHA-256 for the deterministic results/certificate_d38.json; a public GitHub Actions workflow checks both and that certificate bytes are stable.

Reproduce by installing Python and sympy==1.14.0, then running python src/verify_d38.py and python tests/independent_d38.py. The frozen science tag and all predecessor records are included in prior/d37/.

## 6. Scientific scope and references

Orthogonal Euclidean extensions and uniform weighted graph joins are established mathematical techniques. The result is one exactly reproducible synthetic labelled-distance metric and a separately stipulated ideal graph, not a novel universal theorem. No additional physical spacetime dimension, observed particle, measured resonance or independently human peer-reviewed scientific claim is asserted.

References: Gabriel St-Pierre, Polybrutus 37D (2026), DOI 10.5281/zenodo.24413346; frozen scientific source github.com/Topbrutus/polybrutus-d37 v0.1.0; canonical seven-point geometry github.com/Topbrutus/polybrutus.
