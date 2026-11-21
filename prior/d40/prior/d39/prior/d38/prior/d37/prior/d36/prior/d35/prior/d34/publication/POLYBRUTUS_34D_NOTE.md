# Polybrutus 34D: Exact 39-point Euclidean distance geometry and a twenty-eightfold ideal graph mode

**Gabriel St-Pierre** - Independent researcher  
**21 November 2026**  
**Exact mathematical case study / scientific preprint**  
**AI disclosure:** OpenAI ChatGPT (Astra) materially assisted symbolic calculation, two independently implemented verifiers, manuscript writing and controlled release. AI is not a scientific author.

## Abstract

This note presents a specific exact Euclidean labelled squared-distance metric on 39 unit points by conservatively extending the frozen Polybrutus 33D scientific source with one perpendicular unit vertex. All 703 predecessor labelled squared pair distances are identical, while 38 new squared distances equal 2 at unit radius, giving 741 unordered labelled pairs. The exact centered Euclidean Gram matrix is positive semidefinite of rank 34, establishing the minimum affine realization dimension of this complete metric. In a separately stipulated equal-mass ideal complete weighted graph with w_ij=1/d_ij², the Laplacian has eigenvalue 39/2 of exact multiplicity 28 and eigenvalue 19 of multiplicity 4. Two independently implemented exact-rational programs reproduce all pair labels, predecessor conservation, Euclidean ranks and chosen graph modes. These properties do not assert additional physical spacetime dimensions or experimental resonances.

## 1. Frozen predecessor and conservative unit extension

The immutable scientific predecessor is Topbrutus/polybrutus-d33, scientific tag v0.1.0, exact source commit bc78e6c6fe1755a069d8860a0144cbebbb1168a6, Zenodo DOI 10.5281/zenodo.24411569. Its 38 unit vertices determine 703 exact labelled squared distances; the original certificate SHA-256 is 86516b8a6393bf9fdf4e698d6e439d4db783a71956730ac35bbc89e2c8e7b734. All frozen scientific parent sources and ancestor records are preserved unmodified in prior/d33/.

For each old unit point p_i, i=1,...,38, define q_i=(p_i,0) in one additional coordinate dimension. Append a perpendicular norm-one point q39=e35. Every former pair squared distance remains identical. For every old point, ||q39-q_i||²=||q39||²+||q_i||²-2<q39,q_i>=2. Hence the 741 unordered labelled pairs consist of exactly 703 unchanged predecessor values and 38 new squared distances of 2.

## 2. Constructive rational Gram and minimum affine dimension

The seven original canonical unit vertices form a rational rank-three Gram block with inner products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vertices span an independent rank-four Gram block, with diagonal products one and off-diagonal products -1/4. Twenty-seven additional pairwise perpendicular norm-one axes are orthogonal to these former blocks. Thus the 39-by-39 rational Gram matrix G is PSD, has unit diagonal and rank 3+4+27=34.

Define the complete labelled squared-distance matrix D by D_ii=0, D_ij=2-2G_ij for distinct labels. Let H=I-(1/39)11^T. The unique centered Euclidean Gram matrix is B=-(1/2)HDH=HGH. Exact-rational computation gives rank(B)=34. Two original canonical unit vertices are antipodal, so the origin belongs to their affine hull. Since B is determined uniquely by the complete labelled distances, every Euclidean realization of this particular metric needs at least 34 affine dimensions.

## 3. Exact exhaustive labelled distance arithmetic

With radius R=1 define the normalized squared-distance index k=4d²/R². Exactly eight classes occur: k=4 six times, k=7 twice, k=8 710 times, k=9 three times, k=10 ten times, k=12 six times, k=13 twice and k=16 twice. The class counts sum to 741. The exact sum of all unordered squared pair distances equals 5975/4.

Both programs check each of the 703 predecessor labelled pair values against the SHA-256 verified frozen parent certificate; a histogram alone would not protect the correspondence of labels. A deliberate swap of point labels retains the unlabeled histogram but breaks original pair correspondences, and is rejected by the negative controls.

## 4. Stipulated ideal complete weighted graph and eigenspaces

Separately choose a dimensionless ideal model: a complete graph with scalar mass one at each vertex and symmetric weights w_ij=1/d_ij². The new perpendicular unit point joins all 38 old vertices uniformly with edge weight 1/2. Its weighted graph Laplacian therefore has the exact block form:

L34 = [ L33+(1/2)I38        -(1/2)1_38 ]
      [ -(1/2)1_38^T             19  ].

The characteristic polynomial obeys chi39(z)=z(z-39/2)chi38(z-1/2)/(z-1/2). If f6(z)=chi_P7(z)/z denotes the nonzero canonical seven-point polynomial factor, the inherited decomposition reads chi39(z)=z(z-39/2)^28(z-19)^4*f6(z-16).

Two independently implemented exact-rational matrix rank computations give nullity(L34-(39/2)I)=28 and nullity(L34-19I)=4. These are dimensionless properties of a freely selected graph; not experimental sound frequencies, observations of particles, or evidence for physical 34-dimensional spacetime.

## 5. Independent reproducibility and adversarial controls

The primary src/verify_d34.py program constructs exact rational Gram entries, derives all 741 pair labels, verifies all 703 old pair values against the SHA-256 checked frozen parent, and computes exact centered Gram and ideal graph eigenspace ranks. It rejects four deliberately incorrect configurations: one damaged distance, a duplicated unit point, a nonperpendicular tilt and a histogram-preserving point-label swap.

The independent tests/independent_d34.py program reconstructs every distance using Python fractions.Fraction structural rules without importing the primary Gram generator. It independently checks predecessor SHA-256, all old and new pair labels, centered Euclidean rank and graph multiplicities. Both scripts output the identical deterministic SHA-256 for results/certificate_d34.json. The public GitHub Actions workflow reruns them and rejects nondeterministic certificate output.

Install Python and sympy==1.14.0, then run python src/verify_d34.py and python tests/independent_d34.py. The complete frozen D33 science tag resides under prior/d33/ without nested Git metadata.

## 6. Scientific limitations and references

Orthogonal Euclidean embedding and uniform graph Laplacian joining are established mathematical techniques. This is one exact reproducible synthetic metric plus one stipulated graph, not a new universal theorem. No measured resonances, additional physical spacetime dimensions, particle discoveries or independent human peer review are asserted.

References: Gabriel St-Pierre, Polybrutus 33D (2026), DOI 10.5281/zenodo.24411569; frozen source github.com/Topbrutus/polybrutus-d33, v0.1.0; original canonical seven-point geometry github.com/Topbrutus/polybrutus.
