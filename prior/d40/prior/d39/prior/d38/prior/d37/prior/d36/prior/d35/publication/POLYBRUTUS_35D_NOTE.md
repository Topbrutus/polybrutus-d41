# Polybrutus 35D: Exact 40-point Euclidean distance geometry and a twenty-ninefold ideal graph mode

**Gabriel St-Pierre** — Independent researcher  
**21 November 2026**  
**Exact mathematical case study / scientific preprint**  
**AI assistance disclosed:** OpenAI ChatGPT (Astra) materially assisted exact calculations, two independently implemented verification scripts, writing and controlled publication. AI is not a scientific author.

## Abstract

We construct a particular complete labelled Euclidean squared-distance metric on 40 unit points by conservatively adding one perpendicular unit point to frozen Polybrutus 34D. All 741 predecessor labelled squared distances remain unchanged, and 39 new squared distances equal 2 at unit radius, yielding 780 unordered labelled pairs. The centered Euclidean Gram is positive semidefinite of exact rank 35, demonstrating minimum affine realization dimension 35. For a separately stipulated ideal equal-mass complete weighted graph with edge weights w_ij=1/d_ij², the Laplacian has eigenvalue 20 of exact multiplicity 29 and eigenvalue 39/2 of multiplicity 4. Two separately implemented exact-rational programs verify all distance labels, frozen parent conservation, Euclidean rank and chosen graph eigenspaces. These statements do not establish extra physical spacetime dimensions, measured frequencies or particle physics.

## 1. Frozen predecessor and exact orthogonal extension

The immutable predecessor is Topbrutus/polybrutus-d34 at scientific tag v0.1.0, source commit 21385a5875d472af8f035372096a46f5501cee1a, DOI 10.5281/zenodo.24412288. Its complete metric has 39 unit points and 741 labelled unordered squared distances; its certificate SHA-256 is 95274561481b735515a43dd185282cbe6b2424611f06efa5c6e5d100d433cfb6. The full parent tag source and all inherited ancestor files are preserved under prior/d34/ as ordinary files.

Embed each old unit point p_i as q_i=(p_i,0), for i=1,...,39. Append a perpendicular norm-one vector q40=e36. Every old dot product is unchanged, so all 741 old labelled squared distances are identical. By orthogonality each new squared pair distance is ||q40-q_i||²=||q40||²+||q_i||²-2<q40,q_i>=2. There are C(40,2)=780 unordered labelled pairs, split into 741 inherited labels and 39 new pairs of squared distance 2.

## 2. Exact positive-semidefinite Gram construction

The original seven canonical unit vertices form a rational rank-three Gram block with scalar products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered regular-simplex unit vertices form a rank-four Gram block, with diagonal products one and off-diagonal products -1/4. Twenty-eight additional pairwise perpendicular unit axes are orthogonal to the first two blocks. Hence the 40-by-40 rational Gram G is PSD with unit diagonal and exact rank 3+4+28=35.

For the complete labelled squared-distance matrix, define D_ii=0 and D_ij=2-2G_ij for each distinct pair. Let H=I-(1/40)11^T. The unique centered Euclidean Gram is B=-(1/2)HDH=HGH. Exact rational arithmetic yields rank(B)=35. Since two original canonical points are antipodal, the origin belongs to the affine hull. Because the complete labelled squared-distance matrix uniquely determines B, every Euclidean realization of this particular metric has affine dimension at least 35.

## 3. Exhaustive labelled-distance arithmetic

For unit radius R=1, the normalized squared-distance class k=4d²/R² has eight exact frequencies: k=4 occurs 6 times; k=7 occurs 2 times; k=8 occurs 749 times; k=9 occurs 3 times; k=10 occurs 10 times; k=12 occurs 6 times; k=13 occurs twice; k=16 occurs twice. These counts total 780. The exact sum of all unordered squared pair distances is 6287/4.

Conservation is checked for all 741 parent pairs by their labels, not only by an unlabeled histogram. A deliberate relabelling of vertices can retain the histogram while breaking predecessor correspondences; the verification controls reject that case.

## 4. Separately stipulated ideal graph

For one chosen mathematical model, stipulate scalar unit masses and a complete weighted graph with symmetric edge weights w_ij=1/d_ij². The new perpendicular unit point joins all 39 predecessor vertices with weight 1/2. Thus the Laplacian has the uniform-join block identity:

L35 = [ L34+(1/2)I39     -(1/2)1_39 ]
      [ -(1/2)1_39^T          39/2 ].

Its characteristic polynomial satisfies chi40(z)=z(z-20)chi39(z-1/2)/(z-1/2). In the inherited seven-point canonical factorization, this yields chi40(z)=z(z-20)^29(z-39/2)^4*f6(z-33/2), where f6(z)=chi_P7(z)/z is the original nonzero canonical model factor.

Two separately implemented exact-rational Laplacian rank checks establish nullity(L35-20I)=29 and nullity(L35-(39/2)I)=4. These are dimensionless modes of a selected mathematical graph, not observed sound frequencies, measured resonances, particle discoveries or physical spacetime dimensions.

## 5. Reproducibility, independent programs and negative controls

The primary program src/verify_d35.py derives all exact rational Gram entries and 780 labelled squared distances, verifies the 741 old labels against the SHA-256 checked frozen parent certificate, and checks centered Gram and graph ranks. Four intentionally false cases are rejected: a corrupted pair distance, duplicated point, non-orthogonal tilt and histogram-preserving relabelling.

A separately implemented verifier tests/independent_d35.py directly reconstructs every labelled distance using Python fractions.Fraction structural rules, without importing the primary Gram generator. It independently validates the parent certificate SHA-256, old and new pair labels, centered Euclidean rank, and chosen graph nullities. Both emit an identical SHA-256 for results/certificate_d35.json.

For reproduction, install Python and sympy==1.14.0, run python src/verify_d35.py and python tests/independent_d35.py. The CI workflow runs both verifiers on the published scientific source and protects the frozen predecessor certificate and deterministic output.

## 6. Scope and references

Orthogonal Euclidean extension and uniformly joined weighted-graph Laplacians are established mathematical methods. This work is an exact case study of one labelled synthetic Euclidean metric and one chosen ideal graph, not a new universal theorem. It does not demonstrate 35 physical spacetime dimensions, experimentally measured resonances or new particles. Independent human peer review has not taken place.

References: Gabriel St-Pierre, Polybrutus 34D (2026), DOI 10.5281/zenodo.24412288; frozen scientific source github.com/Topbrutus/polybrutus-d34, tag v0.1.0; canonical seven-point geometry github.com/Topbrutus/polybrutus.
