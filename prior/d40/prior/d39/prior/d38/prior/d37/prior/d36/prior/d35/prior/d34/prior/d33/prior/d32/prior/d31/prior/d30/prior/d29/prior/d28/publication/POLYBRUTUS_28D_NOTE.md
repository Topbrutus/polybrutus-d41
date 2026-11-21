# Polybrutus 28D: Exact 33-point Euclidean geometry and a twenty-twofold ideal graph mode

**Gabriel St-Pierre** - Independent researcher  
**20 November 2026**  
**Exact mathematical case study / preprint**  
**AI assistance disclosure:** OpenAI ChatGPT (Astra) materially assisted mathematical verification, independently implemented scripts, manuscript writing, and controlled publication. AI is not an author.

## Abstract

We give an exactly reproducible Euclidean labelled squared-distance metric on 33 unit vertices by conservatively appending one perpendicular point to the frozen Polybrutus 27D configuration. All 496 predecessor labelled squared distances remain unchanged and 32 new squared distances equal 2 at unit radius; there are 528 unordered labelled pairs. An exact centered Euclidean Gram matrix is positive semidefinite and has rank 28, proving the minimal affine dimension of any Euclidean realization of this complete metric. In a separate stipulated unit-mass complete weighted graph with edge weights w(i,j)=1/d(i,j)^2, the Laplacian has eigenvalue 33/2 with exact multiplicity 22 and eigenvalue 16 with multiplicity 4. The metric and these ideal-model properties are independently verified by two exact-rational programs. These findings make no claim about physical spacetime dimensions or measured resonances.

## 1. Immutable predecessor and geometric construction

The scientific predecessor is Topbrutus/polybrutus-d27 source tag v0.1.0, frozen commit a311d63291d9254664340b434c5a13162d0658c1, DOI 10.5281/zenodo.24410593. Its original certificate records 32 unit points and 496 labelled squared distances. Its SHA-256 is e464dab0dc5d6586a94a714f1fee06b93710b3ba2763d710501970e6a3e399c4. The complete frozen predecessor tag is archived unmodified in prior/d27/.

Define q_i=(p_i,0) for each previous point i=1,...,32, and append an orthogonal unit basis point q33=e29. Every former pairwise squared distance remains the same because zero-coordinate embedding preserves scalar products. Every new pair satisfies ||q33-q_i||^2=1+1-2<q33,q_i>=2. Hence all C(33,2)=528 labelled squared distances consist of 496 frozen parent values plus 32 values equal to 2.

## 2. Rational Gram realization and exact affine dimension

The original seven canonical points form a rank-three Gram block with inner products 3*x_i*x_j+y_i*y_j+z_i*z_j. Five centered simplex unit points form a rank-four block: scalar product 1 on the diagonal, -1/4 between distinct vertices. Twenty-one further mutually perpendicular unit axis points form 21 independent coordinate directions with all cross-block products zero. The 33-by-33 Gram matrix G has rational entries, positive-semidefinite factorization and unit diagonal. Its exact rank is 3+4+21=28.

Let D have entries D_ii=0 and D_ij=2-2G_ij for i not equal j. The centered Euclidean distance Gram matrix with H=I-(1/33)11^T satisfies B=-(1/2)HDH=HGH. The canonical first and fifth points are antipodal; the origin is contained in the affine hull. Exact linear algebra proves rank(B)=28. A complete labelled Euclidean distance matrix uniquely determines B, so no Euclidean realization of this metric has affine dimension smaller than 28. This is a statement about a selected geometric configuration, not physical spacetime.

## 3. Complete labelled arithmetic invariants

For unit radius the normalized squared-distance index k=4d^2/R^2 has eight exact families: 4 occurs 6 times; 7 occurs 2 times; 8 occurs 497 times; 9 occurs 3 times; 10 occurs 10 times; 12 occurs 6 times; 13 occurs 2 times; and 16 occurs 2 times. The counts total exactly 528. The sum of all 528 squared pair distances is 4271/4.

Both programs compare each of the 496 original labelled pair values against the frozen D27 JSON certificate, rather than using only the histogram. This guards against a label permutation that preserves the distribution of distances while corrupting the exact labelled structure.

## 4. Stipulated ideal complete weighted graph

For a deliberately chosen mathematical model, assign scalar unit masses to the 33 vertices and complete edge weights w_ij=1/d_ij^2. The new perpendicular point has weight 1/2 to every predecessor. The Laplacian thus has uniform-join form:

L28 = [ L27 + (1/2)I32    -(1/2)1_32 ]
      [ -(1/2)1_32^T             16 ].

Writing chi_n(z) for the characteristic polynomial, the recurrence is chi33(z)=z(z-33/2)chi32(z-1/2)/(z-1/2). The inherited canonical seven-point factorization is chi33(z)=z(z-33/2)^22(z-16)^4*f6(z-13), where f6(z)=chi_P7(z)/z is the nonzero factor from the original seven-vertex chosen graph.

Two independently executed exact-rational Laplacian rank computations give nullity(L28-(33/2)I)=22 and nullity(L28-16I)=4. These are dimensionless spectral values from the stipulated graph, not observed sound, light, physical dimensions or particles.

## 5. Independent reproduction and adversarial tests

The primary script src/verify_d28.py creates the rational Gram matrix, derives all 528 labelled squared distances, compares 496 frozen predecessor pair values, and computes centered Gram and Laplacian ranks using exact symbolic arithmetic. It includes four controlled incorrect configurations: a damaged squared distance, a duplicated point, a nonperpendicular tilt and a relabelling that retains the unlabeled histogram.

The independent script tests/independent_d28.py reconstructs all labelled distances using fractions.Fraction rules without importing the primary Gram generator. It also checks the frozen predecessor certificate checksum, all old labels, the exact centered Gram rank and the graph eigenspace dimensions.

To reproduce, install Python and sympy==1.14.0. Run python src/verify_d28.py followed by python tests/independent_d28.py. Both must print identical SHA-256 values of results/certificate_d28.json. GitHub Actions performs the same checks on committed source and the frozen scientific tag.

## 6. Limits and references

Orthogonal Euclidean extensions and uniformly joined weighted graph Laplacians are known methods. This work is a specific falsifiable mathematical case study, not a novel universal theorem. It does not measure oscillations or demonstrate additional physical spacetime dimensions, and has not undergone independent human peer review.

References: Gabriel St-Pierre, Polybrutus 27D: Exact 32-Point Euclidean Metric and Twenty-Onefold Ideal Graph Modes (2026), https://doi.org/10.5281/zenodo.24410593 . Frozen GitHub source https://github.com/Topbrutus/polybrutus-d27/tree/v0.1.0 . Canonical Polybrutus https://github.com/Topbrutus/polybrutus/blob/main/SPECIFICATION.md .
