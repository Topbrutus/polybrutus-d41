# Polybrutus 30D — Conservative exact Euclidean extension

## Source frozen before this construction

Parent: `Topbrutus/polybrutus-d29@v0.1.0`, scientific commit `cbe179697cc1a287fc0637912ea6a6fd890eeca4`, DOI `10.5281/zenodo.24410725`, exact parent certificate SHA-256 `d718f784bf20dccc1b8c08821466b352316054900fba20705bc28d79763a4731`. Its 34 labelled unit points provide **561 pairwise squared distances** that must be preserved by label, not just histogram.

Embed all old vectors as `q_i=(p_i,0)` and append new perpendicular norm-one vector `q35=e31`. Old scalar products remain fixed, and each new squared pair distance satisfies `||q35-q_i||²=1+1-2<q35,q_i>=2`. Therefore the 35-point metric contains **595 = 561+34** labelled pairs.

## Explicit rational Gram metric and affine dimension

Seven canonical unit vertices form a rank-three Gram block with products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit vertices form a rank-four Gram block, with diagonal 1 and off-diagonal -1/4. Twenty-three additional orthogonal unit-axis points are mutually perpendicular and orthogonal to both blocks. The resulting rational Gram matrix `G` is PSD, has unit diagonal, and `rank(G)=3+4+23=30`.

Let `D_ii=0`, `D_ij=2-2G_ij` for i ≠ j, and `H=I-11^T/35`. Then `B=-(1/2)HDH=HGH` is the exact centered Euclidean Gram matrix. The first and fifth canonical points are antipodal, putting the origin in their affine hull. Exact rational rank computations yield `rank(B)=30`, the minimum affine dimension for any Euclidean realization of this **complete labelled** squared-distance metric.

## Arithmetic invariants at R=1

```text
Vertices                         35
Labelled unordered pairs         595
Conserved D29 labelled pairs     561/561
New pairs of squared distance 2  34
Class k=4*d²/R²                  4:6 7:2 8:564 9:3 10:10 12:6 13:2 16:2
Sum of squared pair distances    4807/4
Affine dimension                 30
Centered Gram rank               30
```

## One selected ideal weighted complete graph

Stipulate scalar unit masses and `w_ij=R²/d_ij²` for each pair. The new point uniformly joins all 34 predecessors at weight 1/2:

```text
L30 = [ L29 + I34/2         -1_34/2 ]
      [ -1_34^T/2               17  ]
```

The uniform-join recurrence for characteristic polynomials is `chi35(z)=z(z-35/2)chi34(z-1/2)/(z-1/2)`. In the inherited canonical factor notation, `chi35(z)=z(z-35/2)^24(z-17)^4*f6(z-14)`, with `f6(z)=chi_P7(z)/z`. The two independent exact-rational verifiers establish `nullity(L30-(35/2)I)=24` and `nullity(L30-17I)=4`. This mathematical spectrum does not measure any physical frequency.
