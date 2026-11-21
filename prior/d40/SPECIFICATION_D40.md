# Polybrutus 40D — Exact Euclidean labelled squared-distance specification

## Immutable scientific parent

The only parent is `Topbrutus/polybrutus-d39@v0.1.0`, frozen scientific commit `a43cd3c67017fb6162005fd521c156cb78795beb`, DOI `10.5281/zenodo.24413712`, original D39 certificate SHA-256 `8cd6731c294537d74c04a3501fe0da76a2cee6d496a2d2d5dfcfb0a26bb9304f`. It contains 44 labelled unit vertices and 946 unordered labelled squared distances.

Embed every predecessor norm-one point as `q_i=(p_i,0)`; append one perpendicular norm-one point `q45=e41`. The 946 old labelled squared distances are unchanged. Each one of 44 new squared distances satisfies `||q45-q_i||²=1+1-2<q45,q_i>=2`. Therefore `C(45,2)=990=946+44` exact labelled unordered pair distances.

## Constructive exact rational Euclidean Gram

Seven canonical unit vertices give a rational rank-three Gram block with products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit vertices form a separate rank-four Gram block with diagonal inner products 1 and off-diagonal products `-1/4`. Thirty-three additional perpendicular unit coordinate axes are mutually orthogonal and orthogonal to both blocks. The full 45-by-45 rational Gram `G` has unit diagonal, is positive semidefinite by construction, and has `rank(G)=3+4+33=40`.

For the complete labelled squared-distance matrix, `D_ii=0` and `D_ij=2-2G_ij` for distinct pairs. Let `H=I-11ᵀ/45`. The uniquely determined centered Euclidean Gram is `B=-(1/2)HDH=HGH`. Exact arithmetic yields `rank(B)=40`. Two original canonical vertices are antipodal, placing the origin in their affine hull. Thus every Euclidean realization of this particular *complete labelled squared-distance metric* requires at least 40 affine dimensions, which the explicit construction attains.

## Exhaustive invariants (R=1)

```text
Norm-one labelled vertices               45
Unordered labelled squared-distance pairs 990
Conserved frozen D39 predecessor pairs    946/946
New squared-distance pairs               44 × 2
Histogram k=4*d²/R²                      4:6 7:2 8:959 9:3 10:10 12:6 13:2 16:2
Sum of all squared distances             7967/4
Minimum affine and centered Gram rank    40
Coordinate storage count                 41
```

## Separately selected ideal complete weighted graph

Independently stipulate scalar unit masses and symmetric complete graph edge weights `w_ij=R²/d_ij²`. The appended vertex joins each old vertex at uniform weight 1/2. Its weighted Laplacian is

```text
L40 = [ L39 + I44/2       -1_44/2 ]
      [ -1_44ᵀ/2               22 ].
```

The exact uniform-join characteristic recurrence is `chi45(z)=z(z-45/2)chi44(z-1/2)/(z-1/2)`. If `f6(z)=chi_P7(z)/z` denotes the inherited canonical seven-point graph's nonzero factor, the closed expression is `chi45(z)=z(z-45/2)^34(z-22)^4*f6(z-19)`. Two separately implemented exact rank computations establish `nullity(L40-(45/2)I)=34` and `nullity(L40-22I)=4`.

These are dimensionless chosen-model eigenvalues, not measured frequencies.
