# Polybrutus 37D — Exact rational labelled Euclidean distance metric

## Frozen scientific parent

The only predecessor is `Topbrutus/polybrutus-d36@v0.1.0`, scientific commit `b1f8be2cdcdcae0d9911d90381a7b1304c84abf0`, DOI `10.5281/zenodo.24413083`. The original D36 certificate SHA-256 is `b189f74f7b38b378fa8f34011660817019c02601e5889378103fd9424207b331` (41 norm-one points and 820 labelled pair distances).

Embed each parent point as `q_i=(p_i,0)` and append one perpendicular unit vector `q42=e38`. This leaves all **820/820** parent pair labels unchanged. Each new pair has squared distance `||q42-q_i||²=1+1-2<q42,q_i>=2`, so `C(42,2)=861=820+41`.

## Constructive Euclidean Gram proof

Seven original canonical unit points form a rank-three rational Gram block, with scalar products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit points form a separate rank-four Gram block, with unit diagonal and distinct inner products `-1/4`. Thirty further unit coordinate axes are pairwise orthogonal and orthogonal to the first two blocks. Hence exact rational Gram matrix `G` is PSD with unit diagonal and `rank(G)=3+4+30=37`.

Let `D_ii=0` and `D_ij=2-2G_ij` for every distinct labelled pair. With `H=I-11ᵀ/42`, the complete labelled squared-distance matrix has the uniquely determined centered Euclidean Gram `B=-(1/2)HDH=HGH`. Exact rational calculation gives `rank(B)=37`. The first and fifth canonical points are antipodal, so the origin is in their affine hull. Consequently every Euclidean realization of this particular complete labelled metric has affine dimension at least 37, and the explicit realization achieves that dimension.

## Exact arithmetic invariants (R=1)

```text
Norm-one labelled vertices            42
Unordered labelled squared pairs      861
Conserved predecessor pair labels     820/820
New squared-distance pairs            41 x 2
Histogram k=4*d²/R²                   4:6, 7:2, 8:830, 9:3, 10:10, 12:6, 13:2, 16:2
Sum of all squared pair distances     6935/4
Affine and centered Gram rank         37
Coordinate storage count              38
```

## Separately chosen ideal weighted graph

Stipulate a complete mathematical graph with unit scalar masses and symmetric weights `w_ij=R²/d_ij²`. The newly appended vertex joins the 41 older vertices uniformly at weight 1/2:

```text
L37 = [ L36 + I41/2      -1_41/2 ]
      [ -1_41ᵀ/2             41/2 ].
```

The exact uniform-join characteristic recurrence is `chi42(z)=z(z-21)chi41(z-1/2)/(z-1/2)`. For the inherited canonical seven-point nonzero factor `f6(z)=chi_P7(z)/z`, the closed expression is `chi42(z)=z(z-21)^31(z-41/2)^4*f6(z-35/2)`. Two independent exact rational matrix-rank calculations give `nullity(L37-21I)=31` and `nullity(L37-(41/2)I)=4`. These values are dimensionless model quantities, not physically measured frequencies.
