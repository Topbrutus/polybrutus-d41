# Polybrutus 28D - Exact conservative orthogonal extension

## Predecessor and invariants

Freeze `Topbrutus/polybrutus-d27@v0.1.0`, scientific commit `a311d63291d9254664340b434c5a13162d0658c1`, parent DOI `10.5281/zenodo.24410593`. Its certificate has SHA-256 `e464dab0dc5d6586a94a714f1fee06b93710b3ba2763d710501970e6a3e399c4`, with **32** unit vertices and **496** labelled squared distances.

Embed the existing 32 points as `q_i=(p_i,0)` and append the perpendicular norm-one point `q33=e29`. All old pair distances are unchanged. Each new pair has `||q33-q_i||²=||q33||²+||q_i||²-2<q33,q_i>=2`. Hence `C(33,2)=528=496+32` labelled squared pair distances.

## Rational Gram construction

The original seven canonical points form a rank-three block with inner products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit points form a rank-four block (diagonal products 1, distinct-pair products -1/4). **Twenty-one** mutually perpendicular unit-axis points are orthogonal to these blocks and to each other. The rational Gram matrix `G` is PSD with unit diagonal and exact rank `3+4+21=28`.

For distinct points `D_ij=2-2G_ij` and `D_ii=0`. Let `H=I-11^T/33`. The centered Euclidean Gram matrix satisfies `B=-(1/2)HDH=HGH`. The first and fifth canonical vertices are antipodal; the origin belongs to the affine hull. Therefore `rank(B)=28`, the exact minimum affine dimension for any Euclidean realization of the full labelled metric.

## Exact arithmetic, R=1

```text
vertices:                           33
labelled squared-distance pairs:    528
old labelled pairs preserved:       496/496
new squared-distance pairs:         32 x 2
k=4*d²/R² histogram:                4:6,7:2,8:497,9:3,10:10,12:6,13:2,16:2
sum of squared distances:           4271/4
affine/centered-Gram rank:          28
```

## One chosen ideal weighted graph

Assign unit scalar masses and complete edge weights `w_ij=R²/d_ij²`. The appended point joins all 32 old points uniformly at weight 1/2:

```text
L28 = [ L27 + I32/2       -1_32/2 ]
      [ -1_32^T/2             16  ]
```

For characteristic polynomials the uniform-join recurrence is `chi33(z)=z(z-33/2)chi32(z-1/2)/(z-1/2)`; with canonical 7-point factor `f6(z)=chi_P7(z)/z`, the inherited form is `chi33(z)=z(z-33/2)^22(z-16)^4*f6(z-13)`. Two exact-rational rank programs establish nullities `22` for `L28-(33/2)I` and `4` for `L28-16I`. No observed physical quantity is inferred.
