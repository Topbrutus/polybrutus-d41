# Polybrutus 33D — Exact conservative orthogonal extension

## Immutable predecessor

Parent: `Topbrutus/polybrutus-d32@v0.1.0`, frozen scientific commit `cb63ed7ab447b13e251d92571fe0f1062db31cca`, DOI `10.5281/zenodo.24411375`, parent certificate SHA-256 `66b7605ab168d99ee33014ed1d61b5fbb3915c32041bdf44dd529c0bb56bb016`. This certifies 37 norm-one unit points and **666** labelled squared-distance pairs.

Embed the predecessor points `q_i=(p_i,0)` and append a new perpendicular norm-one vector `q38=e34`. Each predecessor squared pair distance is unchanged. Every new pair has squared distance `||q38-q_i||²=1+1-2<q38,q_i>=2`, so `C(38,2)=703=666+37` exact labelled pair distances.

## Rational Euclidean Gram and minimum affine rank

Seven original canonical unit vertices have an exact rational rank-three Gram block with products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit points have an independent rank-four Gram block (diagonal products 1 and distinct-pair products -1/4). **Twenty-six** further mutually perpendicular unit-axis vertices add 26 independent directions orthogonal to both blocks. Thus `G` is rational PSD, has unit diagonal, and `rank(G)=3+4+26=33`.

For distinct i,j define `D_ij=2-2G_ij`, with `D_ii=0`. Let `H=I-11ᵀ/38`. The unique Euclidean centered Gram matrix for the complete labelled distance metric is `B=-(1/2)HDH=HGH`. The original first and fifth canonical points are antipodal and their affine hull includes the origin. Direct exact arithmetic gives `rank(B)=33`, so every Euclidean realization of this complete labelled squared-distance metric has affine dimension at least 33.

## Checked invariants, R=1

```text
Vertices                          38
Unordered labelled pairs          703
Conserved parent pairs            666/666
New squared-distance pairs        37 x 2
k = 4*d²/R² histogram             4:6, 7:2, 8:672, 9:3, 10:10, 12:6, 13:2, 16:2
Sum of squared pair distances     5671/4
Affine / centered Gram rank       33
```

## Stipulated mathematical weighted graph

Choose a complete graph with unit scalar masses and `w_ij=R²/d_ij²`, independently of the geometric construction. The new vertex adds uniform coupling 1/2:

```text
L33 = [ L32 + I37/2      -1_37/2 ]
      [ -1_37ᵀ/2           37/2   ]
```

The characteristic-polynomial recurrence is `chi38(z)=z(z-19)chi37(z-1/2)/(z-1/2)`. In the inherited seven-point factor notation, `chi38(z)=z(z-19)^27(z-37/2)^4*f6(z-31/2)`, where `f6(z)=chi_P7(z)/z`. Two exact-rational rank programs check `nullity(L33-19I)=27` and `nullity(L33-(37/2)I)=4`. These are dimensionless theoretical eigenmodes, not observed frequencies.
