# Polybrutus D26 — Exact conservative Euclidean extension

## Immutable parent

Use only `Topbrutus/polybrutus-d25@v0.1.0` and its scientific source commit `0dda9ed43622385d3e5fdddef9336458a32ed5ed`, Zenodo DOI `10.5281/zenodo.24410435`. The canonical parent certificate SHA-256 is `1cad349945e85a0700906762b21c3644ee9188308098b340e99acc525cac4ffd`, with 30 labelled unit vertices and 435 exact squared distances. All parent labels and values are frozen.

Embed the parent points as `q_i=(p_i,0)` for i=1..30. Append a perpendicular unit basis point `q31=e27`. Since each old point has norm one, every new squared distance is `||q31-q_i||²=1+1-0=2`. Hence 435 preserved old pairs and 30 new pairs give `C(31,2)=465`.

## Exact constructive Gram matrix

The seven canonical points form a rank-three block with products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit vertices satisfy diagonal product 1 and off-diagonal product `-1/4`, contributing rank four. Nineteen additional mutually perpendicular unit coordinate vectors have zero scalar products to either preceding block, contributing rank 19. The 31x31 Gram matrix `G` is rational, positive semidefinite, with diagonal ones and `rank(G)=3+4+19=26`.

For distinct points, the exact squared-distance matrix is `D_ij=2-2G_ij` and `D_ii=0`. With `H=I-11^T/31`, the centered Euclidean Gram matrix satisfies `B=-(1/2)HDH=HGH`. The first and fifth core vertices are antipodal, so the origin lies in the affine hull. Consequently `rank(B)=26`, which is the minimum possible affine dimension of any Euclidean realization of the same complete labelled metric.

## Arithmetic invariants

```text
vertices                  31
labelled pair distances   465
parent pairs unchanged    435/435
new squared distances     30 x 2R²
k=4*d²/R²                 4:6 7:2 8:434 9:3 10:10 12:6 13:2 16:2
sum(d²), R=1              3767/4
affine rank               26
centered Euclidean rank   26
```

## Stipulated graph Laplacian

Choose unit scalar masses and a **hypothetical complete weighted graph** with `w_ij=R²/d_ij²`. A new perpendicular unit vertex has weight 1/2 to each of the 30 predecessor vertices. Consequently:

```text
L26 = [ L25 + I30/2      -1_30/2 ]
      [ -1_30^T/2             15 ]
```

The characteristic-polynomial uniform-join recurrence is `chi31(z)=z(z-31/2)chi30(z-1/2)/(z-1/2)`. From the predecessor case, the formal factorization is `chi31(z)=z(z-31/2)^20(z-15)^4*f6(z-12)`, where `f6(z)=chi_P7(z)/z` is the nonzero factor of the canonical seven-vertex model. Both exact-rational rank checks show `nullity(L26-(31/2)I)=20` and `nullity(L26-15I)=4`.

This particular mathematical graph makes **no** claim about new physical spacetime dimensions, measured hertz, or particles.
