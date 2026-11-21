# Polybrutus D24 — Exact conservative orthogonal extension

## Immutable predecessor

Frozen `Topbrutus/polybrutus-d23@v0.1.0`, scientific commit `23f71c9b9932894c502fb67c1326271daa0692e7`, DOI `10.5281/zenodo.24410035`. The original D23 certificate contains 28 norm-one points and 378 labelled squared distances. Its SHA-256 is `a3ef363bb5407ed54b07b3e1002c25d4dbdcbef946bb186862c891076c626724`.

Embed the 28 old unit points as `q_i=(p_i,0)` for `1<=i<=28`, then append one perpendicular unit basis vector `q29=e25`. For each old point, `||q29-q_i||² = ||q29||²+||q_i||² - 2<q29,q_i> =2`. Therefore all old labelled pairs are preserved and `C(29,2)=406=378+28`.

## Exact Gram construction

The original seven canonical vectors are encoded by their scalar products `3*x_i*x_j + y_i*y_j + z_i*z_j` in the three-dimensional core block. Five centered regular-simplex vectors have inner products `1` on the diagonal and `-1/4` between distinct points, contributing rank 4. Seventeen further perpendicular unit axes yield rank 17. Blocks are mutually orthogonal. Thus the complete 29-by-29 Gram matrix `G` has exact rank `3+4+17=24`, is positive semidefinite, and has diagonal ones.

For distinct points define squared distances `D_ij=2-2G_ij`, and let diagonal entries vanish. The centered Euclidean squared-distance Gram matrix is `B=-(1/2)HDH=HGH`, `H=I-11^T/29`. The first and fifth canonical points are antipodal; the origin is in the affine hull. Hence the exact affine rank and minimum dimension of any Euclidean realization of the full labelled distance metric are both **24**.

## Verified invariants, radius R=1

```text
class k=4*d²/R²: 4:6, 7:2, 8:375, 9:3, 10:10, 12:6, 13:2, 16:2
all 406 labelled pair distances: verified
all 378 prior labelled squared distances: unchanged
new squared distances: 28 x 2
sum of squared distances: 3295/4
```

## Ideal graph model

Defining edge weights `w_ij=1/D_ij` gives a complete weighted Laplacian. The new vertex uniformly joins 28 old vertices with weight 1/2. If `L23` denotes the predecessor's model:

```
L24 = [L23 + I28/2     -1_28/2]
      [-1_28^T/2         14   ]
```

The exact uniform-join characteristic recurrence is `chi29(z)=z(z-29/2)chi28(z-1/2)/(z-1/2)`. Substituting the predecessor factor gives `chi29(z)=z(z-29/2)^18(z-14)^4*f6(z-11)`, where `f6(z)=chi_P7(z)/z` is the nonzero factor of the canonical 7-vertex model. Independent rational rank computations establish nullities 18 and 4 for `L24-(29/2)I` and `L24-14I`, respectively.

No physical spacetime dimension or measured resonance is implied by this definition.
