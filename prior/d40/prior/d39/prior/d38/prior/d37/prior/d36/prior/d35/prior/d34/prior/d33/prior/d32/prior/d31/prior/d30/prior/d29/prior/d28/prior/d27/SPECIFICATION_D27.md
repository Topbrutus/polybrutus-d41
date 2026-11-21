# Polybrutus D27 — Exact orthogonal Euclidean extension

## Frozen parent

The only source is `Topbrutus/polybrutus-d26@v0.1.0` with scientific source commit `5334c0fae387338314ea6e0fb8450bb3c1a5fcfa`, DOI `10.5281/zenodo.24410491`. The canonical parent certificate SHA-256 is `1b109ad0e5e7df765a00ca3fc5566d9d3381c96be1dd2e83bfe76c049159707a`, with 31 labelled unit points and 465 exact labelled squared distances.

Embed all 31 predecessor vertices `p_i` as `q_i=(p_i,0)` in one additional storage coordinate. Append a new perpendicular norm-one vertex `q32=e28`. The distances between predecessor vertices are unchanged; `||q32-q_i||²=||q32||²+||q_i||²-2<q32,q_i>=2`. Hence `C(32,2)=496=465+31`.

## Exact Gram and minimal affine dimension

Seven canonical vertices define the original rank-three block, with scalar products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex vertices satisfy norm squared one and pairwise product `-1/4` for different labels, giving a rank-four block. Twenty further unit basis points on mutually perpendicular axes are also perpendicular to the previous blocks. The rational Gram matrix `G` is PSD, with all diagonal entries one and `rank(G)=3+4+20=27`.

For pairwise squared distances, `D_ij=2-2G_ij` (distinct labels) and `D_ii=0`. Let `H=I-11^T/32`. The centered Euclidean Gram matrix obeys `B=-(1/2)HDH=HGH`; two canonical points are antipodal, placing the origin in the affine hull. Its exact rank is 27, proving minimum Euclidean affine realization dimension 27 for the complete labelled metric.

## Exact invariants

```text
Vertices:                          32
Unordered labelled squared pairs:  496
Frozen parent squared pairs:       465/465 preserved
New squared distances:             31 x 2R²
Normalized k = 4d²/R²:             4:6  7:2  8:465  9:3  10:10  12:6  13:2  16:2
Sum of squared distances (R=1):    4015/4
Centered Euclidean affine rank:    27
```

## Separately stipulated weighted graph model

Assign each vertex unit scalar mass and each distinct edge weight `w_ij=R²/d_ij²`. The perpendicular vertex joins all 31 old points with edge weight 1/2, yielding:

```text
L27 = [ L26 + (1/2)I31    -(1/2)1_31 ]
      [ -(1/2)1_31^T          31/2   ]
```

The uniform-join characteristic recurrence is `chi32(z)=z(z-16)chi31(z-1/2)/(z-1/2)`. By the preceding factorization, `chi32(z)=z(z-16)^21(z-31/2)^4*f6(z-25/2)`, with `f6(z)=chi_P7(z)/z` from the canonical seven-vertex model. Both exact-rational rank verifiers give nullity `21` for `L27-16I` and `4` for `L27-(31/2)I`.

These are chosen abstract graph eigenvalues and are not measurements of spacetime, particles or oscillation frequencies.
