# Polybrutus 34D — Exact conservative orthogonal extension

## Frozen scientific predecessor

Use only `Topbrutus/polybrutus-d33@v0.1.0`, scientific source commit `bc78e6c6fe1755a069d8860a0144cbebbb1168a6`, DOI `10.5281/zenodo.24411569`. The immutable parent certificate SHA-256 is `86516b8a6393bf9fdf4e698d6e439d4db783a71956730ac35bbc89e2c8e7b734`. It records 38 unit points and 703 labelled exact squared distances.

Embed the prior 38 unit vertices as `q_i=(p_i,0)` and append one perpendicular norm-one unit point `q39=e35`. Every old labelled squared pair distance is preserved. For each old point, `||q39-q_i||²=1+1-2<q39,q_i>=2`. Therefore `C(39,2)=741=703+38` squared pair distances.

## Exact rational Euclidean Gram and minimum affine rank

Seven canonical unit points span an exact rank-three rational Gram block whose inner products are `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit points span a rank-four block, with diagonal inner products 1 and off-diagonal products `-1/4`. Twenty-seven additional pairwise perpendicular norm-one unit axes are orthogonal to both blocks. Hence the full rational Gram matrix `G` is positive semidefinite, has unit diagonal and exact rank `3+4+27=34`.

Define squared-distance matrix `D_ii=0` and `D_ij=2-2G_ij` for all distinct i,j. The centering projector is `H=I-11ᵀ/39`. The uniquely determined centered Euclidean Gram is `B=-(1/2)HDH=HGH`. The original first and fifth canonical unit points are antipodal and put the origin in the affine hull. Exact rational rank `rank(B)=34` proves that every Euclidean realization of this **complete labelled metric** has minimum affine dimension 34.

## Exact invariants at R=1

```text
Unit vertices                         39
Unordered labelled squared distances  741
Frozen D33 predecessor pair labels    703/703 preserved
New squared distances                 38 x 2
Histogram k=4*d²/R²                   4:6  7:2  8:710  9:3  10:10  12:6  13:2  16:2
Sum of squared pair distances         5975/4
Affine and centered Gram rank         34
```

## Separately specified ideal graph model

Assign scalar unit masses and symmetric complete edge weights `w_ij=R²/d_ij²`, a mathematical model choice. The new vertex joins all 38 predecessor vertices at uniform weight 1/2:

```text
L34 = [ L33 + I38/2        -1_38/2 ]
      [ -1_38ᵀ/2               19  ]
```

Its characteristic-polynomial recurrence is `chi39(z)=z(z-39/2)chi38(z-1/2)/(z-1/2)`. In the inherited seven-point factorization, `chi39(z)=z(z-39/2)^28(z-19)^4*f6(z-16)`, where `f6(z)=chi_P7(z)/z` denotes the nonzero factor of the original seven-point ideal graph. Independently checked exact-rational rank conditions give `nullity(L34-(39/2)I)=28` and `nullity(L34-19I)=4`.

These dimensionless eigenmodes are not physically measured frequencies.
