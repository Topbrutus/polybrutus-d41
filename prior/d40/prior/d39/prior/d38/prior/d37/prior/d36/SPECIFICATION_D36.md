# Polybrutus 36D — Exact rational Euclidean metric specification

## 1. Immutable parent and conservative extension

Take the frozen `Topbrutus/polybrutus-d35@v0.1.0` scientific source, commit `72989a1fd2a7fb6e28522a0f2dcbe168924ca78d`, DOI `10.5281/zenodo.24412435`; its original 40-point certificate SHA-256 is `7390c10f46132faeca0ccd8b52050d5ed899fbd4a078b9dd02407d922a4e4f22`.

Embed all 40 parent unit points `q_i=(p_i,0)`, and append one perpendicular norm-one vector `q41=e37`. All **780** old labelled squared distances are exactly conserved. For every old point `||q41-q_i||²=1+1-2<q41,q_i>=2`. Consequently `C(41,2)=820=780+40` unordered labelled squared pair distances.

## 2. Gram construction and minimum affine dimension

Seven canonical unit points form the rational rank-three Gram block with inner products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit points form a separate rank-four block, with diagonal entries 1 and off-diagonal products `-1/4`. **Twenty-nine** further pairwise orthogonal unit axes add independent directions. Hence the exact rational Gram `G` is positive semidefinite, has unit diagonal, and rank `3+4+29=36`.

Set `D_ii=0` and `D_ij=2-2G_ij` for distinct labels. Let `H=I-11ᵀ/41`. The centered Euclidean Gram matrix uniquely derived from complete squared distances is `B=-(1/2)HDH=HGH`; exact rational rank is **36**. Two original canonical vertices are antipodal, placing the origin in the affine hull. The complete labelled metric therefore has minimum Euclidean affine realization dimension **36**.

## 3. Exhaustive invariants (R=1)

```text
Unit vertices                         41
Unordered labelled squared pairs      820
Conserved parent pair labels          780/780
New squared-distance pairs            40 x 2
Histogram k=4*d²/R²                   4:6, 7:2, 8:789, 9:3, 10:10, 12:6, 13:2, 16:2
Sum of all squared pair distances     6607/4
Affine and centered Gram rank         36
```

## 4. Independently stipulated ideal complete weighted graph

Select an **ideal mathematical model** of scalar unit masses with symmetric complete edge weights `w_ij=R²/d_ij²`. The new point connects to each of the previous 40 at weight 1/2. Its weighted Laplacian is

```text
L36 = [ L35 + I40/2        -1_40/2 ]
      [ -1_40ᵀ/2                20  ].
```

The exact uniform-join characteristic-polynomial recurrence is `chi41(z)=z(z-41/2)chi40(z-1/2)/(z-1/2)`. With `f6(z)=chi_P7(z)/z` the inherited canonical seven-point nonzero factor, `chi41(z)=z(z-41/2)^30(z-20)^4*f6(z-17)`. Exact ranks give `nullity(L36-(41/2)I)=30` and `nullity(L36-20I)=4`.

These eigenvalues are dimensionless mathematical quantities, **not measured Hz**.
