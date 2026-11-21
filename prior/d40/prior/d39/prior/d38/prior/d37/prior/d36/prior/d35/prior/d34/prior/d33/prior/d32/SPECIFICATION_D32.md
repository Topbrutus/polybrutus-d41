# Polybrutus 32D — Exact orthogonal Euclidean completion

## Frozen source and construction

The immutable parent is `Topbrutus/polybrutus-d31@v0.1.0`, scientific commit `0c1f546a6c5442bc83a5bdde506132cce20b6e05`, DOI `10.5281/zenodo.24411160`; the exact 36-point parent certificate SHA-256 is `7c373c50d61ca2964eb832fc6fb8178590715e603f5fdada122575f624f60e8c`.

Embed the 36 old unit points as `q_i=(p_i,0)`, then append one new perpendicular norm-one point `q37=e33`. All 630 parent labelled squared distances are conserved; every new pair satisfies `||q37-q_i||²=1+1-2<q37,q_i>=2`. Thus `C(37,2)=666=630+36`.

## Constructive Gram proof of minimum affine dimension

Seven original canonical unit points define a rank-three rational Gram block with inner products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered simplex unit points define an independent rank-four block with diagonal one and pair products `-1/4`. Twenty-five further mutually perpendicular unit-axis points are orthogonal to these blocks. The exact rational Gram matrix `G` is positive semidefinite with `G_ii=1` and `rank(G)=3+4+25=32`.

Let the complete squared-distance matrix satisfy `D_ii=0` and `D_ij=2-2G_ij` for distinct points. With `H=I-11ᵀ/37` the unique centered Euclidean Gram matrix is `B=-(1/2)HDH=HGH`. The original antipodal core pair puts the origin in the affine hull. An exact rational calculation gives `rank(B)=32`, proving minimal affine Euclidean dimension 32 for this complete labelled metric.

## Exact pair invariants, R=1

```text
Vertices                           37
Unordered labelled pairs           666
Frozen parent labelled pairs       630/630 conserved
New squared-distance pairs         36 x 2
k = 4*d²/R² histogram              4:6, 7:2, 8:635, 9:3, 10:10, 12:6, 13:2, 16:2
Sum of all squared pair distances  5375/4
Centered Gram / affine rank        32
```

## Separately selected ideal complete weighted graph

Stipulate scalar masses one and `w_ij=1/d_ij²` for each complete-graph edge. New point 37 joins the old 36 vertices uniformly with weight 1/2:

```text
L32 = [ L31 + I36/2     -1_36/2 ]
      [ -1_36ᵀ/2             18 ]
```

The uniform-join characteristic recurrence is `chi37(z)=z(z-37/2)chi36(z-1/2)/(z-1/2)`. In the inherited canonical seven-point notation, `chi37(z)=z(z-37/2)^26(z-18)^4*f6(z-15)`, where `f6(z)=chi_P7(z)/z`. Two separately implemented exact-rational rank calculations verify `nullity(L32-(37/2)I)=26` and `nullity(L32-18I)=4`.

These dimensionless graph modes are not experimentally measured frequencies.
