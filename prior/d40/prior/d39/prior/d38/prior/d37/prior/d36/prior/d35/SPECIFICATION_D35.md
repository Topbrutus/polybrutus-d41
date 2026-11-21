# Polybrutus 35D — Exact conservative Euclidean extension

## Scientific predecessor, frozen

`Topbrutus/polybrutus-d34@v0.1.0` at scientific commit `21385a5875d472af8f035372096a46f5501cee1a`, DOI `10.5281/zenodo.24412288`. The original D34 certificate SHA-256 is `95274561481b735515a43dd185282cbe6b2424611f06efa5c6e5d100d433cfb6` (39 unit points; 741 labelled unordered squared distances).

Embed all predecessor points `q_i=(p_i,0)`, and append one perpendicular unit vector `q40=e36`. Every old scalar product and all 741 predecessor labelled distances are unchanged. For each predecessor, `||q40-q_i||²=1+1-2<q40,q_i>=2`. Therefore `C(40,2)=780=741+39` labelled unordered pair distances.

## Euclidean Gram and minimum affine dimension

The canonical seven unit points span a rational rank-three Gram block with scalar products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit vertices span a rank-four Gram block, whose diagonal is one and distinct inner products are `-1/4`. Twenty-eight additional mutually perpendicular unit axes are orthogonal to both blocks. The complete exact rational Gram matrix `G` is PSD, has unit diagonal, and rank `3+4+28=35`.

For the full labelled squared-distance matrix, set `D_ii=0` and `D_ij=2-2G_ij` for distinct labels. With `H=I-11ᵀ/40`, the centered Euclidean Gram matrix is uniquely `B=-(1/2)HDH=HGH`. Exact rational rank is **35**. The first and fifth canonical points are antipodal and their affine hull contains the origin. Hence **35 is the minimum affine Euclidean realization dimension** for this specific complete labelled metric.

## Exhaustive arithmetic invariants at R=1

```text
Vertices                           40
Labelled unordered squared pairs   780
Old labelled pairs preserved       741/741
New squared-distance pairs         39 × 2
Histogram k=4*d²/R²                4:6  7:2  8:749  9:3  10:10  12:6  13:2  16:2
Sum of squared pair distances      6287/4
Affine / centered Euclidean rank   35
```

## Separately selected ideal weighted graph

Choose scalar unit masses and a complete weighted graph with `w_ij=R²/d_ij²`. The appended point joins all predecessor vertices with weight `1/2`:

```text
L35 = [ L34 + I39/2    -1_39/2 ]
      [ -1_39ᵀ/2          39/2 ]
```

The uniform-join characteristic recurrence is `chi40(z)=z(z-20)chi39(z-1/2)/(z-1/2)`. In the inherited seven-point factor notation, `chi40(z)=z(z-20)^29(z-39/2)^4*f6(z-33/2)`, with `f6(z)=chi_P7(z)/z`. Independently calculated exact nullities are `nullity(L35-20I)=29` and `nullity(L35-(39/2)I)=4`. These are dimensionless values of the stipulated graph, not experimentally measured frequencies.
