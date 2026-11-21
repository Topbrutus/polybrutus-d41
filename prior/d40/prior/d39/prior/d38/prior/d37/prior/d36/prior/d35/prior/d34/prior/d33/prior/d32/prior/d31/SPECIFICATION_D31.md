# Polybrutus 31D — Exact Euclidean extension specification

## Frozen parent and conservative embedding

Use only `Topbrutus/polybrutus-d30@v0.1.0` (scientific commit `3c7ce4019a49ff4c54adb9271367aff5c62c1326`), DOI `10.5281/zenodo.24411018`. The frozen parent has 35 unit points, 595 labelled squared pair distances, and certificate SHA-256 `96c2f0693d849b2bab83c7ef4b25f647ecd3ee3a3ae460d3a32c5e9a73f1fb1c`.

Embed every old unit vector as `q_i=(p_i,0)`, and append `q36=e32`, a new perpendicular norm-one unit vector. All 595 prior labelled squared distances are preserved exactly; the 35 additional distances satisfy `||q36-q_i||²=||q36||²+||q_i||²-2<q36,q_i>=2`. Thus `C(36,2)=630=595+35`.

## Constructive Gram and minimum dimension

An exact rational PSD Gram matrix G has three independent canonical directions from seven core unit vertices (products `3*x_i*x_j+y_i*y_j+z_i*z_j`), four independent directions from five centered simplex unit vertices (distinct products `-1/4`), and 24 mutually perpendicular unit-axis points orthogonal to both blocks. Thus `rank(G)=3+4+24=31` and `G_ii=1`.

Define `D_ii=0` and `D_ij=2-2G_ij` for distinct pairs, and the centering matrix `H=I-11ᵀ/36`. Then `B=-(1/2)HDH=HGH`. Since the canonical first and fifth points are antipodal, the origin belongs to their affine hull. The centered Gram has exact rank 31, establishing the minimum Euclidean affine realization dimension for this complete labelled metric.

## Invariants (R=1)

```text
Vertices                         36
Labelled squared-distance pairs  630
Preserved predecessor pairs      595/595
New squared distances            35 × 2
Normalized histogram k=4d²/R²    4:6 7:2 8:599 9:3 10:10 12:6 13:2 16:2
Sum of all squared distances     5087/4
Affine / centered Gram rank      31
```

## Separately stipulated graph

Choose scalar unit masses and a mathematical complete weighted graph with `w_ij=R²/d_ij²`. The new point adds weight 1/2 to every old vertex:

```text
L31 = [ L30 + I35/2       -1_35/2 ]
      [ -1_35ᵀ/2              35/2 ]
```

The exact uniform-join characteristic recurrence is `chi36(z)=z(z-18)chi35(z-1/2)/(z-1/2)`. In the canonical factor notation, `chi36(z)=z(z-18)^25(z-35/2)^4*f6(z-29/2)`, with `f6(z)=chi_P7(z)/z`. Two separately implemented rank computations yield `nullity(L31-18I)=25` and `nullity(L31-(35/2)I)=4`.

The eigenvalues are dimensionless quantities of the selected graph, not measured frequencies.
