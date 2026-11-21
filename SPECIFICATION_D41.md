# Polybrutus 41D — Exact rational Euclidean labelled squared-distance metric

## Frozen scientific predecessor

The only scientific parent is `Topbrutus/polybrutus-d40@v0.1.0`, frozen commit `2a9fa8265b709d6d971afb6c63b7b024011e67a1`, DOI `10.5281/zenodo.24414238`, and original D40 certificate SHA-256 `af361ae471397b6d91beed87fe65c268a7139c02fbc0187400152fb8b202ffca`. It contains 45 labelled norm-one vertices and 990 unordered labelled squared distances.

Embed each predecessor unit point `p_i` as `q_i=(p_i,0)`. Add a perpendicular unit vertex `q46=e42`, using one new coordinate. All 990 frozen labelled squared pair distances remain unchanged. For each old vertex, `||q46-q_i||²=1+1-2<q46,q_i>=2`. Therefore `C(46,2)=1035=990+45` unordered labelled squared pairs are present.

## Exact Gram, positive semidefiniteness and affine rank

Seven canonical norm-one vertices have a rational rank-three Gram block with entries `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex norm-one vertices have a rank-four Gram block with unit diagonal and pairwise off-diagonal products `-1/4`. Thirty-four further mutually perpendicular norm-one axes are orthogonal to these blocks. Thus the rational Gram `G` is PSD, has unit diagonal, and exact rank `3+4+34=41`.

Define `D_ii=0`, `D_ij=2-2G_ij`, and centering projector `H=I-11ᵀ/46`. The unique Euclidean centered Gram is `B=-(1/2)HDH=HGH`, whose exact rational rank is **41**. The original antipodal canonical pair places the origin in the affine hull. Therefore the *specific complete labelled metric* has minimum Euclidean affine realization dimension **41**, not merely a nominal coordinate count.

## Exhaustive unit-radius invariants

```text
Labelled unit points                     46
Unordered labelled squared-distance pairs 1035
Frozen D40 pair labels preserved         990/990
New squared-distance pairs               45 × 2
Histogram k=4*d²/R²                      4:6 7:2 8:1004 9:3 10:10 12:6 13:2 16:2
Sum of all squared distances (R=1)       8327/4
Minimum affine/centered Gram rank        41
Ambient storage coordinates             42
```

## Separately stipulated ideal complete weighted graph

Choose scalar unit masses and symmetric complete graph weights `w_ij=R²/d_ij²`. The new perpendicular vertex joins every one of 45 older vertices at weight 1/2. The block Laplacian is

```text
L41 = [ L40 + I45/2       -1_45/2 ]
      [ -1_45ᵀ/2             45/2 ].
```

The uniform-join recurrence is `chi46(z)=z(z-23)chi45(z-1/2)/(z-1/2)`. With inherited original seven-point nonzero graph factor `f6(z)=chi_P7(z)/z`, the model expression is `chi46(z)=z(z-23)^35(z-45/2)^4*f6(z-39/2)`. Two separately implemented exact rational checks give `nullity(L41-23I)=35` and `nullity(L41-(45/2)I)=4`. All these graph eigenvalues are dimensionless and not measured frequencies.
