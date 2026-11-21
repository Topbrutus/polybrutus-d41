# Polybrutus 39D — Exact Euclidean labelled squared-distance metric

## Frozen scientific parent

Use only `Topbrutus/polybrutus-d38@v0.1.0`, scientific commit `ac8fd112fcd314092f106b210162ba1cc9879e9b`, DOI `10.5281/zenodo.24413606`, original predecessor certificate SHA-256 `3be0e4aaa1487e82f35a97a3521abe70b05992888c550c8fdd95b5d50f306dbc`. That source contains 43 norm-one vertices and 903 labelled pair squared distances.

Embed old norm-one points as `q_i=(p_i,0)` and append one perpendicular norm-one point `q44=e40` (one more storage coordinate). Every old squared pair distance is exactly unchanged. Each of the 43 new pairs satisfies `||q44-q_i||²=1+1-2<q44,q_i>=2`. Thus `C(44,2)=946=903+43` unordered labelled squared distances.

## Explicit rational Gram, Euclidean PSD and affine minimum

The original seven canonical unit points yield a rank-three rational Gram block with products `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit points form a rank-four block, with diagonal one and off-diagonal products `-1/4`. Thirty-two additional unit coordinate axes are pairwise orthogonal and orthogonal to both blocks. Thus the rational Gram `G` is PSD, has unit diagonal, and exact `rank(G)=3+4+32=39`.

For the complete labelled squared-distance matrix define `D_ii=0`, `D_ij=2-2G_ij` for all distinct labels. With `H=I-11ᵀ/44`, the uniquely determined centered Euclidean Gram is `B=-(1/2)HDH=HGH`. Exact arithmetic gives `rank(B)=39`. Two original canonical unit points are antipodal, placing the origin in the affine hull. Therefore **39 is the minimum Euclidean affine realization dimension** of this specific complete labelled metric.

## Exact arithmetic invariants at R=1

```text
Labelled unit points                    44
Unordered squared-distance pairs       946
Predecessor pairs conserved            903/903
New squared-distance pairs             43 x 2
Histogram k=4*d²/R²                    4:6 7:2 8:915 9:3 10:10 12:6 13:2 16:2
Sum of all squared pair distances      7615/4
Affine / centered Euclidean Gram rank  39
Coordinate storage dimension           40
```

## Separately stipulated ideal weighted graph

Independently choose a scalar-unit-mass complete graph with symmetric edge weights `w_ij=R²/d_ij²`. The new vertex joins each of the 43 old vertices uniformly with weight 1/2, giving

```text
L39 = [ L38 + I43/2      -1_43/2 ]
      [ -1_43ᵀ/2             43/2 ] .
```

The exact uniform-join characteristic recurrence is `chi44(z)=z(z-22)chi43(z-1/2)/(z-1/2)`. With inherited nonzero seven-point factor `f6(z)=chi_P7(z)/z`, the closed form is `chi44(z)=z(z-22)^33(z-43/2)^4*f6(z-37/2)`. Both exact-rational implementations verify `nullity(L39-22I)=33` and `nullity(L39-(43/2)I)=4`.

These are dimensionless model eigenvalues, not observed hertz.
