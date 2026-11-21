# Polybrutus 38D — Exact labelled Euclidean distance specification

## 1. Frozen predecessor and conservative extension

Use exclusively `Topbrutus/polybrutus-d37@v0.1.0`, scientific commit `294be7f814f2c8770af806833e98194fb192afff`, DOI `10.5281/zenodo.24413346`. The original predecessor certificate SHA-256 is `f649594103c43927126174dc5d807777bb976bf4feae2e9a1c0ee73aa67f3b2b`, covering 42 norm-one vertices and 861 labelled distances.

Embed the previous vertices `q_i=(p_i,0)` and add a new perpendicular norm-one unit vertex `q43=e39`. All predecessor pair labels retain their original squared distances. For each old unit vertex, `||q43-q_i||²=1+1-2<q43,q_i>=2`. Therefore `C(43,2)=903=861+42` unordered labelled squared-distance pairs.

## 2. Constructive rational Gram and minimum affine rank

The seven canonical unit vertices have an exact rank-three rational Gram block with entries `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit vertices form a rank-four Gram block with unit diagonal and off-diagonal inner products `-1/4`. Thirty-one further pairwise mutually orthogonal unit axes are orthogonal to the two blocks. Consequently the full rational Gram `G` is PSD, has unit diagonal, and rank `3+4+31=38`.

Define `D_ii=0` and `D_ij=2-2G_ij` for every distinct labelled pair. The centering projector is `H=I-11ᵀ/43`; the uniquely determined centered Euclidean Gram is `B=-(1/2)HDH=HGH`. Exact rational computation yields `rank(B)=38`. An original antipodal canonical pair places the origin in the affine hull. Thus every Euclidean realization of this **specific complete labelled metric** requires at least 38 affine dimensions, and the explicit construction attains that minimum.

## 3. Exact distance invariants at R=1

```text
Norm-one vertices                       43
Unordered labelled squared distances   903
Frozen predecessor pairs preserved     861/861
New squared-distance pairs             42 × 2
Histogram k=4*d²/R²                    4:6, 7:2, 8:872, 9:3, 10:10, 12:6, 13:2, 16:2
Sum of squared pair distances          7271/4
Affine and centered Gram rank          38
Coordinate storage count               39
```

## 4. Separately stipulated ideal complete weighted graph

Independently stipulate a complete weighted graph with scalar unit masses and edges `w_ij=R²/d_ij²`. The new unit vertex joins every one of the 42 old vertices at weight `1/2`:

```text
L38 = [ L37 + I42/2       -1_42/2 ]
      [ -1_42ᵀ/2               21  ]
```

The exact uniform-join characteristic-polynomial recurrence is `chi43(z)=z(z-43/2)chi42(z-1/2)/(z-1/2)`. Using the inherited nonzero seven-point model factor `f6(z)=chi_P7(z)/z`, the closed expression is `chi43(z)=z(z-43/2)^32(z-21)^4*f6(z-18)`. Both exact-rational programs verify `nullity(L38-(43/2)I)=32` and `nullity(L38-21I)=4`.

These are dimensionless theoretical eigenspace properties of the chosen graph, not observed hertz.
