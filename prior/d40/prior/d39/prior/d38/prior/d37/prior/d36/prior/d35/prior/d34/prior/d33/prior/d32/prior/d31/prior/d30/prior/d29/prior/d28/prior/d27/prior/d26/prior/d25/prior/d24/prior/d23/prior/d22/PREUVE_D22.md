# Polybrutus D22 — Exact Euclidean metric, Gram and ideal graph spectrum

## 1. Conservative labelled metric

The fixed predecessor has 26 radius-one vectors with 325 exact labelled squared distances. Append one zero coordinate to each and a new orthogonal unit vector, yielding 27 vertices. No old pair changes, and each of 26 added squared pair distances equals 2. Consequently `C(27,2)=351=325+26`. The exact normalized distance histogram for `k=4d²/R²` is `{4:6,7:2,8:320,9:3,10:10,12:6,13:2,16:2}`; the sum of all squared distances is `2855R²/4`.

## 2. Euclidean realizability and exact minimum affine rank

Let `X` be the exact 27-by-23 coordinate matrix, `D` the exact 27-by-27 matrix of **squared** Euclidean distances, and `H=I-11ᵀ/27`. The primary symbolic program verifies exactly, without floating-point approximations:

```
rank([q2-q1,...,q27-q1])=22
rank(B)=22
B=-(1/2)HDH=(HX)(HX)^T.
```

The explicit Gram factorization establishes positive semidefiniteness. A full **labelled** Euclidean squared-distance matrix determines the centered Gram matrix uniquely. Therefore any Euclidean realization of these exact 351 labelled distances has minimum affine dimension **22**. One of 23 storage coordinates is redundant because of the original centered-simplex relation.

## 3. Model-specific weighted Laplacian

For a *chosen mathematical model* of equal scalar masses and complete edge couplings `w_ij=1/d_ij²` at unit radius, the new vertex attaches with uniform weight 1/2:

```
L22 = [ L21+(1/2)I26   -(1/2)1 ]
      [ -(1/2)1ᵀ             13 ].
```

The global uniform zero mode, new differential mode and shifted old zero-sum modes yield an exact recurrence:

```
chi27(z)=z(z-27/2)*chi26(z-1/2)/(z-1/2).
```

The frozen predecessor factorization is `chi26(z)=z(z-13)^15(z-25/2)^4*f6(z-19/2)`, with `f6(z)=chiP7(z)/z`. Substituting gives

```
chi27(z)=z(z-27/2)^16(z-13)^4*f6(z-10).
```

The primary verifier calculates the complete exact 27-by-27 polynomial and proves equality to the recurrence and closed form; `f6(7/2)!=0` and `f6(3)!=0` ensure exact multiplicities **16** at `27/2` and **4** at `13`. The separate rational-metric verifier reconstructs the graph and checks both multiplicities using polynomial derivative tests independently.

## 4. Adversarial controls and limits

Four deliberately false variants are detected: corrupted pair, duplicated point lowering rank, tilted new vector violating perpendicularity, and permuted point labels preserving an unlabelled histogram but not the complete labelled metric.

This is one reproducible specific mathematical geometry and one chosen graph model. Orthogonal Euclidean extensions and uniform graph joins are established techniques, not newly discovered universal mathematics. No physical spacetime dimension, measured vibration or independent human peer review is claimed.
