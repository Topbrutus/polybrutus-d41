# Polybrutus 29D — Exact orthogonal conservation specification

## Frozen predecessor

The sole frozen source is `Topbrutus/polybrutus-d28@v0.1.0`, scientific commit `496e6373abbabfad5c2e99adfcb3994ce8e39531`, DOI `10.5281/zenodo.24410652`. Its original certificate SHA-256 is `d21c21645bc96aa480c840fa24e12e755f6a21d8df8a1d51dac989c7caed1788`: **33** unit points and **528** labelled squared distances.

Embed the 33 old unit points as `q_i=(p_i,0)`, appending a zero coordinate. Create a new perpendicular unit basis point `q34=e30`. All old squared distances remain equal, and for each `i=1,...,33`, `||q34-q_i||²=1+1-2<q34,q_i>=2`. Thus `C(34,2)=561=528+33` labelled pairs.

## Constructive exact Gram matrix

The original seven canonical points span a rank-three block with inner product `3*x_i*x_j+y_i*y_j+z_i*z_j`. Five centered regular-simplex unit points form a rank-four block (self-products 1, distinct-pair products `-1/4`). **Twenty-two** additional mutually perpendicular unit-axis points span 22 independent directions and are orthogonal to both former blocks. Therefore the rational positive-semidefinite Gram matrix `G` has diagonal 1 and `rank(G)=3+4+22=29`.

For `D_ii=0` and `D_ij=2-2G_ij` for distinct points, set `H=I-11ᵀ/34`. The centered Euclidean Gram matrix is `B=-(1/2)HDH=HGH`. Two canonical core points are antipodal, hence the origin lies in the affine hull and the exact rank of `B` is **29**, giving the minimum Euclidean affine embedding dimension of this particular complete labelled distance metric.

## Exact invariant targets

```text
unit vertices                  34
unordered labelled pairs       561
conserved parent pairs         528/528
new squared distances          33 x 2
histogram k=4*d²/R²            4:6, 7:2, 8:530, 9:3, 10:10, 12:6, 13:2, 16:2
sum of squared distances R=1   4535/4
affine rank                    29
centered Gram rank             29
```

## Chosen complete weighted graph

Specify one ideal unit-mass complete weighted graph, `w_ij=R²/d_ij²`. Every old vertex joins the new point with weight 1/2:

```text
L29 = [ L28 + I33/2       -1_33/2 ]
      [ -1_33ᵀ/2            33/2    ]
```

The characteristic-polynomial recurrence is `chi34(z)=z(z-17)chi33(z-1/2)/(z-1/2)`. Equivalently, the inherited canonical factor yields `chi34(z)=z(z-17)^23(z-33/2)^4*f6(z-27/2)`, `f6(z)=chi_P7(z)/z`. Independently computed exact-rational rank gives `nullity(L29-17I)=23` and `nullity(L29-(33/2)I)=4`.

These model eigenvalues are dimensionless and are not observations of vibration, extra physical dimensions or particles.
