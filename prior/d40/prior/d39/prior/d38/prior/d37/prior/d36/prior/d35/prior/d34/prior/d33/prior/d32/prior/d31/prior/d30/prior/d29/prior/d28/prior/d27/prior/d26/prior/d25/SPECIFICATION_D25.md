# Polybrutus D25 — Reproducible exact orthogonal extension

## Frozen predecessor

Use **only** `Topbrutus/polybrutus-d24@v0.1.0`, scientific source commit `16d9f0ee705d15c220bc02e18563a5e5c5c1e5fc`, DOI `10.5281/zenodo.24410102`, with original certificate SHA-256 `c2b1bbfd2c761871e81748fc5621c700443c7a7e7c9daba479d9ff18112c22eb`. Every one of the **406** old labelled squared distances must be preserved exactly.

Embed the 29 old unit points as `q_i=(p_i,0)` and append one new perpendicular unit basis point `q30=e26`. Since every old point has norm one, `||q30-q_i||²=2` for all `1<=i<=29`. Then `C(30,2)=435=406+29` labelled distances.

## Explicit Euclidean Gram construction

Seven canonical points have inner products `3*x_i*x_j + y_i*y_j + z_i*z_j` (three independent directions). Five centered regular-simplex unit vectors satisfy `<s_i,s_j>=1` if `i=j` and `-1/4` otherwise, giving rank four. The remaining **18** unit vectors occupy mutually orthogonal axes and are orthogonal to both prior blocks. Thus Gram rank is exactly `3+4+18=25` and each Gram diagonal entry equals one.

For `D_ij=2-2G_ij` for `i!=j` and `D_ii=0`, the centered Euclidean Gram matrix obeys `B=-(1/2) H D H=H G H`, `H=I-11^T/30`. The first and fifth canonical points are antipodal; therefore the origin is in their affine hull. The full squared-distance matrix has minimal Euclidean affine realization dimension exactly 25.

## Exact invariant targets, verified

- `n=30`, `C(30,2)=435`; all **406** frozen predecessor labels equal their original values.
- **29** new squared distances `2R²`.
- Counts of `k=4d²/R²`: `{4:6,7:2,8:404,9:3,10:10,12:6,13:2,16:2}`.
- Total squared-distance sum at `R=1`: `3527/4`.
- Gram and centered Gram exact rank **25**, positive semidefinite.

## Chosen ideal graph spectrum

For a **mathematically stipulated** complete graph with `w_ij=R²/d_ij²` and unit scalar masses, the new point joins each previous point uniformly with weight `1/2`. Thus:

```text
L25 = [L24 + I29/2    -1_29/2]
      [-1_29^T/2          29/2 ]
```

The uniform-join characteristic recurrence is `chi30(z)=z(z-15)chi29(z-1/2)/(z-1/2)`, equivalently `chi30(z)=z(z-15)^19(z-29/2)^4*f6(z-23/2)`, `f6(z)=chi_P7(z)/z`. Two rational-rank verifiers establish nullities 19 and 4 for `L25-15I` and `L25-(29/2)I` respectively. No laboratory physical quantities are inferred.
