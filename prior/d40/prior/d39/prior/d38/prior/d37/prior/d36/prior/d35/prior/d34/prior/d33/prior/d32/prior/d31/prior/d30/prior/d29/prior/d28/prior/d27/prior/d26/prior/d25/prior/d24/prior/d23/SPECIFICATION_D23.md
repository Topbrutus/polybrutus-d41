# Polybrutus D23 — Exact orthogonal extension specification

## Frozen predecessor

`Topbrutus/polybrutus-d22@v0.1.0` — source commit `4d4e80a7156a118f8f2bebf63aca76ac5f0f06eb`, parent DOI `10.5281/zenodo.24409636`. The frozen 27 unit vectors and all 351 labelled squared distances are immutable.

Embed the old points as `q_i=(p_i,0)` for `1<=i<=27`, and append `q_28=e_24`, a unit vector orthogonal to the entire old 23-coordinate storage space.

For all old points: `||q_28-q_i||²=||q_28||²+||q_i||²-2<q_28,q_i>=2`. Thus every predecessor labelled pair is unchanged and `C(28,2)=378=351+27`.

## Exact constructive Gram model

The seven canonical vectors use the inner product `3*x_i*x_j+y_i*y_j+z_i*z_j` in three orthogonal coordinates. Five centered simplex vertices have `<s_i,s_i>=1` and `<s_i,s_j>=-1/4` for `i!=j`, adding affine rank four. Sixteen remaining points are mutually perpendicular unit axes. All cross-block products vanish. The full 28-by-28 Gram matrix `G` therefore has rank `3+4+16=23`.

Every point has norm one and `D_ij=2-2*G_ij` for `i!=j`, `D_ii=0`. `D` is an **exact squared-distance** matrix. The first and fifth canonical points are antipodal, placing the origin in the affine hull. Hence affine rank is 23. The centered Gram identity `B=-(1/2)H D H=H G H`, with `H=I-11^T/28`, proves Euclidean realizability and positive semidefiniteness.

## Exact histogram and sum

```text
k=4*d²/R² | 4:6  7:2  8:347  9:3  10:10  12:6  13:2  16:2
sum of all 378 squared distances = 3071*R²/4
```

## Ideal graph definition

Define the complete graph weights `w_ij=R²/d_ij²`. At `R=1`, the new point joins each predecessor with weight 1/2:

```text
L23 = [ L22 + I_27/2     -1_27/2 ]
      [ -1_27^T/2         27/2   ]
```

Its characteristic polynomial is `chi28(z)=z*(z-14)*chi27(z-1/2)/(z-1/2)`, equivalently `z*(z-14)^17*(z-27/2)^4*f6(z-21/2)`, where `f6(z)=chi_P7(z)/z` is the canonical seven-vertex graph factor. The chosen graph has eigenvalue 14 exactly 17 times and eigenvalue 27/2 exactly four times.

The spectral model is stipulated and synthetic. No observed physical extra dimension, vibration, or universal new theorem is claimed.
