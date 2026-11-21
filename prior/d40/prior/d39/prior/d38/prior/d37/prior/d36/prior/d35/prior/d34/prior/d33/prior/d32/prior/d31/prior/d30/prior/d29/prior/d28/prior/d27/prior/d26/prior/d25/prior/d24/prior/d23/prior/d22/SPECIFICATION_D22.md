# Polybrutus D22 — Formal exact 27-vertex Euclidean configuration

**Frozen predecessor:** `Topbrutus/polybrutus-d21@v0.1.0`, scientific commit `f17e4d02f7442d18cdafbcbb0b62380c557e4ef2`, DOI `10.5281/zenodo.24409527`.

Take the exact 26 radius-one unit points `p1,...,p26` from the immutable D21 tagged source. They occupy 22 coordinate-storage slots, one of which is redundant because the original centered regular-simplex five-coordinate block imposes one exact linear relation.

Embed `q_i=(p_i,0)` in 23 storage slots for i=1,...,26. Add a perpendicular unit basis point `q27=e23=(0,...,0,1)`. Then for every predecessor point,

```
||q27-q_i||² = ||q27||² + ||q_i||² - 2<q27,q_i> = 2.
```

All 325 prior labelled squared distances are preserved. All 26 new squared distances are 2. Thus `C(27,2)=351=325+26` pairs. Scaling by radius R>0 yields 26 new squared distances equal to `2R²`. The intrinsic affine dimension equals **22**, using 23 storage-coordinate slots.

Normalized squared-distance classes `k=4d²/R²`:

| k | 4 | 7 | 8 | 9 | 10 | 12 | 13 | 16 |
|---|---:|---:|---:|---:|---:|---:|---:|
| pairs | 6 | 2 | 320 | 3 | 10 | 6 | 2 | 2 |

Total squared distance is `2855R²/4`. At R=1 the 23-slot centroid is `(0,0,1/54,0,0,0,0,0,1/27,...,1/27)`, with fifteen final entries 1/27.

## Chosen ideal weighted-graph model

Assign hypothetical unit scalar masses and couplings `w_ij=R²/d_ij²`; these are *defined*, not measured. The new point attaches uniformly to all 26 old points at weight 1/2:

```
L22 = [ L21+(1/2)I26  -(1/2)1 ]
      [ -(1/2)1ᵀ            13  ].
```

With `f6(z)=chiP7(z)/z` the six-nonzero-mode factor of the canonical seven-point graph,

```
chi27(z)=z(z-27/2)*chi26(z-1/2)/(z-1/2)
        =z(z-27/2)^16(z-13)^4*f6(z-10).
```

Eigenvalue `27/2` has exact multiplicity **16** and eigenvalue `13` exact multiplicity **4** *in this specified graph model*. No physical extra dimension or measured resonance follows.
