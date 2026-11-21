# Polybrutus 36D — Exact 41-point Euclidean metric and a 30-fold ideal graph mode

**Scientific author and responsibility:** Gabriel St-Pierre, independent researcher.  
**Material AI assistance:** OpenAI ChatGPT (Astra) assisted exact mathematics, two separately implemented checkers, manuscript writing, and release operations. AI is not a scientific author.

## Exact mathematical results at unit radius R = 1

- **41 labelled unit vertices, 820 unordered squared-distance pairs, minimum affine Euclidean realization dimension 36**.
- Frozen predecessor: [Polybrutus D35 scientific tag v0.1.0](https://github.com/Topbrutus/polybrutus-d35/tree/v0.1.0), scientific commit `72989a1fd2a7fb6e28522a0f2dcbe168924ca78d`, DOI [10.5281/zenodo.24412435](https://doi.org/10.5281/zenodo.24412435).
- All **780/780** predecessor labelled squared distances conserved; a single perpendicular unit vertex adds **40** squared distances equal to **2R²**.
- Exact normalized histogram `k=4d²/R²`: `{4:6,7:2,8:789,9:3,10:10,12:6,13:2,16:2}`; sum of all squared pair distances for R=1: **6607/4**.
- The exact rational centered Euclidean Gram is positive semidefinite of rank **36**.
- For a separately **stipulated ideal complete weighted graph** with scalar unit masses and `w_ij=R²/d_ij²`, the weighted-Laplacian eigenvalue **41/2** has multiplicity **30** and eigenvalue **20** has multiplicity **4**. These dimensionless graph modes are **not measured physical frequencies**.
- Two separately implemented exact-rational verification programs check every one of 820 labels and 780 frozen predecessor values, ranks and eigenspaces. Four primary negative controls pass.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python src/verify_d36.py
python tests/independent_d36.py
```

Both programs compute an identical SHA-256 for `results/certificate_d36.json`. The complete frozen D35 scientific tag including older ancestors is archived as ordinary files under `prior/d35/`; no nested `.git` or submodule.

## Scientific artifacts and limits

- [Exact specification](SPECIFICATION_D36.md), [proof and adversarial checks](PREUVE_D36.md), [scientific provenance](PROVENANCE.md).
- [Scientific PDF preprint](publication/POLYBRUTUS_36D_NOTE.pdf), [Markdown manuscript](publication/POLYBRUTUS_36D_NOTE.md), [certificate of all 820 labels](results/certificate_d36.json).

This describes a particular exact synthetic Euclidean metric and separately chosen graph, using established mathematical techniques. It does **not** claim experimentally verified extra physical spacetime dimensions, measured resonances, a novel universal theorem, observed particles or independent human peer review.

**Licences:** MIT (code), CC BY 4.0 (scientific manuscript).
