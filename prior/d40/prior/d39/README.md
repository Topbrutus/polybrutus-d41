# Polybrutus 39D — Exact 44-point Euclidean metric and a 33-fold ideal graph mode

**Human scientific author and responsibility:** Gabriel St-Pierre, independent researcher.  
**Material AI assistance:** OpenAI ChatGPT (Astra) assisted exact calculations, two separately implemented verification programs, manuscript drafting, and publication. AI is not a scientific author.

## Verified exact results at unit radius R=1

- **44 norm-one labelled points**, **946 unordered labelled squared distances**, minimum affine Euclidean realization dimension **39**.
- The **903/903** predecessor pair labels are unchanged from frozen [Polybrutus D38 v0.1.0](https://github.com/Topbrutus/polybrutus-d38/tree/v0.1.0), scientific commit `ac8fd112fcd314092f106b210162ba1cc9879e9b`, DOI [10.5281/zenodo.24413606](https://doi.org/10.5281/zenodo.24413606).
- One new perpendicular unit point adds **43** squared distances each equal to `2R²`.
- Squared distance histogram `k=4d²/R²`: `{4:6,7:2,8:915,9:3,10:10,12:6,13:2,16:2}`; sum of all squared pair distances at R=1: **7615/4**.
- Exact rational positive-semidefinite centered Gram matrix with rank **39**.
- For a separately **stipulated dimensionless ideal complete weighted graph** with scalar unit masses and `w_ij=R²/d_ij²`, graph Laplacian eigenvalue **22** has multiplicity **33**, while **43/2** has multiplicity **4**. These modes are **not physical resonances or measured frequencies**.
- Two separately implemented exact-rational programs verify 946/946 present pairs, 903/903 frozen predecessor labels, Gram ranks and ideal graph eigenspaces. Four adversarial primary controls pass.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python src/verify_d39.py
python tests/independent_d39.py
```

Both scripts print the same SHA-256 for deterministic `results/certificate_d39.json`. The complete frozen D38 scientific tag and its ancestors are preserved as ordinary files under `prior/d38/`, via `git archive`; no nested `.git` or submodules.

## Scientific resources and limitations

- [Exact specification](SPECIFICATION_D39.md), [constructive proof and adversarial controls](PREUVE_D39.md), [provenance](PROVENANCE.md).
- [Scientific PDF preprint](publication/POLYBRUTUS_39D_NOTE.pdf), [editable Markdown manuscript](publication/POLYBRUTUS_39D_NOTE.md).
- [Certificate listing all 946 labelled squared-distance pairs](results/certificate_d39.json).

**Scope:** A specific exactly reproducible synthetic Euclidean labelled metric and separately stipulated graph model, using established orthogonal embedding and uniform graph joining. No new physical spacetime dimensions, observed particles, experimental resonances, new universal theorem, or independent human peer review claimed.

**Licensing:** code MIT, scientific manuscript CC BY 4.0.
