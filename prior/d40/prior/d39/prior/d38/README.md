# Polybrutus 38D — Exact 43-point Euclidean metric and 32-fold ideal graph mode

**Scientific author and responsibility:** Gabriel St-Pierre, independent researcher.  
**Substantive AI assistance:** OpenAI ChatGPT (Astra) assisted exact arithmetic, separate verification programs, manuscript drafting and controlled publication. AI is not the scientific author.

## Reproducible exact results at radius R=1

- **43 norm-one labelled vertices**, **903 unordered labelled squared-distance pairs**, and minimum Euclidean affine realization dimension **38**.
- The immutable parent is [Polybrutus D37 scientific tag v0.1.0](https://github.com/Topbrutus/polybrutus-d37/tree/v0.1.0), scientific commit `294be7f814f2c8770af806833e98194fb192afff`, DOI [10.5281/zenodo.24413346](https://doi.org/10.5281/zenodo.24413346).
- **861/861** predecessor labelled squared distances are exactly conserved; one new perpendicular unit vertex adds **42** squared distances, all **2R²**.
- Exact normalized histogram `k=4d²/R²`: `{4:6,7:2,8:872,9:3,10:10,12:6,13:2,16:2}`. Sum of all 903 squared distances at R=1: **7271/4**.
- The complete labelled metric has exact rational positive-semidefinite centered Gram rank **38**.
- For a separately **stipulated ideal complete weighted graph** with scalar unit masses and symmetric edges `w_ij=R²/d_ij²`, its Laplacian has eigenvalues **43/2** of multiplicity **32** and **21** of multiplicity **4**. These are dimensionless model eigenmodes, not measured physical frequencies.
- Two separately implemented exact-rational programs check every one of 903 pairs, all 861 conserved parent pair labels, Gram ranks and ideal graph eigenspaces. Four adversarial primary checks pass.

## Reproduce the tests

```bash
python -m pip install sympy==1.14.0
python src/verify_d38.py
python tests/independent_d38.py
```

Both programs emit the same SHA-256 of `results/certificate_d38.json`. The frozen scientific tag D37, with its ancestral sources, is archived under `prior/d37/` as ordinary files, without nested Git metadata or submodules.

## Scientific materials

- [Exact specification](SPECIFICATION_D38.md), [constructive proof and negative controls](PREUVE_D38.md), [provenance and limitations](PROVENANCE.md).
- [Scientific PDF manuscript](publication/POLYBRUTUS_38D_NOTE.pdf) and [manuscript Markdown](publication/POLYBRUTUS_38D_NOTE.md).
- [Deterministic certificate of all 903 labelled squared distances](results/certificate_d38.json).

**Limitations:** This is a reproducible particular synthetic mathematical construction using known Euclidean orthogonal extensions and graph joining. It does not claim physical extra spacetime dimensions, experimentally measured frequencies or particles, a novel universal theorem, or independent human peer review.

**Licences:** code MIT; scientific text CC BY 4.0.
