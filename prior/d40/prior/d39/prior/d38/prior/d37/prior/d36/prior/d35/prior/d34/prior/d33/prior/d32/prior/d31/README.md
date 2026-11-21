# Polybrutus 31D — Exact 36-point Euclidean metric and 25-fold ideal graph modes

**Human scientific author and responsibility:** Gabriel St-Pierre, independent researcher.  
**AI disclosure:** OpenAI ChatGPT (Astra) materially assisted exact mathematical verification, separately implemented checking programs, scientific documentation and release. AI is not the scientific author.  
**Scientific scope:** A specific, reproducible mathematical construction; no measured resonance, observed physical extra dimension, novel universal theorem, or independent human peer review is claimed.

## Reproducible mathematical results (radius R=1)

- **36 distinct unit vertices**, with **630** exact labelled squared distances. A centered Euclidean Gram matrix is positive semidefinite with exact rank **31**, establishing minimum affine realization dimension **31**.
- Frozen predecessor: [Polybrutus D30 v0.1.0](https://github.com/Topbrutus/polybrutus-d30/tree/v0.1.0), scientific commit `3c7ce4019a49ff4c54adb9271367aff5c62c1326`, public DOI [10.5281/zenodo.24411018](https://doi.org/10.5281/zenodo.24411018).
- **595/595** labelled squared distances inherited unchanged, and **35** new squared distances, each exactly `2R²`.
- Distance histogram `k=4d²/R²`: `{4:6,7:2,8:599,9:3,10:10,12:6,13:2,16:2}`. Total squared distance for R=1: **5087/4**.
- In one **deliberately stipulated** complete graph with unit scalar masses and edge weights `w_ij=R²/d_ij²`, the weighted graph Laplacian has eigenvalue **18 of exact multiplicity 25** and **35/2 of exact multiplicity 4**. These are dimensionless graph modes, **not physical hertz**.
- Two separately implemented exact-rational verifiers independently check all 630 labelled pairs, all 595 frozen predecessor pairs, both matrix ranks, graph multiplicities, and negative controls.

## Run exact tests

```bash
python -m pip install sympy==1.14.0
python src/verify_d31.py
python tests/independent_d31.py
```

Both programs reproduce `results/certificate_d31.json` and output the same SHA-256. The complete D30 frozen scientific tag and all inherited sources are archived as ordinary files in `prior/d30/`, without nested Git metadata.

## Scientific artifacts

- [Exact specification](SPECIFICATION_D31.md), [proof and controls](PREUVE_D31.md), [provenance and limitations](PROVENANCE.md).
- [Manuscript PDF](publication/POLYBRUTUS_31D_NOTE.pdf), [manuscript source](publication/POLYBRUTUS_31D_NOTE.md).
- [All 630 labelled squared distances](results/certificate_d31.json).

**Licences:** MIT (code), CC BY 4.0 (scientific text).
