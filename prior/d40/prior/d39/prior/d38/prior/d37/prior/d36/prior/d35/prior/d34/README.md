# Polybrutus 34D — Exact 39-point Euclidean metric and 28-fold ideal graph mode

**Human scientific author and responsibility:** Gabriel St-Pierre, independent researcher.  
**AI assistance:** OpenAI ChatGPT (Astra) materially assisted exact calculations, two separately implemented verification programs, scientific writing, and publication. AI is not the scientific author.  
**Scientific status:** Reproducible, particular synthetic mathematical construction; not evidence of physical spacetime dimensions, measured frequencies, a new universal theorem or independent human peer review.

## Reproducible exact results (unit radius R=1)

- **39 labelled unit points** and **741 unordered labelled squared distances**, with minimum Euclidean affine realization dimension **34**.
- All **703/703** predecessor labelled squared distances conserved from frozen [Polybrutus D33 v0.1.0](https://github.com/Topbrutus/polybrutus-d33/tree/v0.1.0), scientific commit `bc78e6c6fe1755a069d8860a0144cbebbb1168a6`, Zenodo DOI [10.5281/zenodo.24411569](https://doi.org/10.5281/zenodo.24411569).
- A single new perpendicular unit point adds **38** squared distances exactly equal to `2R²`.
- Exact histogram of normalized squared-distance classes `k=4d²/R²`: `{4:6,7:2,8:710,9:3,10:10,12:6,13:2,16:2}`. Sum of all 741 squared pair distances at R=1: **5975/4**.
- Explicit rational positive-semidefinite Gram matrix with centered Gram rank **34**.
- For one independently **stipulated mathematical complete weighted graph** with unit masses and `w_ij=R²/d_ij²`, the Laplacian has eigenvalue **39/2 of exact multiplicity 28** and **19 of multiplicity 4**. These are dimensionless model-specific eigenvalues, not measured hertz.
- The primary and independent exact-rational verification programs check every one of 741 pair labels and all 703 frozen predecessor pairs, centered ranks and graph eigenspaces. Four adversarial primary controls are rejected.

## Reproduce the science

```bash
python -m pip install sympy==1.14.0
python src/verify_d34.py
python tests/independent_d34.py
```

Both programs produce identical SHA-256 of the deterministic `results/certificate_d34.json`. All original frozen D33 tag source files, including ancestors, are preserved under `prior/d33/` as ordinary files without nested Git metadata.

## Scientific resources

- [Formal specification](SPECIFICATION_D34.md), [exact proof and negative controls](PREUVE_D34.md), [provenance](PROVENANCE.md).
- [Scientific preprint PDF](publication/POLYBRUTUS_34D_NOTE.pdf), [editable manuscript Markdown](publication/POLYBRUTUS_34D_NOTE.md).
- [All 741 exact labelled squared distances](results/certificate_d34.json).

**Licensing:** code MIT; scientific manuscript CC BY 4.0.
