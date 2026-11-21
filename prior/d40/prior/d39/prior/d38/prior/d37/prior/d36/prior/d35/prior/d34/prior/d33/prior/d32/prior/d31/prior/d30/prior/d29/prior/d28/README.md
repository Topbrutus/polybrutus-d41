# Polybrutus 28D - Exact 33-point Euclidean metric and chosen ideal graph modes

**Human scientific author:** Gabriel St-Pierre (independent researcher).  
**AI assistance disclosure:** OpenAI ChatGPT (Astra) assisted exact computation, independent verification code, manuscript preparation, and publication; AI is not an author.  
**Scope:** An exact reproducible *mathematical case study*, not evidence for extra physical spacetime dimensions, measured frequencies, a new universal theorem, or completed human peer review.

## Results for unit radius R=1

- **33 distinct unit points**, **528 exact labelled squared pair distances**, and **minimum Euclidean affine dimension 28**.
- Frozen predecessor: [Polybrutus D27, scientific tag v0.1.0](https://github.com/Topbrutus/polybrutus-d27/tree/v0.1.0), source commit `a311d63291d9254664340b434c5a13162d0658c1`, DOI [10.5281/zenodo.24410593](https://doi.org/10.5281/zenodo.24410593).
- **496/496** predecessor labelled squared distances preserved. A perpendicular new unit point contributes **32** additional squared distances, all `2R²`.
- Eight normalized classes `k=4d²/R²`: `{4:6,7:2,8:497,9:3,10:10,12:6,13:2,16:2}`. Sum of all squared pair distances at unit radius: **4271/4**.
- An explicitly constructed exact rational positive-semidefinite Gram matrix has centered rank **28**.
- One **stipulated mathematical complete weighted graph**, unit masses and edge weights `w_ij=R²/d_ij²`, has weighted-Laplacian eigenvalue **33/2 of multiplicity 22**, and eigenvalue **16 of multiplicity 4**. They are dimensionless ideal-model results.
- Two separately written verifiers rederive all 528 labelled distances and compare all 496 predecessor pairs; both independently calculate Gram and graph ranks. Primary verification checks four intentionally invalid cases.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python src/verify_d28.py
python tests/independent_d28.py
```

Both scripts print the identical SHA-256 of the deterministic `results/certificate_d28.json`. The frozen parent scientific tag and its source/certificates are archived as ordinary files under `prior/d27/` (no nested Git submodule).

## Repository evidence

- [Specification](SPECIFICATION_D28.md), [exact proof](PREUVE_D28.md), [provenance](PROVENANCE.md).
- [Scientific manuscript PDF](publication/POLYBRUTUS_28D_NOTE.pdf) and [Markdown source](publication/POLYBRUTUS_28D_NOTE.md).
- [Exact 528 labelled squared distances](results/certificate_d28.json).
- Parent Zenodo DOI: https://doi.org/10.5281/zenodo.24410593.

**Licences:** code MIT; scientific manuscript CC BY 4.0.
