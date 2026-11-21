# Polybrutus 32D — Exact 37-point Euclidean metric and 26-fold ideal graph mode

**Human scientific author:** Gabriel St-Pierre, independent researcher.  
**AI disclosure:** OpenAI ChatGPT (Astra) materially assisted exact mathematics, independent verifying programs, scientific text, and controlled publication; AI is not a scientific author.  
**Scope:** A reproducible, specific synthetic mathematical case study, not evidence for extra physical spacetime dimensions, measured resonances, a new universal theorem, or independent human peer review.

## Exact results for radius R=1

- **37 labelled norm-one vertices**, **666 unordered labelled squared distances**, and **minimum Euclidean affine dimension 32**.
- **630/630** labelled squared distances from frozen [Polybrutus D31 scientific tag v0.1.0](https://github.com/Topbrutus/polybrutus-d31/tree/v0.1.0) are preserved, with scientific commit `0c1f546a6c5442bc83a5bdde506132cce20b6e05` and DOI [10.5281/zenodo.24411160](https://doi.org/10.5281/zenodo.24411160).
- One new perpendicular unit point contributes **36** squared distances all equal to **2R²**.
- Exact normalized classes `k=4d²/R²`: `{4:6,7:2,8:635,9:3,10:10,12:6,13:2,16:2}`; exact sum of squared distances at R=1: **5375/4**.
- Explicit rational positive-semidefinite Gram construction and exact centered Gram rank **32**.
- In a **separately stipulated complete ideal weighted graph**, scalar masses one and `w_ij=R²/d_ij²`, the Laplacian eigenvalue **37/2 has multiplicity 26**, and **18 has multiplicity 4**. These are dimensionless selected-model eigenvalues, not measured physical frequencies.
- Two independently implemented exact-rational programs compare all 666 present labelled pair distances, all 630 frozen predecessor labels, centered affine rank and graph eigenspaces. Four primary adversarial configurations are rejected.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python src/verify_d32.py
python tests/independent_d32.py
```

Both report identical SHA-256 of deterministic `results/certificate_d32.json`. The full frozen parent D31 scientific source, including ancestral sources, is archived under `prior/d31/` as ordinary files without a nested Git submodule.

## Scientific files

- [Exact specification](SPECIFICATION_D32.md), [proof and adversarial checks](PREUVE_D32.md), [provenance and limitations](PROVENANCE.md).
- [Scientific PDF](publication/POLYBRUTUS_32D_NOTE.pdf) and [Markdown manuscript](publication/POLYBRUTUS_32D_NOTE.md).
- [All 666 labelled squared distances](results/certificate_d32.json).

**Licensing:** code MIT; scientific manuscript CC BY 4.0.
