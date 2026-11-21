# Polybrutus 26D — Exact 31-point Euclidean metric and ideal weighted-graph modes

**Scientific author and responsibility:** Gabriel St-Pierre, independent researcher.  
**AI assistance:** OpenAI ChatGPT (Astra) materially assisted exact calculations, independent verification, documentation, and controlled release. AI is not a scientific author.  
**Research scope:** A reproducible mathematical case study. No extra physical spacetime dimension, measured resonance, new universal theorem, or independent human peer review is claimed.

## Exact results for unit radius R=1

- **31 distinct unit points**, **465 exact labelled squared distances**, and minimum Euclidean **affine dimension exactly 26**.
- Parent is the immutable [Polybrutus D25 source tag v0.1.0](https://github.com/Topbrutus/polybrutus-d25/tree/v0.1.0), source commit `0dda9ed43622385d3e5fdddef9336458a32ed5ed` and DOI `10.5281/zenodo.24410435`.
- **435/435 labelled squared-distance pairs from D25 are preserved** and the new perpendicular unit point adds exactly **30** squared distances, all **2** at unit radius.
- Normalized class counts `k=4d²/R²`: `{4:6,7:2,8:434,9:3,10:10,12:6,13:2,16:2}`, giving 465 pairs. The sum of squared distances at R=1 is **3767/4**.
- The centered Euclidean squared-distance Gram matrix is positive semidefinite with exact rank **26**, establishing minimal affine realization dimension 26.
- In **one stipulated ideal graph** with unit masses and complete edge weights `w_ij=R²/d_ij²`, the weighted Laplacian has eigenvalue **31/2 of exact multiplicity 20** and **15 of exact multiplicity 4**. These are dimensionless model eigenvalues; not experimentally measured frequencies.
- Two independently coded verifiers reproduce all labels, predecessor pairs, Gram ranks, and spectral modes. The primary verifier rejects four invalid variants, including a relabelling that preserves the unlabeled histogram.

## Reproduction

```bash
python -m pip install sympy==1.14.0
python src/verify_d26.py
python tests/independent_d26.py
```

Both programs print the same SHA-256 of the deterministic `results/certificate_d26.json`. Frozen full D25 tag sources, including the original exact certificate and older provenance files, are copied under `prior/d25/` with **no nested .git directory or Git submodule**.

## Scientific files

- [Exact specification](SPECIFICATION_D26.md) and [proof](PREUVE_D26.md).
- [Provenance and limitations](PROVENANCE.md).
- [Mathematical preprint PDF](publication/POLYBRUTUS_26D_NOTE.pdf) and [manuscript source](publication/POLYBRUTUS_26D_NOTE.md).
- [465 labelled squared distances and model certificate](results/certificate_d26.json).
- Frozen parent DOI: https://doi.org/10.5281/zenodo.24410435.

**Licences:** source code MIT; scientific manuscript CC BY 4.0.
