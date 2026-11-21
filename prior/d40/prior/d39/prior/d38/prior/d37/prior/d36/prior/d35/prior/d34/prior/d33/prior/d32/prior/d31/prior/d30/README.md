# Polybrutus 30D — Exact 35-point Euclidean metric and selected ideal graph modes

**Human scientific author and responsibility:** Gabriel St-Pierre, independent researcher.  
**AI assistance disclosure:** OpenAI ChatGPT (Astra) materially assisted exact calculations, two independent verification scripts, documentation and publication. AI is not a scientific author.  
**Scientific status:** Reproducible mathematical case study, not independent human peer review, an experimentally measured frequency or a claim of extra physical spacetime dimensions.

## Exact mathematical results (radius R = 1)

- **35 labelled unit vertices**, **C(35,2) = 595 exact squared pair distances**, and minimum Euclidean affine dimension **30**.
- **561/561 predecessor labelled squared distances conserved**, from the immutable [Polybrutus D29 scientific tag v0.1.0](https://github.com/Topbrutus/polybrutus-d29/tree/v0.1.0), commit `cbe179697cc1a287fc0637912ea6a6fd890eeca4`, DOI [10.5281/zenodo.24410725](https://doi.org/10.5281/zenodo.24410725).
- A new perpendicular unit point adds **34 squared distances**, each **2 R²**.
- Eight classes `k=4d²/R²`: `{4:6,7:2,8:564,9:3,10:10,12:6,13:2,16:2}`. Sum of all squared distances at R=1: **4807/4**.
- Explicit rational Gram construction; centered Euclidean Gram is **positive semidefinite of exact rank 30**.
- For a **separately stipulated** complete ideal weighted graph with unit masses and edges `w_ij=R²/d_ij²`, the Laplacian has eigenvalue **35/2 of exact multiplicity 24** and **17 of exact multiplicity 4**. These are dimensionless, model-dependent eigenvalues, not experimentally measured resonances.
- Two independently implemented exact-rational verifiers check all 595 present pair labels, all 561 predecessor pair labels, the centered rank and the graph eigenspace dimensions. The primary verifier checks four adversarial incorrect cases.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python src/verify_d30.py
python tests/independent_d30.py
```

Both scripts print the same SHA-256 for the deterministic `results/certificate_d30.json`. The **entire frozen D29 scientific tag**, including ancestor evidence, is archived under `prior/d29/` as ordinary files, without nested Git metadata.

## Source and proof

- [Mathematical specification](SPECIFICATION_D30.md), [proof](PREUVE_D30.md), and [provenance](PROVENANCE.md).
- [PDF scientific preprint](publication/POLYBRUTUS_30D_NOTE.pdf) and [editable Markdown source](publication/POLYBRUTUS_30D_NOTE.md).
- [All 595 labelled squared distances](results/certificate_d30.json).
- Frozen parent DOI: https://doi.org/10.5281/zenodo.24410725.

**Licences:** code MIT; scientific preprint CC BY 4.0.
