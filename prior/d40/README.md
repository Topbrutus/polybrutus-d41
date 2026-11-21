# Polybrutus 40D — Exact 45-point Euclidean metric and a 34-fold ideal graph mode

**Human scientific author and responsibility:** Gabriel St-Pierre, independent researcher.  
**Material AI assistance:** OpenAI ChatGPT (Astra) assisted exact arithmetic, two separately implemented verification programs, scientific manuscript and controlled publication. AI is not a scientific author.

## Verified exact mathematical properties at unit radius R=1

- **45 norm-one labelled vertices**, **990 unordered labelled squared distances**, minimum affine Euclidean realization dimension **40**.
- **946/946** original predecessor pair labels are unchanged from the frozen [Polybrutus D39 scientific tag v0.1.0](https://github.com/Topbrutus/polybrutus-d39/tree/v0.1.0), scientific commit `a43cd3c67017fb6162005fd521c156cb78795beb`, DOI [10.5281/zenodo.24413712](https://doi.org/10.5281/zenodo.24413712).
- One newly appended perpendicular norm-one vertex contributes **44** new squared distances, each equal to **2R²**.
- Exact normalized squared-distance class histogram `k=4d²/R²`: `{4:6,7:2,8:959,9:3,10:10,12:6,13:2,16:2}`. Exact total squared-distance sum for R=1: **7967/4**.
- The complete exact rational centered Euclidean Gram matrix is positive semidefinite with rank **40**.
- In one **separately stipulated ideal complete weighted graph** with scalar unit masses and symmetric weights `w_ij=R²/d_ij²`, the weighted Laplacian has eigenvalue **45/2** of multiplicity **34**, and **22** of multiplicity **4**. These dimensionless model eigenmodes are **not experimentally measured physical frequencies**.
- Two separately implemented exact-rational programs verify all 990 labels, all 946 predecessor pair values, Euclidean ranks and stipulated graph eigenspaces. Four primary negative controls pass.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python src/verify_d40.py
python tests/independent_d40.py
```

Both verification programs print the same SHA-256 for deterministic `results/certificate_d40.json`. The complete frozen predecessor D39 scientific tag, including its ancestors, is archived under `prior/d39/` as ordinary files via `git archive` (no nested Git repository or submodules).

## Scientific materials and limits

- [Exact specification](SPECIFICATION_D40.md), [constructive proof and adversarial controls](PREUVE_D40.md), [provenance](PROVENANCE.md).
- [Scientific PDF preprint](publication/POLYBRUTUS_40D_NOTE.pdf), [Markdown manuscript](publication/POLYBRUTUS_40D_NOTE.md), [certificate of all 990 labelled distances](results/certificate_d40.json).

This is one exactly reproducible particular synthetic Euclidean metric and one stipulated graph, using known orthogonal Euclidean extension and uniform graph-join techniques. It does **not** claim additional physical spacetime dimensions, measured resonances, particles, a novel universal theorem, or independent human peer review.

**Licenses:** software MIT; scientific manuscript CC BY 4.0.
