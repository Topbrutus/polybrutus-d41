# Polybrutus 41D — Exact 46-point Euclidean metric and a 35-fold ideal graph mode

**Responsible human scientific author:** Gabriel St-Pierre, independent researcher.  
**Material AI assistance:** OpenAI ChatGPT (Astra) assisted exact calculations, two separately implemented checking programs, manuscript and controlled publication. AI is not a scientific author.

## Verified exact mathematics for unit radius R=1

- **46 labelled norm-one vertices**, **1,035 unordered labelled squared-distance pairs**, minimum Euclidean affine realization dimension **41**.
- **990/990** predecessor pair labels are conserved from frozen [Polybrutus D40 v0.1.0](https://github.com/Topbrutus/polybrutus-d40/tree/v0.1.0), scientific commit `2a9fa8265b709d6d971afb6c63b7b024011e67a1`, DOI [10.5281/zenodo.24414238](https://doi.org/10.5281/zenodo.24414238).
- One newly appended perpendicular unit vertex contributes **45** new squared distances, each equal to `2R²`.
- Exact normalized histogram `k=4d²/R²`: `{4:6,7:2,8:1004,9:3,10:10,12:6,13:2,16:2}`; exact sum at R=1: **8327/4**.
- Exact rational positive semidefinite centered Euclidean Gram rank **41**.
- In a **separately stipulated ideal complete weighted graph** with scalar unit masses and edges `w_ij=R²/d_ij²`, its weighted Laplacian has eigenvalue **23** with multiplicity **35** and eigenvalue **45/2** with multiplicity **4**. These are dimensionless *model* eigenmodes, **not measured physical frequencies**.
- Two separately implemented exact-rational programs check all 1,035 labels, 990 frozen parent pair values, Euclidean ranks and selected graph eigenspaces; four principal adversarial controls pass.

## Reproducible verification

```bash
python -m pip install sympy==1.14.0
python src/verify_d41.py
python tests/independent_d41.py
```

The two programs produce an identical SHA-256 for deterministic `results/certificate_d41.json`. The complete frozen D40 source tag and older ancestors are archived as ordinary files under `prior/d40/` using `git archive`, with no nested `.git` or submodules.

## Official DOI and verified public release

- **Official DOI:** [10.5281/zenodo.24414346](https://doi.org/10.5281/zenodo.24414346).
- **Public Zenodo record:** https://zenodo.org/records/24414346, all 3 public files checked for exact sizes and MD5 against the original bytes.
- **Frozen science tag:** [v0.1.0](https://github.com/Topbrutus/polybrutus-d41/tree/v0.1.0), source commit `17cc53498fd96d2e50fbdfa63df6705a4c062299`.
- **Tag CI:** [GitHub Actions 44882862876](https://github.com/Topbrutus/polybrutus-d41/actions/runs/44882862876), success.
- **Durable Memoria checkpoint:** https://github.com/Topbrutus/topbrutus_one/blob/main/memoria/checkpoints/POLYBRUTUS_D41_PUBLICATION_2026-11-21.md .
- **Publications registry:** https://github.com/Topbrutus/topbrutus_one/blob/main/memoria/projects/POLYBRUTUS_X72_PUBLICATIONS_INDEX.md .

## Scientific materials and limitations

- [Exact specification](SPECIFICATION_D41.md), [proof and adversarial tests](PREUVE_D41.md), [provenance](PROVENANCE.md).
- [Two-page PDF scientific preprint](publication/POLYBRUTUS_41D_NOTE.pdf), [Markdown manuscript](publication/POLYBRUTUS_41D_NOTE.md) and [complete 1,035-pair exact certificate](results/certificate_d41.json).

This describes one reproducible *synthetic* Euclidean metric and one separately chosen graph using established perpendicular embedding and uniform graph joining. It does **not** establish physical extra spacetime dimensions, observed particles, measured resonance, a new universal theorem or independent human peer review.

**Licenses:** code MIT; scientific manuscript CC BY 4.0.
