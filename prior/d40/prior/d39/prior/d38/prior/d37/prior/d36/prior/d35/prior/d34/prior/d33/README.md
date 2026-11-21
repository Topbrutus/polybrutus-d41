# Polybrutus 33D — Exact 38-point Euclidean metric and 27-fold ideal graph mode

**Human scientific author:** Gabriel St-Pierre, independent researcher.  
**AI assistance disclosure:** OpenAI ChatGPT (Astra) materially assisted exact mathematics, independent checking programs, the scientific manuscript and publication. AI is not a scientific author.  
**Scope:** A particular reproducible synthetic mathematical case study; no experimentally observed extra physical spacetime dimensions, measured resonances, new universal theorem, or independent human peer review is claimed.

## Exact results at radius R = 1

- **38 distinct norm-one vertices**, **703 unordered labelled squared-distance pairs**, and minimum Euclidean affine dimension **33**.
- Frozen predecessor: [Polybrutus D32 scientific tag v0.1.0](https://github.com/Topbrutus/polybrutus-d32/tree/v0.1.0), scientific commit `cb63ed7ab447b13e251d92571fe0f1062db31cca`, public DOI [10.5281/zenodo.24411375](https://doi.org/10.5281/zenodo.24411375).
- **666/666 predecessor labelled squared distances preserved** exactly, and a new perpendicular unit point adds **37 squared distances**, all equal to **2R²**.
- Distance classes `k=4d²/R²`: `{4:6,7:2,8:672,9:3,10:10,12:6,13:2,16:2}`. Exact sum of the 703 squared pair distances at R=1: **5671/4**.
- Exact rational positive-semidefinite centered Euclidean Gram rank: **33**.
- In a **separately stipulated ideal complete weighted graph** with unit scalar masses and edges `w_ij=R²/d_ij²`, the weighted-Laplacian eigenvalue **19 has exact multiplicity 27**, and **37/2 has exact multiplicity 4**. These are dimensionless chosen-model modes, not experimentally measured hertz.
- Two independently implemented exact-rational programs check all **703/703** pair labels, all **666/666** conserved predecessor pairs, Gram ranks and graph eigenspaces. Four controlled incorrect cases are rejected by the primary verifier.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python src/verify_d33.py
python tests/independent_d33.py
```

Both programs print an identical SHA-256 for the deterministic `results/certificate_d33.json`. The full frozen D32 scientific tag, with inherited ancestral records, is preserved as ordinary files under `prior/d32/`; no nested Git metadata or submodule.

## Scientific files

- [Exact specification](SPECIFICATION_D33.md), [proof and controls](PREUVE_D33.md), [provenance and limitations](PROVENANCE.md).
- [Scientific PDF](publication/POLYBRUTUS_33D_NOTE.pdf) and [manuscript source](publication/POLYBRUTUS_33D_NOTE.md).
- [All 703 labelled squared distances](results/certificate_d33.json).

**Licences:** software MIT; manuscript CC BY 4.0.
