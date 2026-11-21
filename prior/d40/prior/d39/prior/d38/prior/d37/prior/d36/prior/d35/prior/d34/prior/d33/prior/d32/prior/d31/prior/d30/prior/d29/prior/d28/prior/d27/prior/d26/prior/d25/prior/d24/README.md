# POLYBRUTUS 24D — Exact 29-vertex Euclidean metric and an eighteenfold ideal graph mode

**Human research author and scientific responsibility:** Gabriel St-Pierre, independent researcher.  
**AI disclosure:** OpenAI ChatGPT (Astra) substantially assisted mathematics, exact verification, documentation and publication. AI is not an author.  
**Scientific status:** Specific reproducible mathematical case study. No experimental extra physical dimension, measured vibration or novel universal theorem claimed. No independent human peer review.

## Proven exact mathematical results for R=1

- **29 distinct unit vertices**, **406 labelled squared pairwise distances**, and minimum Euclidean affine embedding dimension exactly **24**.
- The frozen **Polybrutus 23D v0.1.0** parent provides 28 unit points with 378 labelled pairwise squared distances. All **378/378** are preserved, and a new perpendicular unit point contributes **28** squared distances each exactly **2**.
- Eight normalized distance classes `k=4d²/R²`: `{4:6,7:2,8:375,9:3,10:10,12:6,13:2,16:2}`; squared-distance sum **3295/4**.
- Exact positive-semidefinite centered Gram matrix of rank **24**.
- **Specified ideal complete weighted graph**, `w_ij=R²/d_ij²`, with equal unit scalar masses: eigenvalue **29/2 of exact multiplicity 18** and eigenvalue **14 of exact multiplicity 4**. These graph eigenvalues are dimensionless, not measured frequencies.
- Independent exact verification scripts using SymPy and Python `fractions.Fraction`; first script rejects four deliberately wrong variants. Both reproduce the same SHA-256 of the labelled certificate.

## Reproducibility

```bash
python -m pip install sympy==1.14.0
python src/verify_d24.py
python tests/independent_d24.py
```

The original D23 source, both original verifiers and the original D23 labelled certificate remain copied from the immutable tag under `prior/d23/`; predecessors are not modified.

## Source and scientific documents

- [Specification](SPECIFICATION_D24.md), [proof](PREUVE_D24.md), [provenance](PROVENANCE.md)
- [Scientific note PDF](publication/POLYBRUTUS_24D_NOTE.pdf), [note text](publication/POLYBRUTUS_24D_NOTE.md)
- [Exact labelled certificate](results/certificate_d24.json)
- Frozen predecessor [D23 source v0.1.0](https://github.com/Topbrutus/polybrutus-d23/tree/v0.1.0) and [DOI 10.5281/zenodo.24410035](https://doi.org/10.5281/zenodo.24410035)

**Licence:** code MIT, scientific manuscript CC BY 4.0.
