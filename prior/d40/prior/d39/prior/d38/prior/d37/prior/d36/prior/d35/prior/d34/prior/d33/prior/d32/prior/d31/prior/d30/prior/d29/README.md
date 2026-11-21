# Polybrutus 29D — Exact 34-vertex Euclidean labelled metric and ideal weighted-graph modes

**Human scientific author:** Gabriel St-Pierre, independent researcher.  
**AI assistance:** OpenAI ChatGPT (Astra) assisted exact symbolic mathematics, independent checking programs, manuscript drafting and controlled release; AI is not a scientific author.  
**Research status:** A reproducible, specific mathematical case study. No measured physical frequencies, added physical spacetime dimensions, new universal theorem or independent human peer review is claimed.

## Exact results, unit radius R=1

- **34 distinct unit points** and **561 labelled squared distances** with Euclidean minimum affine dimension exactly **29**.
- Frozen predecessor: [Polybrutus 28D v0.1.0](https://github.com/Topbrutus/polybrutus-d28/tree/v0.1.0), scientific source commit `496e6373abbabfad5c2e99adfcb3994ce8e39531`, DOI [10.5281/zenodo.24410652](https://doi.org/10.5281/zenodo.24410652).
- Exactly **528/528 parent labelled squared distances preserved**, and one perpendicular unit point adds **33** new squared distances, each `2R²`.
- Eight normalized squared-distance classes `k=4d²/R²`: `{4:6,7:2,8:530,9:3,10:10,12:6,13:2,16:2}`. They total 561.
- Sum of squared distances for unit radius: **4535/4**.
- Explicit exact-rational Gram representation with positive-semidefinite centered rank **29**.
- A deliberately **stipulated ideal complete weighted graph** with unit scalar masses and `w_ij=R²/d_ij²` has Laplacian eigenvalue **17 of multiplicity 23** and **33/2 of multiplicity 4**, both checked by independent exact-rational rank computations. These are dimensionless chosen-model modes.
- Two independently implemented scripts verify the complete labelled metrics and ranks, including 528/528 conserved parent labels. The primary checker rejects four adversarial variants.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python src/verify_d29.py
python tests/independent_d29.py
```

The programs produce the identical SHA-256 for `results/certificate_d29.json`. The immutable predecessor tag sources and ancestral evidence are archived under `prior/d28/` without any embedded `.git` metadata or Git submodule.

## Scientific files

- [Specification](SPECIFICATION_D29.md), [exact proof](PREUVE_D29.md), [provenance](PROVENANCE.md).
- [Scientific manuscript PDF](publication/POLYBRUTUS_29D_NOTE.pdf), [editable manuscript Markdown](publication/POLYBRUTUS_29D_NOTE.md).
- [Exact 561 labelled squared distances](results/certificate_d29.json).
- Frozen parent DOI: https://doi.org/10.5281/zenodo.24410652.

**Licence:** source software MIT; scientific manuscript CC BY 4.0.
