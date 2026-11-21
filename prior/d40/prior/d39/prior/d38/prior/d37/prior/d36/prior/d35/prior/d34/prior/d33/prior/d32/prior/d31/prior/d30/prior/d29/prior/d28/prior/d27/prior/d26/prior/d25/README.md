# Polybrutus 25D — Exact 30-point Euclidean metric and chosen ideal graph spectrum

**Human author and scientific responsibility:** Gabriel St-Pierre, independent researcher.  
**AI assistance:** OpenAI ChatGPT (Astra) contributed exact symbolic calculations, independently implemented verification, research writing and release preparation; AI is not an author.  
**Scientific scope:** Reproducible mathematical case study, not a discovery of extra physical dimensions, a measured resonance, or a new universal theorem. No independent human peer review has occurred.

## Exact results at unit radius

- **30 unit points**, `C(30,2)=435` exact labelled squared distances, and **Euclidean affine dimension exactly 25**.
- Frozen predecessor [Polybrutus 24D v0.1.0](https://github.com/Topbrutus/polybrutus-d24/tree/v0.1.0) contributes **406/406 preserved labelled squared distances**; the perpendicular new unit point adds **29 squared distances equal to 2**.
- Eight classes for `k=4d²/R²`: `{4:6,7:2,8:404,9:3,10:10,12:6,13:2,16:2}`. Sum of all squared distances at `R=1`: **3527/4**.
- Rational Gram construction, positive semidefinite centered Gram matrix of exact rank **25**.
- **Explicitly chosen ideal unit-mass weighted complete graph**, `w_ij=R²/d_ij²`: eigenvalue **15 with multiplicity 19** and **29/2 with multiplicity 4**. These are dimensionless graph eigenvalues, not measured frequencies.
- Two separately implemented verification programs agree on all labelled distances, predecessor certificate, exact Gram rank and Laplacian multiplicities. The primary script additionally rejects four invalid perturbations.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python src/verify_d25.py
python tests/independent_d25.py
```

Both scripts print an identical SHA-256 of `results/certificate_d25.json`. The frozen predecessor source (including D24 certificate and its verification programs) is copied under `prior/d24/`.

## Scientific contents

- [Specification](SPECIFICATION_D25.md) and [proof](PREUVE_D25.md)
- [Provenance and limitations](PROVENANCE.md)
- [Scientific manuscript PDF](publication/POLYBRUTUS_25D_NOTE.pdf) and [source Markdown](publication/POLYBRUTUS_25D_NOTE.md)
- [All 435 exact labelled squared distances](results/certificate_d25.json)
- Frozen parent DOI: https://doi.org/10.5281/zenodo.24410102

**Licensing:** MIT source code; CC BY 4.0 scientific note.
