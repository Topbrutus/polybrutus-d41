# Polybrutus 27D — Exact 32-point Euclidean metric and an ideal graph spectrum

**Scientific author and responsibility:** Gabriel St-Pierre, independent researcher.  
**AI assistance disclosure:** OpenAI ChatGPT (Astra) materially assisted exact mathematics, independent verification programs, manuscript drafting, and publication. AI is not a scientific author.  
**Scientific scope:** Reproducible mathematical case study. No observed extra physical dimension, measured resonance, new universal theorem or independent human peer review is claimed.

## Exact mathematical results, radius R=1

- **32 unit vertices**, **496 labelled squared distances**, and Euclidean minimum affine dimension exactly **27**.
- Conservative extension of the frozen [Polybrutus D26 v0.1.0](https://github.com/Topbrutus/polybrutus-d26/tree/v0.1.0) (source commit `5334c0fae387338314ea6e0fb8450bb3c1a5fcfa`, DOI [10.5281/zenodo.24410491](https://doi.org/10.5281/zenodo.24410491)).
- All **465/465** frozen predecessor labelled squared distances are preserved; one new perpendicular unit vertex adds **31** squared distances, all equal to **2R²**.
- **Eight normalized distance classes** for `k=4d²/R²`: `{4:6,7:2,8:465,9:3,10:10,12:6,13:2,16:2}`; all 496 pairs are accounted for.
- Sum of all unordered squared distances at `R=1`: **4015/4**.
- Exact centered Euclidean Gram rank **27**; positivity proven by explicit rational Gram construction.
- One separately **stipulated ideal complete weighted graph** `w_ij=R²/d_ij²` with unit scalar masses has eigenvalue **16 of exact multiplicity 21** and eigenvalue **31/2 of multiplicity 4**. They are dimensionless model quantities, not measured frequencies.
- Two independently written exact verifiers compare all 496 current labels and 465 frozen pairs, Gram rank, and graph eigenspace dimensions. Four adversarial negative configurations are rejected by the primary program.

## Reproducibility

```bash
python -m pip install sympy==1.14.0
python src/verify_d27.py
python tests/independent_d27.py
```

Both emit the same SHA-256 for the deterministic `results/certificate_d27.json`. All files from the frozen D26 scientific tag were packaged without Git metadata under `prior/d26/`.

## Scientific contents

- [Exact specification](SPECIFICATION_D27.md) and [proof](PREUVE_D27.md).
- [Provenance and limitations](PROVENANCE.md).
- [Scientific note PDF](publication/POLYBRUTUS_27D_NOTE.pdf) and [Markdown source](publication/POLYBRUTUS_27D_NOTE.md).
- [Exact labelled distance certificate](results/certificate_d27.json).
- Frozen parent DOI: https://doi.org/10.5281/zenodo.24410491.

**Licences:** MIT for software, CC BY 4.0 for scientific note.
