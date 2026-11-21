# POLYBRUTUS 23D — Exact 28-point Euclidean metric and an ideal graph mode

**Scientific responsibility and human author:** Gabriel St-Pierre (independent researcher).  
**AI disclosure:** OpenAI ChatGPT (Astra) assisted symbolic computation, independent verification scripts, documentation, and release orchestration. AI is not a scientific author.  
**Status:** Exact mathematical case study, not externally human-peer-reviewed, and no physical extra dimension or measured resonance is claimed.

## Exact results for unit radius R=1

- **28 labelled unit vertices**, **378 exact pairwise squared distances**, and **affine embedding dimension exactly 23**.
- The immutable [Polybrutus D22 v0.1.0](https://github.com/Topbrutus/polybrutus-d22/tree/v0.1.0), source commit `4d4e80a7156a118f8f2bebf63aca76ac5f0f06eb`, supplies all **351/351 preserved predecessor labelled squared distances**.
- One new perpendicular unit vertex creates **27** new squared distances, all exactly **2**.
- Normalized classes `k=4*d²/R²`: `{4:6,7:2,8:347,9:3,10:10,12:6,13:2,16:2}`; squared distance sum **3071/4**.
- The centered Euclidean squared-distance Gram matrix is positive semidefinite by exact construction and has rank **23**.
- **Chosen theoretical graph**, equal unit masses with edge weights `w_ij=R²/d_ij²`: Laplacian eigenvalue **14 of exact multiplicity 17**, eigenvalue **27/2 of exact multiplicity 4**; these are dimensionless graph numbers, not observed frequencies.
- Two separately written programs verify the 378 labelled rational distances, the 351 predecessor distances, Gram rank, and spectral multiplicities. Primary negative controls: **4/4 PASS**.

## Reproduction

```bash
python -m pip install sympy==1.14.0
python src/verify_d23.py
python tests/independent_d23.py
```

Both scripts print the same SHA-256 of `results/certificate_d23.json`, with no floating-point approximations. GitHub Actions runs both. Parent sources and original D22 certificate are preserved under `prior/d22/` without editing the parent release.

## Publication and evidence

- Specification: [SPECIFICATION_D23.md](SPECIFICATION_D23.md)
- Proof: [PREUVE_D23.md](PREUVE_D23.md)
- Provenance and scientific limits: [PROVENANCE.md](PROVENANCE.md)
- Scientific paper: [publication/POLYBRUTUS_23D_NOTE.pdf](publication/POLYBRUTUS_23D_NOTE.pdf)
- Labelled distance certificate: [results/certificate_d23.json](results/certificate_d23.json)
- Frozen parent DOI: https://doi.org/10.5281/zenodo.24409636

**Licences:** MIT for code; CC BY 4.0 for scientific manuscript. All source citations refer to original canonical Polybrutus where appropriate. No novel universal theorem is asserted.
