# Polybrutus 35D — Exact 40-point Euclidean metric and a twenty-ninefold ideal graph mode

**Human scientific author and responsibility:** Gabriel St-Pierre, independent researcher.  
**Substantive AI assistance:** OpenAI ChatGPT (Astra) assisted exact arithmetic, distinct verifier implementations, manuscript and release. AI is not the scientific author.

## Verified exact results (unit radius R=1)

- **40 unit vertices; 780 unordered labelled squared distances; minimum affine Euclidean dimension 35.**
- **741/741** predecessor labelled squared distances conserved from the frozen [Polybrutus D34 v0.1.0 scientific tag](https://github.com/Topbrutus/polybrutus-d34/tree/v0.1.0), commit `21385a5875d472af8f035372096a46f5501cee1a`, DOI [10.5281/zenodo.24412288](https://doi.org/10.5281/zenodo.24412288).
- Exactly **39** new squared distances, each `2R²`.
- Normalized squared-distance classes `k=4d²/R²`: `{4:6,7:2,8:749,9:3,10:10,12:6,13:2,16:2}`. Exact sum of squared distances: **6287/4**.
- Exact rational positive-semidefinite centered Euclidean Gram rank **35**.
- In a separately **stipulated** complete weighted graph with scalar unit masses and `w_ij=R²/d_ij²`, weighted Laplacian eigenvalues **20 (multiplicity 29)** and **39/2 (multiplicity 4)**. These are dimensionless graph-model properties, not measured physical frequencies.
- Two independently implemented exact-rational verification programs check 780/780 present labels, 741/741 original labels, Gram ranks and model eigenspaces. The primary program runs four negative controls.

## Reproduce

```bash
python -m pip install sympy==1.14.0
python src/verify_d35.py
python tests/independent_d35.py
```

Both verifiers produce the same SHA-256 for `results/certificate_d35.json`. All files from the immutable D34 scientific source and older ancestors are copied into `prior/d34/` as ordinary files via `git archive`, without a nested Git repository or gitlinks.

## Scientific source and limits

- [Mathematical specification](SPECIFICATION_D35.md); [constructive proof and negative tests](PREUVE_D35.md); [provenance](PROVENANCE.md).
- [Scientific PDF manuscript](publication/POLYBRUTUS_35D_NOTE.pdf); [editable manuscript](publication/POLYBRUTUS_35D_NOTE.md).
- [Deterministic certificate listing all 780 labelled distances](results/certificate_d35.json).

**Scope:** This is a reproducible specific synthetic mathematical construction using established orthogonal extension and graph-join techniques. It is **not** evidence for additional physical spacetime dimensions, measured resonances, new particles, a universal new theorem, or independent human peer review.

**Licensing:** code MIT; scientific text CC BY 4.0.
