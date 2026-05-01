[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![License: CC-BY 4.0](https://img.shields.io/badge/Data%20%26%20Manuscript-CC--BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Status](https://img.shields.io/badge/Status-Scaffolding%20%2B%20Pilot-orange)](pilot/calibration_report.md)

# BioEpistasis — Developmental epistasis in single-cell RNA-seq

**Biological sister project to [epistasis-transformer-heads](https://github.com/mool32/epistasis-transformer-heads): does the synthetic-lethal/redundancy regime found in trained transformer heads have a developmental analog in cell differentiation?**

Theodor Spiro | [ORCID 0009-0004-5382-9346](https://orcid.org/0009-0004-5382-9346) | tspiro@vaika.org

📋 **Status:** Scaffolding + active pilot calibration. Pre-registration is **not yet locked**; results below are preliminary scoping work, not findings.
🧬 **Sister project (ML side):** [epistasis-transformer-heads](https://github.com/mool32/epistasis-transformer-heads) — the source of the testable hypothesis (78% of significant top-30 pairs in Pythia 410M show ε_loss > 0)
🧪 **Methodology:** [`methodology/observational_epistasis_limits.md`](methodology/observational_epistasis_limits.md), [`methodology/perturbseq_scoping.md`](methodology/perturbseq_scoping.md)
🔬 **Sign convention:** [`design_notes.md`](design_notes.md) §1 (locked, READ FIRST before interpreting any ε)

---

## What this is

The sister project [epistasis-transformer-heads](https://github.com/mool32/epistasis-transformer-heads) found that **78% of significant top-30 pairs in Pythia 410M show ε_loss > 0** — joint ablation hurts more than additive, synthetic-lethal-like in Costanzo's terminology. This repo asks whether the same regime appears in cell differentiation.

**Direction:** ML observation → biological hypothesis. Use Schiebinger 2019 reprogramming scRNA-seq data (and Norman 2019 perturb-seq as a true-perturbation comparator) to test whether differentiated cells accumulate the synthetic-lethal/redundancy phenotype among their highest-variance genes, paralleling how trained ML heads accumulate the same phenotype.

This is **the second half of the broader DFE-universality program:**

- [Paper 1, arXiv:2604.10571](https://arxiv.org/abs/2604.10571) — DFE shape across AI architectures
- [Paper 2 (LLM track) — functional-differentiation-dfe](https://github.com/mool32/functional-differentiation-dfe) — single-ablation DFE in Pythia 410M
- [Paper 2 (biology track) — clonal-crystallization-aging](https://github.com/mool32/clonal-crystallization-aging) — bone-marrow vs Pythia in (Gini, eff_N) plane
- [epistasis-transformer-heads](https://github.com/mool32/epistasis-transformer-heads) — second-order signature (epistasis) in ML
- **This repo** — second-order signature (epistasis) in biological development

## Status — pilot, not paper

> **Scaffolding only. Pilot first, then pre-registration based on pilot calibration verdict. No pre-reg lock until pilot sanity checks pass.**

The repository at this stage contains:

- The locked sign convention and design rationale
- Methodology scoping notes (observational vs perturbational designs)
- Pilot iterations on Schiebinger / Paul / Norman datasets
- A pilot calibration verdict file
- No claimed findings; results are pilot diagnostics that gate pre-reg lock

## Repository structure

```
├── design_notes.md                       # Final design + locked sign convention
├── methodological_findings.md            # Methodology decisions captured during pilot
├── methodology/
│   ├── observational_epistasis_limits.md # Why scRNA-seq alone may not yield ε reliably
│   └── perturbseq_scoping.md             # Norman 2019 / Replogle perturb-seq feasibility
├── pilot/
│   ├── 01_schiebinger_d8_10pairs.ipynb   # First pilot, 10 pairs at day 8
│   ├── 05_norman_iter5.ipynb             # Norman perturb-seq iteration
│   ├── iter3_paul15/                     # Paul 2015 cohort iteration
│   ├── iter4_synth_null/                 # Synthetic null calibration
│   ├── run_bio_pilot.py                  # Pilot runners (one per iteration)
│   ├── run_bio_pilot_iter3.py
│   ├── run_bio_pilot_iter4.py
│   ├── build_pilot_notebook.py           # Notebook builders
│   ├── build_norman_iter5_notebook.py
│   ├── calibration_report.md             # Pilot calibration writeup (in progress)
│   └── pilot_calibration_verdict.json    # Machine-readable verdict
├── README.md
└── LICENSE
```

## Foundational discipline (locked)

- **Sign convention** for ε is explicit in [`design_notes.md`](design_notes.md) §1. **ε_loss > 0 ↔ ε_fitness < 0 ↔ synthetic-lethal/redundancy.** Confused conventions caused interpretation errors in the sister Epistasis project that propagated for several rounds before being caught. Do not skip this section.
- **Pilot calibration verdict precedes pre-reg lock.** No analysis claiming biological findings is run until the pilot demonstrates the apparatus measures something stronger than noise.
- **Pre-reg lock precedes any decision-defining analysis.** Versioned and locked alongside the working draft, matching the practice in `analyses/*.LOCKED.md` of the sister ML project.
- **All numerical results bootstrap-quantified.** No point estimates.

## Citation

The work is at the pilot stage; there is no preprint to cite. If referencing the design or sign convention from `design_notes.md`, please cite the sister ML preprint that motivates the hypothesis:

```bibtex
@misc{spiro2026epistasisheads,
  author = {Spiro, Theodor},
  title  = {Epistasis mapping in transformer attention heads: cross-model multi-checkpoint ablation interactions},
  year   = {2026},
  url    = {https://github.com/mool32/epistasis-transformer-heads}
}
```

## Contact

Theodor Spiro — tspiro@vaika.org

## License

- **Code** (`pilot/*.py`, `pilot/*.ipynb`): MIT (see [LICENSE](LICENSE))
- **Data** (none committed at this stage): N/A
- **Design notes** (`design_notes.md`, `methodology/*.md`, `methodological_findings.md`): CC-BY 4.0
