# BioEpistasis — Developmental epistasis in single-cell RNA-seq

Sister project to `DFE research/Epistasis/` (transformer attention-head
epistasis). Tests whether the synthetic-lethal/redundancy regime
identified in trained transformer heads — 78 % of significant top-30
pairs in Pythia 410M show ε_loss > 0 (joint hurts more than additive,
synthetic-lethal-like in Costanzo's terminology) — has a
developmental analog in cell differentiation.

**Direction:** ML observation → biological hypothesis. Use Schiebinger
2019 reprogramming scRNA-seq data to test whether differentiated cells
accumulate synthetic-lethal/redundancy phenotype among their highest-
variance genes, paralleling how trained ML heads accumulate the same
phenotype.

## Status

Scaffolding only. **Pilot first**, then pre-registration based on pilot
calibration verdict. No data commits, no pre-reg lock until pilot
sanity checks pass.

## Layout

```
BioEpistasis/
├── README.md                     this file
├── design_notes.md               sign convention + final design
├── pilot/
│   ├── 01_schiebinger_d8_10pairs.ipynb
│   └── calibration_report.md     after pilot
└── analyses/
    └── pre_registration_v1.md    after pilot calibration PASS
```

## Foundational discipline (locked)

- **Sign convention** for ε is explicit in `design_notes.md` section 1.
  ε_loss > 0 ↔ ε_fitness < 0 ↔ synthetic-lethal/redundancy. Confused
  conventions caused interpretation errors in the sister Epistasis
  project that propagated for several rounds before being caught.
  Do NOT skip this section.
- Pilot calibration verdict precedes pre-reg lock.
- Pre-reg lock precedes any decision-defining analysis.
- All numerical results bootstrap-quantified.
