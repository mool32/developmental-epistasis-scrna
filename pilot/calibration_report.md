# BioEpistasis pilot — calibration report (2026-04-26)

## TL;DR

**Verdict: FAIL** (sanity 1 calibration). Pilot identified two
methodological issues that block the originally-designed 2×2 tertile
stratification on sparse marker-gene scRNA-seq:

1. **Tertile cuts are degenerate on sparse, categorical hormone HVG.**
   With > 33 % cells at zero expression, `q_low` quantile = 0, so the
   "low" boolean (`expr ≤ 0`) inflates and overlaps with "high" (`expr
   ≥ 0`, which includes zeros). Corners overlap and bootstrap measures
   deterministic intersections, not statistical noise.

2. **Original Schiebinger primary dataset is unreachable** (404 on
   `broadinstitute.github.io/wot/example_data`). Pilot pivoted to
   Bastidas-Ponce E15.5 pancreas via `scvelo.datasets.pancreas()`,
   which has discrete endocrine cell types (alpha/beta/delta/...) —
   the fundamentally different statistical regime where tertile cuts
   on hormone markers are meaningless.

These are content findings: the pilot surfaced design assumptions that
don't survive contact with this data type, BEFORE any pre-registration
was locked. Pilot did its job.

## What ran

Local laptop, `python3 /tmp/run_bio_pilot.py`. Compute time ~3 min
(after initial 50 MB pancreas download). Two iterations:

- **Iteration 1.** Marker-score outcome (`-mean(beta markers)`),
  median-corner statistic. Sanity 1 FAIL (|z|_p95 = 0.38), Sanity 2
  FAIL (|z| = 2.37 on Ins1/Ins2 — outcome contaminated by test genes
  themselves), Sanity 3 PASS (ρ = 0.736), Sanity 4 PASS by sample size.
- **Iteration 2.** Cluster-identity outcome (`1 - is_beta`),
  mean-corner statistic. Sanity 1 STILL FAIL (|z|_p95 = 0.31), Sanity
  2 PASS (|z| = 5.94, signal recovered after decoupling outcome),
  Sanity 3 PASS strongly (ρ = 0.951), Sanity 4 PASS.

The Iteration 1 → Iteration 2 jump on Sanity 2 (2.37 → 5.94)
demonstrates that **outcome decoupling matters**: marker-score-based
outcomes are *circular* when the marker genes themselves are tested.
This generalises to any future Schiebinger run — outcome must use
genes EXCLUDED from the test set.

## Why Sanity 1 still fails

Iteration 2 fixed the median→mean issue, but |z|_p95 = 0.31 ≪ 1.96
remains. Diagnosis: tertile cuts on sparse hormone HVG produce
degenerate corners. With many cells at zero expression, the q_low
quantile collapses to 0, and the two boolean masks (`a_low = expr ≤
0`, `a_high = expr ≥ 0`) are not disjoint — `a_high` includes all
zero-expression cells too. Corners HH and LL share most of their
zero-expression mass.

The bootstrap of differences-of-means between near-identical corners
gives tiny SE on tiny ε. Z = ε / SE has small numerator and small
denominator, dominated by floating-point and quantile-tie effects.

This is **not a bug in the bootstrap or in the design philosophy** —
it's the wrong statistical primitive for sparse categorical-marker
data.

## Three paths forward

### Path A — Find Schiebinger MEF→iPSC data elsewhere

The original primary dataset has continuous variation along the
reprogramming trajectory (early-MEF, mid-trajectory, late-iPSC, plus
off-trajectory cells). Tertile cuts are meaningful there. Sources:

- Broad Single Cell Portal SCP1010 (requires Broad/Google login)
- GEO GSE122662 (raw counts, ~30 GB — heavy preprocessing)
- Author's lab figshare? unknown
- Hugging Face datasets? unknown

This restores the originally-designed test. **Recommended if data
location is found.**

### Path B — Switch primary to a different continuous-trajectory dataset

`scanpy.datasets` includes:
- `paul15()` — myeloid progenitor differentiation (~2700 cells,
  continuous progenitor → MEP/GMP fates). Built-in, no download.
- `moignard15()` — early hematopoietic lineage (~3900 cells).

Both have continuous variation, no categorical marker dominance.
Methodology identical to Schiebinger plan. Top HVGs are TFs and
metabolic genes, not categorical markers — tertile cuts will be
meaningful.

### Path C — Adapt 2×2 design for categorical/sparse data

Replace tertile cuts with one of:

1. **Median split on nonzero subset.** `a_high = expr > median(expr[expr>0])`,
   `a_low = expr == 0`. Skip cells with `0 < expr < median`. Loses
   sample but keeps disjoint corners.
2. **Binary expressed/not.** `a_high = expr > 0`, `a_low = expr == 0`.
   Cleanest for hormone-like genes.
3. **Top-decile vs bottom-decile** instead of terciles. Requires denser
   genes; might still degenerate for ribosomal hormones.

Path C keeps the pivot to pancreas data but changes the statistical
primitive. Most permissive but least principled.

## Recommendation

**Path B** (switch primary to `paul15`) is the cleanest immediate move:
- Continuous trajectory data, original design assumption holds
- Built into scanpy, no download
- Same biological flavour as Schiebinger (lineage commitment)
- Allows pilot to PASS without re-engineering

**Path A** parallel: search for Schiebinger data through alternative
mirrors. If found, switch to Schiebinger as primary; paul15 becomes
secondary replication (parallel role to Bastidas-Ponce in the original
plan).

Path C reserved for the case where neither A nor B yields PASS — at
that point the categorical-data adaptation is its own methodological
contribution, but pre-reg becomes harder to justify.

## Files

- `/Users/teo/BioEpistasis_pilot/pilot_calibration_verdict.json`  
  Full numerical record from Iteration 2.
- `/Users/teo/BioEpistasis_pilot/top10_pairs.parquet`  
  45 top-10 HVG pair ε measurements (both outcomes).
- `/Users/teo/BioEpistasis_pilot/sanity3_outcome_consistency.png`  
  Marker vs distance scatter (ρ = 0.951; both outcomes agree).
- `/tmp/run_bio_pilot.py`  
  Pilot script, ready to retarget at a new dataset/preprocessing.

## Discipline note

This is **the second time in this project that a pilot caught a
methodological problem before a pre-reg locked**. The first was the
Phase 1 mean-vs-zero ablation finding in the sister Epistasis project,
which forced a switch to mean ablation before the production scan.
Both validate the discipline of *pilot before pre-reg* — both would
have produced wrong-direction results had we skipped the pilot.

The two findings should be added to `BioEpistasis/design_notes.md`
section 8 (discipline lessons) and to any future cross-project
methodology summary.
