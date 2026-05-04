# Pre-registration v1 — Functional epistasis on Norman 2019 K562 Perturb-seq

**Status: DRAFT.** To lock: review, freeze copy as
`bio_preregistration_v1_norman2019.LOCKED.md`, tag
`bio_prereg_v1_locked`, commit hash recorded both in this file's
header and in `bio_norman_verdict.json` produced by full scan.

**Locked-on-commit:** `<commit-hash-after-lock>`

---

## Conceptual frame

The sister project `mool32/epistasis-transformer-heads` established
that pairwise epistasis among top-30 attention heads in trained
transformers (Pythia 410M, OLMo-2 1B) is dominated by the
synthetic-lethal/redundancy regime: ε > 0 in loss space (joint
ablation hurts more than additive prediction) for 78% (Pythia) and
57% (OLMo) of significant pairs. This matches Costanzo 2010 *Science*
yeast genetic-interaction signatures (60-70% synthetic-sick/lethal).

This pre-registration tests whether **direct experimental
intervention** in human cells reproduces the same regime. Norman 2019
K562 combinatorial CRISPRi provides cells classified by guide-RNA
barcode rather than natural expression — soft correlation between
class indicator and outcome collapses to ~0 (verified in pilot iter 5,
C3 = 0.0003), breaking the soft-correlation × √n ceiling that bounds
observational designs (Constraint 4 in
`methodological_findings.md`).

**Why this matters.** The cross-substrate universality claim of the
ML programme currently relies on:
- Statistical match with **published** yeast literature (Costanzo
  2010 numbers as reference, no direct measurement)
- Calibration-scale match on Norman 2019 pilot (top-3 by |z|, n = 20
  pairs, exploratory)

A locked pre-registered full scan over all 131 testable Norman pairs,
with explicit pre-reg decision rules, converts this from "statistical
agreement with literature + small calibration sample" to "full-scale
direct intervention measurement under pre-registered locked criteria".

**Rules 1-6 binding** (inherited from ML pre-reg tradition v3+):
1. Direction predicted in advance: ε > 0 dominant. Negative direction
   = FAIL.
2. Single primary test, single decision.
3. Numeric thresholds pre-locked.
4. Null is legitimate.
5. No post-hoc reformulation.
6. Pre-reg commit hash baked into verdict JSON.

**Pilot iter 5 PASS conditions inherited:**
- C3 SOFT correlation (≈ 0.0003 in pilot) expected to hold at full
  scan; verified in `bio_norman_verdict.json`.
- C4 dominance (5.33 vs 0.11 = 50× margin in pilot) expected to hold;
  verified.
- S1 calibration (|z|_p95 = 1.79 in pilot) expected to hold under
  permutation null over the full pair set; verified.
- S2 sensitivity (curated pair |z| > 3) expected to be exceeded by
  the strongest pairs in the full scan; reported descriptively.

**Constraint 5 applied.** Per `methodological_findings.md`, binary
cluster outcome is saturating in Perturb-seq pair tests and
disqualified as a sanity-gate denominator. **Single primary outcome
design adopted.**

---

## 1. Dataset + cells

- **Source:** Norman et al. 2019, *Science* 365:786-793, GEO
  GSE133344. Mirror used: scPerturb harmonized release on Zenodo
  (`NormanWeissman2019_filtered.h5ad`, 698 MB, MD5
  `c870e6967d91c017d9da827bab183cd6`, DOI `10.5281/zenodo.10044268`).
- **Cell line:** K562 (chronic myeloid leukemia). One cell line, no
  cross-line comparison in this pre-reg.
- **Perturbation:** CRISPRi (dCas9-KRAB knockdown — gradient, not
  full knockout, structurally analogous to ML mean ablation).
- **Cells:** 111,445 total after scPerturb filtering.
  - Control (NT guide): 11,855
  - Single perturbations: 99,358
  - Pair perturbations: distributed across ~150 pre-curated pairs

## 2. Test set

**131 testable pairs.** A pair (A, B) is testable if all four classes
have ≥ 30 cells:
- WT (control, NT guide)
- A only (guide for gene A, NegCtrl0 for B)
- B only (guide for gene B, NegCtrl0 for A)
- AB pair (guides for both A and B)

Identified in pilot iter 5; locked in this pre-reg as the full
analysis set.

**No post-hoc filtering of pairs.** Every testable pair as defined
above enters the analysis. Excluded pairs (those with < 30 cells in
any class) are reported in verdict JSON but not analyzed.

## 3. Methodology

### 3.1 Outcome (single primary, locked)

**PCA distance to control-cluster centroid in 50-PC space.**

Procedure:
1. Standard scanpy preprocessing: filter genes (min_cells = 20),
   `normalize_total(target_sum = 1e4)`, `log1p`.
2. Highly variable genes: `cell_ranger` flavor, top 2000.
3. PCA: 50 components, seed 42.
4. Neighbors: 15 nearest, seed 42.
5. Leiden clustering: resolution 1.0, seed 42. Identify control-
   dominant cluster as the one in which control cells form the modal
   class.
6. Compute control centroid in 50-PC space as mean of cells in
   control-dominant cluster.
7. For every cell: outcome = Euclidean distance to control centroid.

**Why distance, not pseudotime or fate:** distance is monotone with
cell-state perturbation magnitude, continuous (not saturating per
Constraint 5), and uses the full HVG-derived state representation
rather than a single trajectory. Pseudotime requires a directional
choice (which fate?); fate probability requires fate-class
specification. Distance is neutral.

**Rejected outcomes** (per Constraint 5): cluster identity (binary,
saturates in pair tests). May appear in `bio_norman_verdict.json`
for descriptive reference but does not gate verdict.

### 3.2 Epistasis statistic (locked, loss-space convention)

For each pair (A, B):
```
Δ_A  = mean(outcome | A only)  − mean(outcome | control)
Δ_B  = mean(outcome | B only)  − mean(outcome | control)
Δ_AB = mean(outcome | A+B)     − mean(outcome | control)
ε    = Δ_AB − Δ_A − Δ_B
```

**Sign convention:** ε > 0 = synthetic-lethal/redundancy direction
(joint perturbation more deleterious than additive prediction).
ε < 0 = suppression/buffering direction.

### 3.3 Bootstrap SE

Resample cells **within each class** with replacement, 1000 bootstrap
iterations, seed = 42. SE = std of bootstrap distribution. z = ε / SE.
Significance gate: |z| > 3.

### 3.4 Statistical regularization (none)

No post-hoc effect-size cuts, no transformation of outcome before
comparison. Raw means in raw outcome units.

## 4. Primary pre-registered test

**Statistic.** Fraction of significant pairs (|z| > 3) with positive
ε:
```
frac_synth = #{pairs : |z| > 3 AND ε > 0} / #{pairs : |z| > 3}
```

**Predicted direction.** **frac_synth > 0.55** (synthetic-lethal /
redundancy regime dominant), matching Costanzo 2010 yeast (60-70%)
and Pythia/OLMo Tier 1 (78%, 57%).

**Null permutation test.** Pool all 131 pairs' (ε, SE, z); for each of
10,000 permutations, randomly assign sign to each ε (Bernoulli p=0.5),
recompute frac_synth_null. Two-sided p = fraction of permutations
with frac_synth_null ≥ frac_synth_obs. Locked permutation seed:
20260504.

**Four-tier decision (locked).**

| Tier | frac_synth | permutation p | Interpretation |
|---|---|---|---|
| **PASS** | > 0.55 | < 0.01 | Synthetic-lethal regime dominant in human Perturb-seq, matching ML and yeast |
| **PARTIAL** | 0.50 < frac_synth ≤ 0.55 | < 0.05 | Slight bias toward synthetic-lethal but not regime-dominant |
| **WEAK** | 0.45 ≤ frac_synth ≤ 0.50 | < 0.10 | Symmetric distribution |
| **FAIL_REVERSED** | < 0.45 | n/a | Suppression direction dominant — biology DIVERGES from ML and yeast |

The pre-registered prediction is **PASS**. Any other outcome is
content:
- PARTIAL: human Perturb-seq shows attenuated bias compared to yeast/ML
- WEAK: regime is symmetric in human cells
- FAIL_REVERSED: human cells show suppression-dominant regime,
  contradicting both yeast and ML — this would be a major content
  finding requiring re-interpretation of cross-substrate universality
  claim

## 5. Methodology gates (mandatory pre-flight)

Before reporting primary verdict, all five gates must hold:

| Gate | Threshold | Pilot iter 5 value | Notes |
|------|-----------|-------------------|-------|
| **G1 calibration** | \|z\|_p95 on permutation null ∈ [1.5, 2.5] | 1.79 | Verified pre-data on full pair set |
| **G2 sensitivity** | At least one curated pathway pair shows \|z\| > 3 | CBL/CNN1: 11.27 | Pre-specified curated list |
| **G3 sample size** | Min class size ≥ 30 across all 131 pairs | 59 | By construction (filter applied) |
| **G4 SOFT correlation** | \|corr(class_indicator, outcome)\| < 0.05 | 0.0003 | Independence of guide barcode from PCA distance by construction |
| **G5 dominance** | median \|z\|_sig > max \|corr\| × √n_eval | 5.33 vs 0.11 | Pilot ratio 50× |

If any gate fails on full scan, verdict is HOLD pending investigation;
primary test not reported until gate resolved.

## 6. Mandatory secondary tests

All computed and reported in verdict; none gate the primary verdict.

### 6.1 Magnitude — median |ε|_significant

Median |ε| of significant pairs (|z| > 3). Reported alongside
direction. No threshold; descriptive only.

### 6.2 Distribution shape

AIC fit of |ε|/SE distribution to Gaussian / Laplace / Student-t.
Predicted: Student-t (heavy-tailed, matching Pythia/OLMo Tier 1
shape).

### 6.3 Pathway coherence

For pairs both in same KEGG / GO BP pathway: median |ε|. Compare to
cross-pathway pairs. Report Mann-Whitney U one-sided p (same >
cross). Tests the operon analog observed in ML (same-layer 4.5–7×
cross-layer in Pythia/OLMo).

### 6.4 Cross-substrate magnitude comparison

For all four substrates with available frac_synth measurements:
- Pythia 410M Tier 1: 0.78
- OLMo 1B Phase 4: 0.57
- Yeast Costanzo 2010: ~0.60-0.70 (literature range)
- **Norman 2019 K562:** verdict value

Report 4-substrate comparison table. No statistical test (literature
vs measurement comparison is descriptive).

### 6.5 Top-K identification (descriptive only)

Top-30 pairs by |ε|/SE in Norman scan. Report (gene_a, gene_b) list
with |z|. No claims about specific gene identities; this serves
follow-up work (mechanism, pathway annotation).

## 7. Compute

CPU-only, Colab CPU runtime sufficient.

- Per-pair compute: ~0.5 s (4-class means + 1000 bootstrap)
- 131 pairs × 0.5 s = 66 s for sequential scan
- Plus permutation null: 10,000 × ~0.05 s = 500 s
- Plus secondary analyses: ~5 minutes
- Plus pre-flight gates: ~2 minutes

**Total: ~10-15 minutes wall-clock** on Colab CPU. Compute is not the
bottleneck; design and pre-reg discipline are.

## 8. Artifacts

Produced by full scan and committed to `data/analysis/norman/`:
- `norman_pairs.parquet` — 131 rows, schema: gene_a, gene_b,
  class_n_ctrl, class_n_a, class_n_b, class_n_ab, delta_a, delta_b,
  delta_ab, epsilon, epsilon_se, z_score, same_pathway (bool).
- `norman_pairs_perbatch.npz` — per-bootstrap-iteration ε arrays
  for re-analysis without re-scanning.
- `bio_norman_verdict.json` — pre-reg commit hash, primary
  frac_synth + permutation p, all five gates' values, all
  secondaries, decision tier.
- `bio_norman_headline.png` — primary figure: ε distribution
  (significant pairs), sign breakdown, cross-substrate comparison
  table.

## 9. What is NOT pre-registered

- Cell-cycle / batch / UMI-depth confound corrections. Norman dataset
  is from one cell line, one experiment, so batch is constant.
  Cell-cycle correlation reported in secondary 6.5 if relevant but not
  used in primary.
- Mechanistic claims about specific pairs. Top-K table is descriptive.
- Cross-cell-line replication. Norman is K562-only here. Other cell
  lines (e.g., RPE1, A375) require separate pre-reg.
- Joint-mean ablation alternative (versus independent means). Locked
  to independent means in pilot.

If primary **FAILS_REVERSED**, the cross-substrate universality claim
of the ML programme weakens substantially:
- ML and yeast both show synthetic-lethal-dominant regime
- Human Perturb-seq showing suppression-dominant regime would
  indicate the regime is NOT substrate-independent in the strict
  sense — depends on system type
- This becomes the headline finding (negative result for universality;
  positive content for "substrate matters")

If primary **PARTIAL** or **WEAK**, the regime is preserved
qualitatively but attenuated quantitatively. Reportable as "regime
holds in direction but not in strength of bias".

If primary **PASS** (predicted), the four-substrate universality
claim is established at full-scale measurement: trained transformers,
evolved yeast, intervened human cells all show synthetic-lethal /
redundancy regime dominance.

---

*Draft 2026-05-04. To lock: review, copy to `…LOCKED.md`, tag
`bio_prereg_v1_locked`. One-shot. No rescue.*
