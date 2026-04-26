# BioEpistasis — Design notes

Final design after the 2026-04-26 design discussion. Captures what is
locked, what is open, and the sign-convention safeguard placed at the
top of the document explicitly to prevent the labeling error that
propagated through the sister Epistasis project's pre-regs v1/v2/v3.

---

## 1. Sign convention (READ FIRST. Locked.)

ε is an empirical quantity in **loss space** with **additive null**:

    ε_loss = Δ_AB − Δ_A − Δ_B
    where Δ_X = L_X − L_baseline ≥ 0 for deleterious perturbations

The biological literature (Costanzo et al. 2010 "The Genetic Landscape
of a Cell" and successors) uses **fitness space** with **multiplicative
null**:

    ε_fitness = f_AB − f_A · f_B
    where f_X = fitness of mutant X relative to WT, f_X ≤ 1

The two conventions have **opposite signs** for the same biological
phenomenon:

| Phenomenon | Single-mutant phenotype | Joint phenotype | ε_fitness | ε_loss |
|------------|-------------------------|-----------------|-----------|--------|
| **Synthetic-lethal / redundancy** (Costanzo) | mild | catastrophic | **< 0** | **> 0** |
| **Suppression / buffering / true compensation** | deleterious | partially rescued | **> 0** | **< 0** |

Concrete worked example (two functionally redundant attention heads):
- Δ_A = 0.01 (B compensates when A alone is ablated)
- Δ_B = 0.01 (A compensates when B alone)
- Δ_AB = 0.10 (no one compensates when both are ablated)
- ε_loss = 0.10 − 0.01 − 0.01 = **+0.08 → ε > 0**
- This is the canonical *synthetic-lethal* phenotype in Costanzo's
  framework. The English word "compensatory" in some ML papers refers
  to this same redundancy phenotype, **but that creates ambiguity**
  with Costanzo's "compensation" / "suppression" which means the
  *opposite* sign.

**Rule for this project.** Always use **Costanzo terminology in the
loss-space convention**:
- ε_loss > 0  →  "synthetic-lethal / redundancy phenotype"
- ε_loss < 0  →  "suppression / buffering"

Pythia Tier 1 finding (78 % significant pairs ε_loss > 0) IS biology
parallel — most epistasis is synthetic-lethal-direction, as in
Costanzo. The earlier "reversed_compensatory_dominant" label in the
sister project's `tier1_verdict.json` is a *terminological* error;
the underlying numbers are correct.

---

## 2. Conceptual framing

**Hypothesis (parallel direction).** Both ML training and biological
differentiation drive top functional units toward synthetic-lethal /
redundancy dominance. Reasoning:

- *ML*: training accumulates redundant heads in the residual stream;
  top-K functional heads form tightly coupled groups whose joint
  knockout is catastrophic relative to single knockouts.
- *Biology*: differentiation establishes paralog systems and
  master-regulator hubs; top-variance genes in differentiated states
  acquire critical functional dependencies that manifest as
  synthetic-lethal / redundancy phenotype.

**Direction predicted.** `fraction(ε_loss < 0)` decreases with
differentiation depth; equivalently, `fraction(ε_loss > 0)` (the
synthetic-lethal/redundancy fraction) increases.

This is the same numerical prediction as the original "Mapping B" but
**reframed**: the parallel arises through *both* substrates accumulating
synthetic-lethal regime via specialization, not through one being
compensatory and the other not.

---

## 3. Methodology (locked for pilot, candidate for pre-reg)

### 3.1 What 2×2 expression stratification actually measures

We do **not** measure causal ablation. Observational scRNA-seq data
permits only:

    stratify cells by natural expression of A and B  →  compare outcomes

This is the **conditional expectation** analog. It mirrors **mean
ablation** in ML (which replaces a head's output with its dataset-mean,
i.e. reads out the "what if this signal were at its expected value
across the population" question). Both compare cell/forward-pass
outcomes between "this signal is at its low natural level" vs "high
natural level", not "this signal is intervened to zero".

The pre-reg must state this explicitly to head off the predictable
"you are measuring conditional expectation, not ablation" objection.
The reframe: this project tests whether the *same conditional-
expectation analog of epistasis* shows the regime predicted by ML.

### 3.2 Outcome (locked)

**Primary:** WOT-derived trajectory probability,
`log P(reach iPSC by t+1 | current state)`. Continuous, defined
naturally for cells off the main trajectory (low probability rather
than meaningless distance), and standard in the Schiebinger pipeline.

**Secondary (consistency):** Euclidean distance to the iPSC cluster
centroid in PCA space (top 50 PCs). Run as cross-check of outcome
choice. Disagreement between primary and secondary means outcome is
fragile and outcome must be reconsidered before pre-reg.

### 3.3 Gene set (locked)

**Top-30 highly variable genes** (HVG) by `scanpy.pp.highly_variable_
genes(flavor='seurat_v3')` at the **terminal** timepoint (analog of
Phase 2B final-checkpoint top-K). Same 30 genes tracked across
timepoints in the trajectory test.

For the pilot: top-10 HVG (smaller set for calibration speed).

### 3.4 ε computation per pair

- Pair (A, B) chosen
- All cells in the timepoint binned by tertile of expression of A
  *and* of B → **9 bins**, of which we use the 4 corners:
  - high_A high_B  (analog of "WT")
  - low_A  high_B  (analog of "A perturbed")
  - high_A low_B   (analog of "B perturbed")
  - low_A  low_B   (analog of "both perturbed")
- Outcome φ in each corner = median of the outcome variable
- Δ_A  = φ(low_A high_B)  − φ(high_A high_B)
- Δ_B  = φ(high_A low_B)  − φ(high_A high_B)
- Δ_AB = φ(low_A low_B)   − φ(high_A high_B)
- ε    = Δ_AB − Δ_A − Δ_B   (loss-space convention; sign per §1)

**Bootstrap SE.** Resample cells *within each corner* with replacement
1000 times; recompute ε; SE = std of bootstrap distribution.
**z-score** = ε / SE. Significance gate (initial pilot): |z| > 3.

**Sample size guard.** Each corner must contain ≥ 30 cells; if any
corner has < 30, the pair is skipped and recorded as
"underpowered_in_quadrant".

### 3.5 Architectural baseline (mandatory secondary)

50 random gene pairs at each timepoint, sampled from the full HVG pool
with seed=42. Median |ε|_T2 → contemporaneous architectural baseline.
Direct analog of the sister-project Phase 2A.

This is an **autonomous finding** — even if the trajectory primary
fails, observing nonzero architectural epistasis baseline in biology
parallels the residual-stream baseline in ML and is reportable on
its own.

### 3.6 Confound controls

- **Cell cycle.** Repeat all primary analyses after regressing out
  cell-cycle scores (S, G2M phases via `scanpy.tl.score_genes_cell_
  cycle`). Substantive change → cell cycle confound is serious.
- **Batch.** If the timepoint draws from multiple batches, repeat
  per-batch and report consistency.
- **Total UMI.** Match cells by total-UMI quantile across the four
  corners to control for sequencing depth.

### 3.7 Multiple testing

- Per-pair: BH FDR control, q < 0.05.
- Trajectory primary: a single Spearman ρ across timepoints applied
  to fraction(ε_loss < 0)(t). No multiple testing on the primary.

---

## 4. Predictions (locked)

### 4.1 Primary

Spearman ρ between timepoint index and fraction(ε_loss < 0):

    ρ < 0  (fraction of suppression-direction pairs DECREASES as cells
            differentiate; equivalently, synthetic-lethal/redundancy
            fraction increases)

**Decision rule.**
- **PASS:** ρ < 0  AND  permutation p < 0.05  AND  fraction drops
  ≥ 0.10 from earliest to latest timepoint.
- **FAIL_OPPOSITE:** ρ > 0 with p < 0.05 (biology develops
  *more* suppression-direction pairs over differentiation, opposite
  to ML).
- **FAIL_NULL:** |ρ| < threshold or p ≥ 0.05 → no detectable trajectory.

All three outcomes are content. PASS supports the universality claim
across substrates. FAIL_OPPOSITE means biology and trained ML diverge,
which is itself worth reporting.

### 4.2 Secondaries (mandatory)

- **Same-pathway enrichment.** Pairs of genes sharing a KEGG / GO BP
  term should show stronger |ε| than cross-pathway pairs (operon
  analog of the sister project's same-layer enrichment finding).
- **Architectural baseline at each timepoint.** Predict ε_T2 ≠ 0
  (gene regulatory network coupling, ribosome competition, shared
  TFs). Magnitude and stability across timepoints reported.
- **Magnitude scaling.** Predict |ε| roughly proportional to
  |Δ_A · Δ_B| (product of single effects), as observed in the sister
  project at top-K scale.
- **Distribution shape.** AIC fit |ε|/SE to Gaussian / Laplace /
  Student-t. Predict Student-t (heavy-tailed), parallel to Pythia
  Tier 1 result.

---

## 5. Pilot scope (locked, runs first)

### 5.1 Dataset and timepoint

- **Schiebinger 2019** reprogramming MEF→iPSC scRNA-seq (~250k cells,
  18 timepoints, days 0 → 18).
- Pilot uses **day 8 only** (mid-trajectory). One cluster /
  pre-existing batch label if multi-batch.
- Bastidas-Ponce pancreatic dataset is *not* in pilot — reserved for
  potential replication after pre-reg PASS.

### 5.2 Gene set

- HVG via `scanpy.pp.highly_variable_genes(flavor='seurat_v3')` on
  the day-8 subset.
- **Top-10 HVG** → 45 pairs (small for pilot speed, scales to top-30
  for full run).

### 5.3 Outcomes computed (both, for cross-check)

- WOT log P(reach iPSC by t+1).
- Euclidean distance to iPSC centroid in 50-PC PCA space.

### 5.4 Sanity checks (calibration verdict)

The pilot **must** pass all four to permit pre-reg drafting. Any single
FAIL → halt, reconsider design before locking.

1. **Bootstrap calibration on random pairs.** 50 random gene pairs.
   The empirical distribution of |z| should be approximately N(0, 1).
   - PASS: empirical 95th percentile ∈ [1.5, 2.5] (theoretical 1.96).
   - FAIL: outside this range → bootstrap is mis-calibrated; cell
     resampling assumption (i.i.d. cells within corner) does not
     hold and a different SE estimator is needed.

2. **Sensitivity on a known same-pathway pair.** Pick a
   pre-registered curated pair: two ribosomal proteins
   (`Rps3`, `Rps5`) or two glycolysis enzymes (`Eno1`, `Pgk1`).
   - PASS: |z| > 3 on the pilot pair.
   - FAIL: data resolution insufficient at this depth; consider
     stricter expression cutoffs (extreme deciles instead of
     terciles) or different outcome measure.

3. **WOT vs distance consistency on top-10 pairs.** Compute ε for the
   45 top-10 HVG pairs under both outcomes. Pearson ρ between the two
   ε vectors.
   - PASS: ρ > 0.7 (outcomes agree on epistasis directionality).
   - FAIL: outcome choice is fragile; primary must be redefined.

4. **Sample size per quadrant.** For each top-10 pair, all four
   corners (low_A low_B, etc.) contain ≥ 30 cells.
   - PASS: ≥ 30 in all corners across all 45 pairs.
   - FAIL: tertile split with 10 genes leaves rare quadrants
     under-sampled at this cell count; switch to median split or
     widen tertile thresholds.

### 5.5 What is NOT in pilot

- Multiple timepoints. Day 8 only.
- Trajectory primary. Computed on full pre-reg run, not pilot.
- Bastidas-Ponce replication.
- Cell-cycle / batch / UMI confound controls. Pilot only.
- Pre-registration drafting.

### 5.6 Compute

Schiebinger day-8 has ~14 k cells. HVG + 45 pairs × 1000 bootstrap
runs is comfortably handled by Colab CPU runtime in <1 h. No GPU
needed.

---

## 6. Decision tree after pilot

| Pilot result | Action |
|--------------|--------|
| All 4 sanity PASS | Draft pre-reg v1, lock, run full Schiebinger trajectory (all timepoints, top-30 HVG). |
| Sanity #1 FAIL (calibration) | Try cluster-bootstrap (cells grouped by batch/cluster). If still FAIL, switch to Mann-Whitney based ε null instead of bootstrap-z. |
| Sanity #2 FAIL (sensitivity) | Switch tertile → decile (10 % vs 90 %) split. Re-run pilot. |
| Sanity #3 FAIL (outcome consistency) | Investigate which outcome is closer to expected behaviour on curated pair from #2; pick that as primary, drop the other. |
| Sanity #4 FAIL (sample size) | Switch tertile → median split. Re-run pilot. |
| Multiple FAILs | Reconsider whether 2×2 stratification is the right design at all on this data scale. Revisit Variant 3 (manifold deviation) as primary. |

---

## 7. After pre-reg PASS, downstream

If Schiebinger trajectory primary PASSES:
- Replicate on Bastidas-Ponce pancreatic development (separate pre-reg).
- Cross-organism: zebrafish or fly developmental atlases.
- Cross-substrate paper: ML, yeast (Costanzo), mammalian development
  all show synthetic-lethal/redundancy regime in top functional units.
  Universality of epistasis structure.

---

## 8. Discipline lessons recorded for cross-project reuse

1. **Sign convention must be stated and verified before any
   interpretation.** This project's section 1 was placed at the top
   of the design document, immediately, because the sister Epistasis
   project's pre-regs propagated a sign-convention error through three
   locked documents and eight sessions before being caught.
   Numerical decision rules survived; narrative interpretation did
   not.

2. **Stratification ≠ ablation.** Conditional-expectation analog of
   ablation, properly named, is more defensible than calling it
   ablation and then defending the gap.

3. **Pilot before pre-reg.** Locking decision rules on a design that
   has not been calibrated risks locking around assumptions that the
   data may not support (bootstrap calibration, sample size, outcome
   sensitivity).

---

*Initial version 2026-04-26. To update: edit, append a "Revision N"
section preserving the prior content, and note in the README. Do not
silently rewrite locked numerical decisions once pilot calibration
PASSes.*
