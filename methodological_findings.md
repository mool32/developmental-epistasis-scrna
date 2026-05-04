# BioEpistasis — Methodological findings (locked design rules)

This document catalogues design constraints discovered through pilot
iterations BEFORE any pre-registration is locked. Each constraint must
be verified by an explicit sanity check in any future pilot or full
pre-reg run.

Document is append-only. New constraints emerge from new failure modes;
existing constraints are not retracted.

---

## Constraint 1 — Outcome must be functionally decoupled from test genes (HARD)

### Refinement after Iteration 3 (2026-04-26): hard vs soft circularity

Iteration 3 on paul15 surfaced a finer distinction: there are **two
kinds** of outcome–gene coupling, with different consequences.

**HARD circularity (blocking).** Outcome computed *directly* from a
test gene's expression. Iteration 1's `outcome = -mean(beta markers)`
was hard-circular when testing `(Ins1, Ins2)` — those genes WERE the
outcome. 2×2 stratification then conflates "expression of test gene"
with "outcome value", producing artefactual ε.
**Verification:** the outcome formula must not mention any test gene
by name. PASS = formula and test-set are disjoint.

**SOFT correlation (transparent, not blocking).** Outcome is a
trajectory/cell-state metric (pseudotime, KNN-cluster identity,
lineage probability) that happens to correlate with lineage-marker
test genes because *both* reflect the same hidden differentiation
state — neither is computed from the other directly. Iteration 3
showed `|corr(outcome_pseudotime, Car1)| = 0.48` arose this way:
pseudotime came from the full-HVG neighborhood graph; Car1 was simply
a marker of the lineage that pseudotime tracks.

Soft correlation biases the null toward zero ε (cells with low Car1
are systematically "less differentiated" in outcome, which is what
the stratification measures), but the bias is **bounded by the
correlation magnitude**, not catastrophic. Real ε signals on lineage
TFs in Iteration 3 reached |z| = 20+ — orders of magnitude above any
plausible spurious contribution from |corr| ≤ 0.5.

**Updated rule.**
- HARD circularity: outcome formula mentions any test gene → FAIL.
  This is non-negotiable; Iteration 1's ε = −2.4 on Ins1/Ins2 vs ε =
  +5.9 in Iteration 2 (cluster outcome, no formula overlap) is the
  evidence.
- SOFT correlation: report `|corr(outcome, gene)|` for every test
  gene in the verdict JSON. Threshold > 0.5 warrants an explicit
  caveat in interpretation. Threshold > 0.8 indicates the outcome is
  effectively *redefined* by the test set and warrants picking a
  different outcome (e.g., cell-cycle phase, batch label) at the cost
  of biological interpretability.



**Source:** Pilot Iteration 1 on Bastidas-Ponce pancreas (2026-04-26),
Sanity 2 FAILed (|z| = 2.37 on Ins1/Ins2 paralog pair) when
`outcome = -beta_marker_score` and beta_marker_score included Ins1,
Ins2. Iteration 2 swapped to cluster-identity outcome (decoupled from
any single gene's expression) and Sanity 2 jumped to |z| = 5.94 — same
biology, signal recovered after decoupling.

**Failure mode.** When the outcome is a function of the test genes'
expression, the 2×2 stratification creates artificial correlations
between corner labels and outcome values that have nothing to do with
biological epistasis.

**Rule.** Outcomes must be computed from a feature set DISJOINT from
the test gene set. Acceptable outcomes:

- Pseudotime (`sc.tl.dpt`) — derived from cell graph topology, not
  any single gene
- KNN-cluster identity (binary, "is in target cluster") — discrete
- Lineage probability from external optimal-transport computation
  (CellRank, WOT) — uses cell trajectories, not direct gene expression
- TF panel score where the TFs are explicitly NOT in the HVG/test set

**Mandatory sanity check (any pilot or pre-reg).** For every test
gene g in the chosen test set:

    |Pearson correlation(outcome, expression_g)| < 0.3

If any test gene exceeds this threshold, either drop the gene from the
test set, or pick a different outcome.

---

## Constraint 2 — 2×2 tertile design requires continuous expression variation

**Source:** Pilot Iteration 2 on Bastidas-Ponce pancreas (2026-04-26),
Sanity 1 FAIL (|z|_p95 = 0.31, target 1.96) caused by degenerate
tertile cuts on sparse categorical hormone HVG (Sst, Ghrl, Ins1,
Ins2, Gcg, ...). With > 33 % cells at zero expression, q_low quantile
collapses to 0, and the boolean masks `a_low (expr ≤ 0)` and
`a_high (expr ≥ 0)` are not disjoint — `a_high` includes all zero
cells. Corners HH, LL, HL, LH overlap, bootstrap measures
deterministic mass intersections rather than statistical noise.

**Failure mode.** The 2×2 design implicitly assumes test genes have
continuous expression variation across cells. Categorical markers
(genes that are essentially "on" in one cell type and "off"
elsewhere) violate this assumption — the tertile boundary collapses
into a binary expressed/not split.

**Rule.** Test gene set must have continuous variation. Specifically:

- `q_low(expression_g) > 0` for every test gene g, OR
- Fraction of nonzero-expressing cells > 50 % for every test gene g

This excludes hormone-style markers (Ins, Gcg, Sst) and any other
gene that is silent in most cells. Suitable test gene sets emerge
naturally from continuous-trajectory datasets (Schiebinger, paul15,
moignard15) where top HVGs are TFs and metabolic genes — broadly
expressed with cell-state-dependent variation.

**Mandatory sanity check.** Compute q_low for each test gene
in the chosen set; assert > 0. If any gene fails, either remove it
from the test set, or switch to an alternative stratification design
(binary expressed/not, or top-decile vs bottom-decile on positive
subset only).

**Datasets that satisfy Constraint 2** (verified or expected):

- `scanpy.datasets.paul15()` — myeloid lineage TFs, expected continuous
- `scanpy.datasets.moignard15()` — early hematopoietic, expected continuous
- Schiebinger 2019 MEF→iPSC mid-trajectory (when data accessible) —
  continuous reprogramming intermediates

**Datasets that violate Constraint 2** (verified):

- `scvelo.datasets.pancreas()` Bastidas-Ponce E15.5 — hormone-marker
  dominated, top-10 HVG categorical. Adaptable via Path C
  (binary or top-decile design) but not under standard tertile.

---

## Constraint 3 — HARD vs SOFT circularity (refinement of Constraint 1)

**Source:** Pilot Iteration 3 on Paul15 (2026-04-26). Iter 3 ran with
pseudotime as outcome (graph-derived, not direct gene formula) but
showed |corr(outcome_pseudotime, top-10 HVG)| up to 0.48. This is
correlation through shared latent state (lineage), not formula
inclusion.

**HARD circularity (blocking).** Outcome formula directly references
test gene expression. Iter 1 case. Always disqualifies.

**SOFT correlation (transparent, conditional).** Outcome and test gene
share information about latent cell state but neither computed from
the other. Iter 3 case. Bias is bounded by correlation magnitude.

**Resolution.** Constraint 1 became HARD-only blocking. SOFT correlation
is reported in verdict JSON, threshold > 0.5 warrants caveat,
threshold > 0.8 forces outcome change.

---

## Constraint 4 — Soft-correlation × √n dominance ceiling

**Source:** Pilot Iteration 4 on Paul15 + synthetic null (2026-04-26).

**Failure mode.** On observational scRNA-seq with cell-state outcomes,
SOFT correlation between test-gene expression and outcome amplifies
linearly with √n. With max |corr| = 0.48 and n ≈ 2700, the
amplification term reaches |z| ≈ 25 — comparable to or exceeding the
strongest observed real signals (median |z|_significant = 8). Real
biology cannot be statistically distinguished from the soft-correlation
contribution.

**Quantitative criterion.** `median |z|_sig > max |corr| × √n_eval`
must hold to claim signal dominates spurious correlation.

**Resolution.** This is the CORE limitation of observational designs
with cell-state outcomes. **Direct experimental intervention
(Perturb-seq) breaks the ceiling** because guide barcode is
independent of expression by experimental construction → C3 SOFT
correlation collapses to ~0 → C4 dominance margin opens up.

This finding motivated the pivot from observational scRNA-seq
(Bastidas-Ponce, Paul15) to Perturb-seq (Norman 2019).

---

## Constraint 5 — Binary outcomes saturate in Perturb-seq pair tests

**Source:** Pilot Iteration 5 on Norman 2019 K562 Perturb-seq
(2026-04-28). Verdict at
`pilot/iter5_norman/iter5_verdict.json`.

**Failure mode.** With binary cluster-identity outcome (1 = "not in
control cluster", 0 = "in control cluster"), single-perturbation cells
saturate near outcome 1: a single guide is enough to push the cell out
of the control transcriptional state, so single-class-A and single-
class-B mean outcomes are ≈ 1. Pair-class AB cannot be pushed further
than ≈ 1. Therefore Δ_AB ≈ Δ_A ≈ Δ_B → ε(cluster) saturates
artificially in the suppression direction (ε < 0 mechanically).

This is structurally analogous to a saturating dose-response: real
synthetic-lethal biology may be present but the binary readout cannot
detect it because the readout has hit ceiling.

**Empirical signature.** S3 outcome-consistency gate (Pearson ρ
between ε(continuous) and ε(binary) > 0.7) FAILED on Norman pilot
calibration: |ρ| = 0.198. Cluster ε saturated negative while distance
ε grew positive across the same pairs.

**Worked example (Norman pilot).** CBL/CNN1 pair:
- ε(distance, continuous) = +1.34, z = +11.27 (synthetic-lethal direction)
- ε(cluster, binary) = -0.015, z = -0.56 (saturating, no signal)

Same biology; different outcomes; opposite measured directions because
binary outcome saturates.

**Rule.** For Perturb-seq combinatorial pair tests:
- **Primary outcome must be continuous** (PCA distance to reference
  centroid, pseudotime, fate probability, perturbation-phenotype
  magnitude).
- **Secondary outcome (for outcome-consistency cross-check) must
  also be continuous.** Binary cluster identity is permissible as
  *descriptive* report (e.g., "fraction of pair-class cells outside
  control cluster") but not as a sanity-gate denominator.
- If only binary outcome is available, the S3 outcome-consistency
  gate is suspended and the verdict is reported as single-outcome
  with explicit note.

**Verification check.** For any new outcome under consideration:
compute mean outcome value in single-class-A vs pair-class-AB cells
across a range of pairs. If mean(AB) and mean(A) are within 5% for
the majority of pairs, outcome is saturating and disqualified.

**Datasets / designs that satisfy Constraint 5** (verified or expected):
- Norman 2019 with PCA distance as primary: continuous, monotone
  with cell-state distance, PASSES.
- Replogle 2022 perturbation-phenotype magnitude metric: continuous
  by design, expected to PASS.
- Adamson 2016 UPR-pathway score: continuous, expected to PASS.

**Datasets / designs that violate Constraint 5** (verified):
- Norman 2019 with binary cluster identity: saturates in single
  perturbations, FAILS for pair tests.
- Any binary "responder/non-responder" outcome where single
  perturbation already saturates the threshold.

---

## Future use

Constraints 1-5 become explicit sanity gates in pre-reg drafts. Their
verification reports go into the verdict JSON of every pilot/pre-reg
run. Failure to verify any constraint blocks the run unless explicitly
addressed (e.g., HARD vs SOFT split).

The discipline rule: pilot first, gaps surface, they get documented
here as design rules, then pre-reg locks under those rules. **Five
iterations have surfaced five distinct constraints.** Each was real,
each would have produced wrong-direction or non-detectable results
had pre-reg locked without it.

---

*2026-04-26 / 2026-05-04 (C5 added). Append-only. Updated after every
pilot iteration that surfaces a new constraint.*
