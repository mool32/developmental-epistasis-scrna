# Observational epistasis on scRNA-seq: limits of 2×2 stratification

*Technical note. Synthesis of pilot iterations 1–4 (2026-04-26) on
Bastidas-Ponce E15.5 pancreas and Paul15 myeloid lineage data. Internal
documentation; ready for standalone publication or as Part 1 of a
follow-up study using Perturb-seq.*

---

## Abstract

We sought to test whether the synthetic-lethal/redundancy regime
observed in transformer attention heads (Pythia 410M Tier 1: 78 % of
significant top-30 pairs ε > 0 in loss space) has a developmental
analog in cell differentiation. The proposed instrument was 2×2
expression-tertile stratification on observational scRNA-seq data —
the conditional-expectation analog of mean ablation.

Four pilot iterations on two scanpy-resident datasets surfaced four
methodological constraints. Under the cumulative constraint set, the
design fails to lock pre-registration on currently-tested datasets.
Specifically, the dominance criterion `median |z| > max |corr| × √n`
fails because soft outcome–gene correlation, an unavoidable property
of cell-state outcomes computed from gene expression, can numerically
account for the magnitude of observed signals.

This is a methodology contribution, not a closing of the question.
The constraints define when 2×2 stratification IS valid and motivate
a pivot to Perturb-seq data, where genetic perturbations are
experimentally introduced rather than read off natural variation.

---

## 1. Question

Pythia 410M Tier 1 (sister project): 78 % of top-30 attention-head
pair epistasis values ε_loss > 0, i.e. dominantly synthetic-lethal /
redundancy phenotype. This matches Costanzo 2010's biological prior
(yeast double-mutants are most often synthetic-sick/lethal, after
sign-convention translation between fitness and loss spaces).

Question: do mammalian developmental cells, in their natural
expression heterogeneity, show the same regime among their
top-variance genes? If yes, the synthetic-lethal regime is a
universal property of differentiated functional units across
substrates. If no, ML and biology diverge — interesting either way.

## 2. Method (proposed)

For a gene pair (A, B) and a continuous outcome φ:

- Stratify cells by expression of A (low/high tertiles) × B
  (low/high tertiles) → four corners HH / HL / LH / LL.
- Compute φ-mean per corner.
- Δ_A = φ(LH) − φ(HH), Δ_B = φ(HL) − φ(HH),
  Δ_AB = φ(LL) − φ(HH).
- ε = Δ_AB − Δ_A − Δ_B (loss-space additive null).
- Bootstrap SE via resampling cells within each corner.
- Significance: |z| = |ε|/SE > 3.

Sign convention is fixed: ε > 0 = synthetic-lethal / redundancy
(Costanzo's *negative epistasis* in fitness space); ε < 0 =
suppression / buffering (Costanzo's *positive epistasis*).

This is the conditional-expectation analog of mean ablation in ML:
both compare cell/forward-pass outcomes between low and high values
of a signal, not against a counterfactual zeroed intervention.

## 3. Pilot iterations

### Iteration 1 — Bastidas-Ponce pancreas, marker-score outcome, median statistic

Outcome: `−mean(beta-cell markers Ins1, Ins2, Pdx1, Mafa, ...)`.

- Sanity 1 (calibration on random pairs): FAIL, |z|_p95 = 0.38
  (target [1.5, 2.5]).
- Sanity 2 (sensitivity on curated pair Ins1, Ins2): FAIL, |z| = 2.37
  (target > 3).
- Sanity 3 (outcome consistency): PASS, ρ = 0.74.
- Sanity 4 (sample size): PASS by sample size, but *illusory* — most
  of the cells are zero-expressing for hormone HVGs, inflating corner
  counts deceptively.

**Finding 1 — HARD CIRCULARITY.** Outcome (`−beta-marker score`)
explicitly *includes* Ins1, Ins2 in its formula. When testing the
(Ins1, Ins2) pair, the 2×2 corners are partitioned by the same gene
expressions that drive the outcome. The signal is autocorrelation,
not biology.

### Iteration 2 — Bastidas-Ponce pancreas, cluster-identity outcome, mean statistic

Two changes from Iter 1:
- Outcome: `(cluster ≠ Beta)` binary — formally disjoint from any
  single test gene's expression.
- Statistic: mean per corner (not median) — bootstrap on median is
  unstable on sparse scRNA data because resampled medians flip
  discretely between mass points.

- Sanity 2 jumped from |z| = 2.37 → |z| = 5.94. Outcome decoupling
  RECOVERED the signal. Confirms Finding 1.
- Sanity 1 still FAILS, |z|_p95 = 0.31.

**Finding 2 — TERTILE DEGENERACY ON SPARSE DATA.** Top-10 HVG in E15.5
pancreas are categorical hormone markers (Sst, Ghrl, Ins1, Ins2, Gcg,
Npy, ...). Each is expressed in one endocrine cell type and silent
elsewhere — > 33 % of cells at expression zero. Tertile threshold
q_low collapses to 0, then `a_low := expr ≤ 0` and
`a_high := expr ≥ 0` overlap (zeros in both), and corners HH/LL/HL/LH
have largely shared cells. Bootstrap measures mass-overlap
deterministics, not statistical noise.

The 2×2 design implicitly assumes test-gene expression varies
*continuously*. On categorical lineage-marker datasets, the design is
not applicable in its tertile form. Workarounds (binary expressed/not,
top-decile vs bottom-decile on positive subset, etc.) exist but lose
the standard tertile interpretation.

### Iteration 3 — Paul15 myeloid lineage, pseudotime outcome, mean statistic

Switched primary to Paul15 (continuous trajectory: progenitor → MEP /
GMP fates), top HVGs are TFs and lineage markers with continuous
variation. Outcome: diffusion pseudotime (`sc.tl.dpt`), derived from
the full neighborhood graph, not from any single gene.

- Sanity 1 (calibration): PASS, |z|_p95 = 2.23. Continuous data fixes
  the Iter-2 issue.
- Sanity 2: PASS, |z| = 20.6 on Car2/Klf1 (erythroid TFs).
- Sanity 3: failed under literal `ρ > 0.7`, but |ρ| = 0.93 — the
  *sign-flip* between pseudotime (↑ with differentiation) and
  distance-to-terminal (↓ with differentiation) is an artefact of
  outcome construction, not biology.
- Sanity 4: 9/45 pairs underpowered (min corner < 30) due to
  gene-gene co-expression among top-10 (Mpo, Prtn3, Ctsg are all
  granulocyte markers).
- Sanity 5 (NEW, outcome decoupling): max |corr(outcome, gene)| = 0.48
  in top-10. Pseudotime correlates with lineage markers.

**Finding 3 — HARD vs SOFT CIRCULARITY.** Iter-1's circularity
(outcome formula explicitly mentions test gene) is *hard* — it
breaks the test. Iter-3's correlation is *soft*: outcome and test
gene share information about a hidden differentiation state, but
neither is computed from the other. Pseudotime depends on the
neighborhood graph built from many HVGs; no single test gene drives
it.

The two cases require different treatment:
- HARD: disqualifying. Outcome formula must not mention any test gene.
- SOFT: report-only, not gating, *provided* signal magnitude
  exceeds spurious-correlation amplification (Constraint 4 below).

### Iteration 4 — Synthetic null gate (final pre-reg attempt)

User-specified additional gate: 100 synthetic null pairs constructed
by permuting one gene's expression across cells. Three checks:

- N1: `|mean(ε_null)| < 0.1 × |mean(ε_observed)|`.
  PASS, mean = 2.4e-4 vs threshold 8.5e-3.
- N2: median SE_null in [0.5, 2.0] × median SE_observed.
  PASS, ratio = 0.67.
- N3: `|z_null|` distribution close to folded N(0,1) by KS test.
  FAIL, KS D = 0.38, p = 8.7e-14.

Plus dominance criterion (Constraint 4): median |z| of significant
pairs must exceed `max |corr(outcome, gene)| × √n_eval`.
- max |corr| = 0.476, √n = 52.2, threshold = 24.9.
- median |z|_significant = 8.03.
- FAIL: observed signal magnitude is below the worst-case spurious
  amplification threshold.

**Finding 4 — SOFT CORRELATION × √n MATCHES OBSERVED SIGNAL
MAGNITUDE.** With n = 2730 cells and max outcome–gene correlation
= 0.48, pure marginal-correlation amplification can produce |z|
values up to ≈ 25. The observed median significant |z| = 8 is well
below this ceiling, so cannot be claimed to dominate
spurious-correlation contribution.

N3's failure is structural to permutation-of-B null construction:
permuting B's values redistributes cells across high-B / low-B
randomly within each A-stratum, but the high-B and low-B subsets are
random splits of the *same* underlying group. Bootstrap (resampling
within each corner independently) over-estimates SE relative to the
true across-permutation variance. This is a bootstrap-method
limitation under permutation null specifically; bootstrap on real
random-pair nulls (Iter-3 Sanity 1) calibrates correctly.

## 4. Constraint synthesis

| # | Constraint | Source | Status |
|---|------------|--------|--------|
| 1 | Outcome formula must not mention test genes (HARD circularity) | Iter 1 | PASS-able, pre-conditional |
| 2 | Test-gene expression must vary continuously across cells | Iter 2 | PASS-able, dataset-conditional |
| 3 | Outcome–gene correlation reported (SOFT, transparent) | Iter 3 | always reportable |
| 4 | Signal dominance: median \|z\| > max \|corr\| × √n_eval | Iter 4 | failed on Paul15 |

Constraints 1–3 are pre-conditions: any conforming
dataset/outcome/test-set passes all three. Constraint 4 is the
*minimum-effect-size* gate: it requires that real biology produces
signal strength exceeding the spurious-correlation ceiling. Failure
on Constraint 4 is empirical, not methodological.

## 5. When 2×2 stratification IS valid

Combining the four constraints, the design is valid when *all* of:

- Outcome is computed from a feature set disjoint from the test set
  (Constraint 1).
- Test genes have continuous expression (q_low > 0 or, on
  scaled/normalized data, q_low ≠ q_high) — Constraint 2.
- Test-gene–outcome correlations are reported in the verdict for
  transparency (Constraint 3).
- Signal magnitude clears `max |corr| × √n` floor (Constraint 4).

In practice on **observational scRNA-seq with cell-state outcomes**
(pseudotime, fate probability, cluster identity), Constraint 4 is
the binding constraint: cell-state outcomes inevitably correlate
softly with lineage-marker test genes through the shared latent
differentiation signal, and the soft correlation × √n term grows
faster than typical second-order interaction signal magnitudes.
We did not identify an observational scRNA-seq dataset in the
scanpy-resident set on which Constraint 4 passes.

## 6. Path forward

**Direct experimental control breaks the soft-correlation barrier.**
In Perturb-seq (CRISPR-based pooled scRNA-seq), gene perturbations
are experimentally introduced; cells with a given perturbation are
identified by guide-RNA barcode, not by natural expression
stratification. The test becomes:
- Compare outcome between perturbation classes (genuine intervention
  contrast).
- Pair perturbations: cells receiving guides for both A and B
  (combinatorial Perturb-seq, e.g. Norman 2019) provide direct ε
  measurement.
- Outcome can be cluster-state or fate-probability without inducing
  spurious soft-correlation, because cell membership in pair classes
  is determined by guide barcodes (independent of expression).

**Constraint 4 in Perturb-seq:** soft correlation between perturbation
identity and outcome is design-controlled (perturbation efficiency,
not baseline expression heterogeneity). The √n amplification
collapses to noise scale rather than correlation × √n. Real ε signals
can clear the resulting (much lower) noise floor.

A scoping report (`perturbseq_scoping.md`) assesses Norman 2019,
Replogle 2022, and other candidate datasets for feasibility.

## 7. Implications

This pilot programme demonstrates:

1. The 2×2 stratification design is a valid methodological framework
   *under specific data and outcome conditions*. Constraints 1–3 are
   in-principle satisfiable; Constraint 4 is the practical limiter.

2. Observational scRNA-seq with cell-state outcomes reaches
   Constraint 4 floor. Tested datasets (Bastidas-Ponce E15.5 pancreas,
   Paul15 myeloid lineage) failed in distinct ways but converged on
   the same underlying issue: outcome and test-gene expression share
   information through latent cell state.

3. Future bio-epistasis work should use experimental perturbation
   data (Perturb-seq, single-cell CRISPR-screen variants) rather than
   observational stratification.

4. Discipline: pilot-before-pre-reg is essential. All four findings
   surfaced at pilot stage, before any decision-locking pre-reg was
   committed.

## 8. Connection to ML programme

The ML side of the parallel programme uses *deliberate ablation*
(zero or mean replacement of attention-head outputs), which IS a
direct intervention — analog of Perturb-seq, not of observational
stratification. The methodological asymmetry between ML and
observational biology is intrinsic: ML intervenes by construction,
biology requires CRISPR to match.

The methodological contribution of this document is therefore to
locate *where* the bio analog of ML mean ablation can be measured.
Not in observational data; in Perturb-seq. The cross-substrate
universality claim — synthetic-lethal/redundancy regime in top
functional units across both systems — remains testable, conditional
on Perturb-seq feasibility.

---

## References

- Costanzo M, et al. (2010) The Genetic Landscape of a Cell. *Science*
  327:425. — biological epistasis sign convention.
- Schiebinger G, et al. (2019) *Cell* 176. — wot reprogramming
  dataset (initially planned, URL since removed).
- Bastidas-Ponce A, et al. (2019) *Development* 146. — pancreas
  E15.5 dataset (Iter 1–2).
- Paul F, et al. (2015) *Cell* 163. — myeloid lineage dataset
  (Iter 3–4).
- Replogle JM, et al. (2022) *Cell* 185. — large-scale Perturb-seq
  reference for scoping.
- Norman TM, et al. (2019) *Science* 365. — combinatorial
  Perturb-seq, candidate primary dataset.
- Pythia 410M Tier 1 verdict (sister project, tag `tier1_pass`):
  cross-reference for ML side of the universality claim.

---

*2026-04-26. Status: internal documentation. Update if Iter 5+ on
Perturb-seq surfaces additional constraints.*
