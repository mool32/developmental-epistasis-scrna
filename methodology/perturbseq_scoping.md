# Perturb-seq dataset scoping for BioEpistasis pivot

*2026-04-26. Internal scoping report. Assesses candidate Perturb-seq
datasets against the 2×2 epistasis design requirements identified in
`observational_epistasis_limits.md`. Output: feasibility verdict per
dataset and overall recommendation.*

---

## Goal

Identify a Perturb-seq dataset on which the 2×2 epistasis design is
applicable WITHOUT the observational-data limitations (HARD/SOFT
circularity, Constraint 4 dominance failure). The candidate dataset
must satisfy:

1. **Combinatorial perturbations.** Either explicit pair perturbations
   (cells receiving guides for both A and B simultaneously), or
   sufficient density of single-perturbation cells to construct
   counterfactual joint-knockout estimates via interaction modelling.
2. **Cell counts adequate for bootstrap SE.** ≥ ~30 cells per
   perturbation pair class (matching the Iter-4 S4 floor).
3. **Outcome variable computable.** Cluster identity, fate probability,
   distance from a reference state — analog of the ML loss.
4. **Accessible.** Public download, standard format (h5ad / mtx / GEO).

---

## Dataset 1 — Norman et al. 2019 *Science*

**Reference.** Norman TM, et al. (2019) "Exploring genetic interaction
manifolds constructed from rich single-cell phenotypes." *Science*
365:786–793.

**Source.** GEO GSE133344. Accompanying GitHub:
[github.com/thomasmaxwellnorman/perturbseq_GI](https://github.com/thomasmaxwellnorman/perturbseq_GI)

**Why this matters most.** This is the *gold standard combinatorial
Perturb-seq study*: it explicitly delivers (single, single, pair)
perturbations of curated gene pairs. The paper's own analysis is
genetic-interaction inference — exactly the question we're asking.

**Key parameters:**
- Cell line: K562 (chronic myeloid leukaemia, well-characterised).
- Perturbation method: CRISPRi (dCas9-KRAB) — knocks down
  transcription, doesn't fully knock out → *gradient* of perturbation
  efficiency that mirrors mean-ablation more than zero-ablation.
- Genes perturbed: 287 unique genes × ~150 pairwise combinations,
  curated for putative interaction relevance. Pair coverage is
  designed for genetic-interaction analysis, not full combinatorial.
- Cell count: ~110,000 cells across all conditions, ~150 cells per
  pair condition.
- Format: published as h5 / loom / available via Norman et al. python
  package.

**Suitability for our 2×2 analog:**
- Direct match: cells assigned to "perturbed A only" / "perturbed B
  only" / "perturbed both" / "control" by guide barcodes. Maps cleanly
  onto our HH/HL/LH/LL semantics.
- Outcome: cluster identity in latent space (paper uses UMAP+Leiden);
  also cell-cycle phase, fitness via cell viability proxy.
- ~150 cells per pair class is *tight* but sufficient for bootstrap
  with the Iter-4 S4 floor of 30 (margin of 5×).
- Perturbation gradient (CRISPRi knockdown, not full knockout) is a
  *better* analog of ML mean ablation than CRISPR-knockout would be
  — both are "remove signal, don't fully zero".
- Genes are pre-curated by paper authors for interaction relevance.
  Selection bias is a feature, not a bug, for our purposes — analog
  of choosing top-K heads by |Δ| in ML.

**Feasibility:** **HIGH.** Dataset is purpose-built for our question,
public, well-documented, modest cell count, with the right
perturbation modality (CRISPRi knockdown).

**Preprocessing effort:** ~3–5 days. The Norman team's own python
package handles loading and preprocessing; integration with our 2×2
analysis pipeline is ~1 day on top.

---

## Dataset 2 — Replogle et al. 2022 *Cell*

**Reference.** Replogle JM, et al. (2022) "Mapping information-rich
genotype-phenotype landscapes with genome-scale Perturb-seq." *Cell*
185:2559–2575.

**Source.** Figshare:
[plus.figshare.com/articles/dataset/_/20029387](https://plus.figshare.com/articles/dataset/_/20029387).
Mirror on GitHub-hosted analysis repo.

**Key parameters:**
- Cell line: K562 (same as Norman, allows cross-study comparison).
- Perturbation method: CRISPRi.
- Genes perturbed: ~9,800 (genome-scale single-gene). **No paired
  combinations** in the published dataset.
- Cell count: ~2.5 million across all conditions, ~250 cells per
  single-gene perturbation.
- Format: h5ad (~30 GB compressed, ~100 GB unpacked).

**Suitability:**
- *Single-gene only* — the design lacks the pair-class structure we
  need for direct ε measurement. ε can in principle be inferred via
  *interaction-on-singles* statistical models (each cell carries
  one guide, pair effects estimated from main effects + interaction
  term), but that requires an additional modelling layer and is not
  the same instrument as the 2×2 design.
- Cell counts massive — single-gene power is excellent.
- Outcome: paper proposes a "perturbation phenotype" continuous
  metric; usable as outcome.

**Feasibility for OUR 2×2 design:** **LOW.** The dataset is
single-gene-only. The 2×2 design assumes joint-perturbation
observability.

**Alternative use:** valuable for *single-gene* single-ablation
DFE measurement (analog of Phase 2B in ML side). Could complement
Norman 2019 if both are used: Norman for pairs, Replogle for
single-gene Δ ranking.

**Preprocessing effort:** 5–7 days for full ingestion (large data,
nontrivial format). Reduced if used only for single-gene Δ.

---

## Dataset 3 — Adamson et al. 2016 (UPR Perturb-seq)

**Reference.** Adamson B, et al. (2016) "A multiplexed single-cell
CRISPR screening platform enables systematic dissection of the
unfolded protein response." *Cell* 167:1867–1882.

**Source.** GEO GSE90546.

**Parameters:**
- Cell line: K562.
- Perturbation method: CRISPRi.
- Genes: ~91 single-gene perturbations focused on the UPR pathway.
  Some pair combinations (~40 pairs, low coverage).
- Cell count: ~65,000 cells.

**Suitability:** **MEDIUM-LOW.** Pair coverage is small (40 pairs)
and confined to one pathway. Can serve as a *replication* dataset
for any Norman-2019 finding within the UPR domain, but not as
primary.

**Preprocessing:** ~3 days.

---

## Dataset 4 — Dixit et al. 2016 (original Perturb-seq)

**Reference.** Dixit A, et al. (2016) "Perturb-seq: dissecting
molecular circuits with scalable single-cell RNA profiling of pooled
genetic screens." *Cell* 167:1853–1866.

**Source.** GEO GSE90063.

**Parameters:** ~24 single-gene perturbations, ~20,000 cells. No
combinatorial design.

**Suitability:** **LOW.** Too few genes. Single-gene only. Useful as
a method reference, not a primary dataset.

---

## Dataset 5 — Gasperini et al. 2019 (enhancer Perturb-seq)

**Reference.** Gasperini M, et al. (2019) "A Genome-wide Framework
for Mapping Gene Regulation via Cellular Genetic Screens." *Cell*
176:377–390.

**Parameters:** ~5,000 enhancer perturbations, ~250,000 cells.
Single-perturbation only, focus on enhancer–gene mapping.

**Suitability:** **NOT APPLICABLE.** Enhancer-perturbation, not
gene-knockdown. Different question.

---

## Comparison

| Dataset | Pair design | Cells/pair class | Genes (pairs) | Effort | Feasibility |
|---------|-------------|------------------|---------------|--------|-------------|
| **Norman 2019** | **explicit** | ~150 | 287 (~150 pairs) | 3–5 d | **HIGH** |
| Replogle 2022 | single-gene only | n/a | 9800 (0 pairs) | 5–7 d | LOW (for 2×2) |
| Adamson 2016 | sparse pairs | ~100 | 91 (~40 pairs) | 3 d | medium-low |
| Dixit 2016 | single-only | n/a | 24 (0 pairs) | 2 d | low |
| Gasperini 2019 | enhancer | n/a | n/a | n/a | n/a |

---

## Recommendation

**Primary dataset: Norman 2019.** It is purpose-built for the question
we want to ask (genetic interactions in pooled scRNA-seq); its
~150 cells per pair condition clears the Iter-4 sample-size floor
with margin; CRISPRi knockdown is the cleaner ML-mean-ablation analog
than full knockout would be; and its 150-pair design directly
populates the Tier-1-equivalent test set. **Replogle 2022 is the
secondary** dataset — for single-gene Δ ranking (analog of Phase 2B
in ML side), to complement the pair design from Norman.

**New pilot iteration on Norman 2019** is required regardless,
because the design space is genuinely new: every constraint surfaced
in Iter 1–4 is observational-data-specific and may or may not apply
to perturbation data. Expect:

- Constraint 1 (HARD circularity): satisfied by construction —
  cells classified by guide barcode, not gene expression.
- Constraint 2 (continuous variation): N/A — perturbation is binary
  (guide present / absent), the 2×2 corners are perturbation classes
  not expression bins.
- Constraint 3 (SOFT correlation): negligible — guide barcode is
  independent of post-perturbation expression by experimental design.
- Constraint 4 (dominance): expected to PASS, because soft-correlation
  amplification term collapses (max |corr| in this regime is
  determined by experimental noise, not cell-state coupling).

**New constraints expected from Perturb-seq pilot:**

- Perturbation efficiency variability across cells (CRISPRi knockdown
  is incomplete; some cells may retain target expression).
- Guide misassignment in pooled experiments (~5–10 % typical).
- Cell-cycle and other state confounds.

These are different from observational issues, addressable via
existing Perturb-seq analysis tools (Mixscape for guide-call quality,
hierarchical models for perturbation efficiency).

## Plan

If Teo approves the Norman pivot:

1. **Pilot iteration 5 (Norman 2019).** Replicate Iter-4 sanity gates
   against the new design. Expected to PASS Constraints 1–4. New
   Perturb-seq-specific constraints likely to surface; document them
   in `methodological_findings.md` as Constraints 5+.
2. **Pre-reg v1 (Norman primary).** Lock primary statistic
   (analog of ratio test in ML Tier 1, possibly recast as effect-
   size-ratio with proper experimental error model).
3. **Replication via Replogle (singles) or Adamson (UPR pairs)** as
   secondary, after primary verdict.

**Timeline.** 1 week for Norman pilot + scoping, 1 week for pre-reg
draft + lock, 1–2 weeks for primary scan and analysis. Total ≈ 4
weeks to first verdict.

**Compute.** Norman 2019 fits comfortably on Colab CPU. No GPU
needed for ε computation; bootstrap is the bottleneck and is
embarrassingly parallel. Estimate: ~2 hours per full pair scan.

---

## Open questions for Teo

1. **Outcome definition on Norman.** Two principled choices:
   (a) cluster identity in 50-PC space (binary: cell in target
   transcriptional state) — robust, interpretable; (b) per-cell
   "perturbation phenotype" magnitude as defined in the Norman paper
   — more sensitive but introduces paper-specific statistic.
   Recommend (a) for first pilot.

2. **Pair selection.** Use Norman's pre-curated pairs (~150 pairs)
   or sub-select top-30 by single-perturbation magnitude (analog of
   top-K in ML)? Recommend top-30: matches ML protocol exactly,
   reduces multiple-testing burden, focuses on pairs most likely to
   show real ε.

3. **Sign-convention application.** Norman 2019 internally uses
   their own GI-score formulation, which is *fitness-space*. Our
   convention is loss-space. Translation is mechanical (negate
   sign), but verdict reports must be explicit about which space we
   report in.

4. **Single-gene extension via Replogle.** If pair-scan PASSes on
   Norman, do we also run the single-gene Δ scan on Replogle as
   replication-of-Phase-2B? Replogle is heavier (100 GB+) but
   provides independent single-gene scaffold.

---

*Status: scoping complete, pivot decision pending. Awaiting Teo
approval to switch primary to Norman 2019.*
