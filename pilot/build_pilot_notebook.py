"""
Generate pilot/01_schiebinger_d8_10pairs.ipynb.

Pilot calibration of the BioEpistasis 2×2 design on Schiebinger day-8.
Sanity checks gate the pre-registration draft.

Outcome strategy in pilot:
- Primary  : iPSC-marker score (mean expression of pluripotency markers)
- Secondary: distance to iPSC-marker-high cells in PCA space
Both are PROXIES; full pre-reg run will use proper WOT trajectory
probability. Pilot's job is to calibrate bootstrap, sample size, and
outcome-consistency, not to produce trajectory results.

Build with:
    python pilot/build_pilot_notebook.py
"""

from __future__ import annotations

import json
import os

NB_PATH = os.path.join(os.path.dirname(__file__),
                       "01_schiebinger_d8_10pairs.ipynb")


def md(src: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": src}


def code(src: str) -> dict:
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": src}


cells: list[dict] = []


cells.append(md(r"""# BioEpistasis pilot — Schiebinger day 8, top-10 HVG, 45 pairs

**Locked methodology (per `design_notes.md`).**

- Sign convention: ε in **loss space**, additive null. ε > 0 ↔
  synthetic-lethal/redundancy; ε < 0 ↔ suppression/buffering.
- 2×2 stratification on tertile bins of A and B expression
  (analog of ML mean ablation = conditional-expectation operation).
- Outcome (pilot): iPSC-marker score. Secondary cross-check: PCA
  distance to iPSC-marker-high cells. Full pre-reg uses WOT log P.
- Bootstrap SE via cell-resampling within each corner, n=1000.
- Sanity gates (all four MUST pass to proceed to pre-reg):
  1. Random-pair bootstrap calibration: |z| ~ N(0,1) on null pairs.
  2. Curated same-pathway pair sensitivity: |z| > 3.
  3. Outcome consistency (marker vs distance): Pearson ρ > 0.7.
  4. Sample size per quadrant ≥ 30 across all 45 top-10 pairs.

Compute: ~30 min on Colab CPU. No GPU needed.
"""))


cells.append(md("""## 1. Setup + dependencies"""))
cells.append(code(r"""!pip install -q scanpy==1.10.3 anndata==0.10.9 \
                    pyarrow==16.1.0 matplotlib 2>&1 | tail -3"""))


cells.append(md("""## 2. Imports + sign-convention reminder

This block prints the sign convention to stdout before any
computation, so the notebook record always begins with the convention
locked in plain sight."""))

cells.append(code(r"""import os, time, json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt

print('=' * 70)
print('SIGN CONVENTION (locked, see design_notes.md §1):')
print('  ε > 0 ↔ synthetic-lethal / redundancy phenotype')
print('         (joint hurt MORE than additive — the dominant biology pattern')
print('          per Costanzo 2010, AND the dominant Pythia top-30 pattern)')
print('  ε < 0 ↔ suppression / buffering / true rescue (Costanzo "compensation")')
print('=' * 70)

SEED = 42
N_BOOT = 1000
np.random.seed(SEED)

# Output directory (Drive if Colab, else local)
try:
    from google.colab import drive
    drive.mount('/content/drive', force_remount=False)
    OUT_DIR = '/content/drive/MyDrive/BioEpistasis_pilot'
except (ImportError, ModuleNotFoundError):
    OUT_DIR = os.path.expanduser('~/BioEpistasis_pilot')
os.makedirs(OUT_DIR, exist_ok=True)
print(f'OUT_DIR: {OUT_DIR}')"""))


# ─────────────────────────────────────────────────────────────────────────────
# Data
# ─────────────────────────────────────────────────────────────────────────────

cells.append(md("""## 3. Download + load Schiebinger day-8 subset

The wot tutorial dataset (`ExprMatrix.h5ad`, ~1 GB) covers all 18
timepoints. We load it lazily and filter to day 8 immediately to avoid
holding the full matrix in memory."""))

cells.append(code(r"""DATA_DIR = os.path.join(OUT_DIR, 'data')
os.makedirs(DATA_DIR, exist_ok=True)
H5AD = os.path.join(DATA_DIR, 'ExprMatrix.h5ad')

if not os.path.exists(H5AD):
    print('downloading Schiebinger ExprMatrix.h5ad (~1 GB) ...')
    !wget -q -O {H5AD} https://broadinstitute.github.io/wot/example_data/ExprMatrix.h5ad
    print('done.')

size_mb = os.path.getsize(H5AD) / 1e6
print(f'ExprMatrix.h5ad: {size_mb:.0f} MB')
assert size_mb > 100, 'h5ad download incomplete — re-run cell'"""))

cells.append(code(r"""adata_full = sc.read_h5ad(H5AD)
print(f'full AnnData: {adata_full}')
print(f'obs columns: {list(adata_full.obs.columns)}')

# Strategy: try day info in adata.obs first; if absent, download cell_days.txt
# (correct wot tutorial filename — NOT days.txt).
day_col_candidates = [c for c in adata_full.obs.columns if 'day' in c.lower()]
if day_col_candidates:
    day_col = day_col_candidates[0]
    print(f'day column found in obs: {day_col}')
else:
    CELL_DAYS = os.path.join(DATA_DIR, 'cell_days.txt')
    if not os.path.exists(CELL_DAYS) or os.path.getsize(CELL_DAYS) < 1024:
        print('downloading cell_days.txt ...')
        !wget -q -O {CELL_DAYS} https://broadinstitute.github.io/wot/example_data/cell_days.txt
    size_kb = os.path.getsize(CELL_DAYS) / 1024
    print(f'cell_days.txt: {size_kb:.0f} KB')
    assert size_kb > 1, 'cell_days.txt download failed — check URL'

    # File format: header row "id\tday", then tab-separated rows
    days_df = pd.read_csv(CELL_DAYS, sep='\t', index_col=0)
    print(f'days_df rows: {len(days_df)}, columns: {list(days_df.columns)}')
    adata_full.obs = adata_full.obs.join(days_df, how='left')
    day_col = days_df.columns[0]
    print(f'day column joined into obs: {day_col}')

unique_days = sorted(adata_full.obs[day_col].dropna().unique())
print(f'unique days ({len(unique_days)}): {unique_days[:8]} ... {unique_days[-3:]}')

# Filter to day 8
day_target = 8.0
mask_d8 = adata_full.obs[day_col] == day_target
print(f'day-8 mask: {int(mask_d8.sum())} cells')
adata = adata_full[mask_d8].copy()
print(f'subset AnnData: {adata}')
del adata_full"""))


cells.append(md("""## 4. Standard preprocessing + HVG selection (seurat_v3)"""))
cells.append(code(r"""# Schiebinger ExprMatrix is already log-normalized in the wot tutorial,
# but seurat_v3 HVG flavor expects RAW counts. We use the cell_ranger
# flavor instead, which works on log-normalized data and is the next-best
# choice for this dataset.
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=200)
print(f'HVGs flagged: {adata.var["highly_variable"].sum()}')

# Top-10 HVG ranked by dispersion_norm (seurat_v3 analog)
hvg_df = adata.var[adata.var['highly_variable']].copy()
hvg_df = hvg_df.sort_values('dispersions_norm', ascending=False)
TOP10 = hvg_df.head(10).index.tolist()
print(f'top-10 HVGs: {TOP10}')"""))


cells.append(md("""## 5. Outcome computation: iPSC-marker score (primary, pilot)"""))

cells.append(code(r"""# Pluripotency markers used in Schiebinger's analysis.
# Names cap-cased to match the mouse gene symbols in the dataset.
IPSC_MARKERS = ['Pou5f1', 'Nanog', 'Sox2', 'Klf4', 'Esrrb', 'Zfp42']
present = [g for g in IPSC_MARKERS if g in adata.var_names]
print(f'iPSC markers present: {present}')
assert len(present) >= 3, 'too few iPSC markers found — check gene name casing'

# Marker score = mean log-normalized expression of present markers
sc.tl.score_genes(adata, gene_list=present, score_name='ipsc_score', random_state=SEED)

# Outcome convention: HIGHER outcome = "WORSE" (analog of loss).
# Define outcome = - ipsc_score (lower marker → higher outcome → "loss").
adata.obs['outcome_marker'] = -adata.obs['ipsc_score']
print(f'outcome_marker: mean={adata.obs["outcome_marker"].mean():.3f}, '
      f'std={adata.obs["outcome_marker"].std():.3f}')"""))


cells.append(md("""## 6. Outcome computation: PCA distance to iPSC-high (secondary)"""))
cells.append(code(r"""# PCA on the day-8 subset
sc.pp.pca(adata, n_comps=50, random_state=SEED)

# iPSC-high cells = top quartile by ipsc_score
ipsc_high_mask = adata.obs['ipsc_score'] >= adata.obs['ipsc_score'].quantile(0.75)
ipsc_centroid = adata.obsm['X_pca'][ipsc_high_mask].mean(axis=0)

# Euclidean distance to centroid in 50D PCA
diff = adata.obsm['X_pca'] - ipsc_centroid[None, :]
adata.obs['outcome_distance'] = np.linalg.norm(diff, axis=1)
print(f'outcome_distance: mean={adata.obs["outcome_distance"].mean():.3f}, '
      f'std={adata.obs["outcome_distance"].std():.3f}')
print(f'iPSC-high cells (top quartile): {ipsc_high_mask.sum()}')"""))


# ─────────────────────────────────────────────────────────────────────────────
# Core: 2×2 stratification + bootstrap ε
# ─────────────────────────────────────────────────────────────────────────────

cells.append(md("""## 7. 2×2 stratification + bootstrap ε helper

For pair (A, B): tertile-cut both, take the four corner quadrants,
median outcome per corner, ε = Δ_AB − Δ_A − Δ_B.

Bootstrap resamples cells within each corner with replacement, recomputes
ε, std of the bootstrap distribution = SE."""))

cells.append(code(r"""def stratify_2x2(adata, gene_a, gene_b, outcome_col,
                  q_low=1/3, q_high=2/3):
    '''Return (cells_HH, cells_HL, cells_LH, cells_LL) DataFrames where
    H/L are tertile bins of expression. None if any corner has 0 cells.'''
    if gene_a not in adata.var_names or gene_b not in adata.var_names:
        return None
    expr_a = adata[:, gene_a].X
    expr_b = adata[:, gene_b].X
    if hasattr(expr_a, 'toarray'):
        expr_a = expr_a.toarray()
    if hasattr(expr_b, 'toarray'):
        expr_b = expr_b.toarray()
    expr_a = np.asarray(expr_a).flatten()
    expr_b = np.asarray(expr_b).flatten()

    a_low_thr  = np.quantile(expr_a, q_low)
    a_high_thr = np.quantile(expr_a, q_high)
    b_low_thr  = np.quantile(expr_b, q_low)
    b_high_thr = np.quantile(expr_b, q_high)

    out = adata.obs[outcome_col].values

    # Boolean masks
    a_high = expr_a >= a_high_thr
    a_low  = expr_a <= a_low_thr
    b_high = expr_b >= b_high_thr
    b_low  = expr_b <= b_low_thr

    return {
        'HH': out[a_high & b_high],
        'HL': out[a_high & b_low],
        'LH': out[a_low  & b_high],
        'LL': out[a_low  & b_low],
    }


def bootstrap_epsilon(corners, n_boot=1000, seed=42):
    '''Bootstrap ε for a 2×2 stratification.

    Δ_A  = median(LH) − median(HH)   (A "perturbed" — expression LOW)
    Δ_B  = median(HL) − median(HH)
    Δ_AB = median(LL) − median(HH)
    ε    = Δ_AB − Δ_A − Δ_B
    '''
    n = {k: len(v) for k, v in corners.items()}
    if min(n.values()) < 10:
        return None  # under-sampled

    rng = np.random.default_rng(seed)

    def stat(c):
        return (np.median(c['LL']) - np.median(c['HH'])) \
             - (np.median(c['LH']) - np.median(c['HH'])) \
             - (np.median(c['HL']) - np.median(c['HH']))

    eps_obs = stat(corners)

    boots = np.empty(n_boot)
    for i in range(n_boot):
        resampled = {k: rng.choice(v, size=len(v), replace=True)
                     for k, v in corners.items()}
        boots[i] = stat(resampled)

    se = float(boots.std(ddof=1))
    z = eps_obs / se if se > 0 else 0.0
    return {
        'epsilon':  float(eps_obs),
        'se':       se,
        'z':        float(z),
        'n_corners': n,
    }"""))


cells.append(md("""## 8. Run on top-10 HVG → 45 pairs (both outcomes)"""))

cells.append(code(r"""from itertools import combinations

PAIRS_TOP10 = list(combinations(TOP10, 2))
print(f'top-10 pair count: {len(PAIRS_TOP10)}')

results_top10 = []
for (A, B) in PAIRS_TOP10:
    rec = {'gene_a': A, 'gene_b': B}
    for outcome_name, col in [('marker', 'outcome_marker'),
                              ('distance', 'outcome_distance')]:
        corners = stratify_2x2(adata, A, B, col)
        if corners is None:
            rec[f'{outcome_name}_skip'] = True
            continue
        boot = bootstrap_epsilon(corners, n_boot=N_BOOT, seed=SEED)
        if boot is None:
            rec[f'{outcome_name}_skip'] = True
            continue
        rec[f'{outcome_name}_eps']  = boot['epsilon']
        rec[f'{outcome_name}_se']   = boot['se']
        rec[f'{outcome_name}_z']    = boot['z']
        rec[f'{outcome_name}_nHH']  = boot['n_corners']['HH']
        rec[f'{outcome_name}_nHL']  = boot['n_corners']['HL']
        rec[f'{outcome_name}_nLH']  = boot['n_corners']['LH']
        rec[f'{outcome_name}_nLL']  = boot['n_corners']['LL']
    results_top10.append(rec)
df_top = pd.DataFrame(results_top10)
df_top.to_parquet(os.path.join(OUT_DIR, 'top10_pairs.parquet'))
print(df_top.head().to_string(index=False))"""))


# ─────────────────────────────────────────────────────────────────────────────
# Sanity checks
# ─────────────────────────────────────────────────────────────────────────────

cells.append(md("""## 9. Sanity check 1 — bootstrap calibration on random pairs

50 random gene pairs (sampled from full gene pool, NOT HVG). Expect
|z| ~ N(0,1) under null of no real epistasis. PASS if 95th percentile
of |z| is in [1.5, 2.5] (theoretical 1.96)."""))

cells.append(code(r"""rng_random = np.random.default_rng(SEED)
all_genes = adata.var_names.tolist()
n_random = 50
random_pairs = []
seen = set()
while len(random_pairs) < n_random:
    i, j = rng_random.integers(0, len(all_genes), size=2)
    if i == j: continue
    a, b = all_genes[int(i)], all_genes[int(j)]
    if a == b or (a, b) in seen or (b, a) in seen: continue
    seen.add((a, b))
    random_pairs.append((a, b))

random_zs = []
for (A, B) in random_pairs:
    corners = stratify_2x2(adata, A, B, 'outcome_marker')
    if corners is None: continue
    boot = bootstrap_epsilon(corners, n_boot=N_BOOT, seed=SEED)
    if boot is None: continue
    random_zs.append(boot['z'])
random_zs = np.array(random_zs)
abs_z_p95 = float(np.percentile(np.abs(random_zs), 95))
print(f'random pairs measured: {len(random_zs)}')
print(f'|z| 95th percentile  : {abs_z_p95:.2f}  (target [1.5, 2.5])')
SANITY1_PASS = 1.5 <= abs_z_p95 <= 2.5
print(f'Sanity #1 (calibration): {"PASS" if SANITY1_PASS else "FAIL"}')"""))


cells.append(md("""## 10. Sanity check 2 — curated same-pathway pair sensitivity

Pre-registered curated pair: two ribosomal proteins. Both expressed
broadly, both in shared pathway (translation). Predict |z| > 3
(detectable epistasis on this proxy outcome)."""))

cells.append(code(r"""# Try a few candidate ribosomal pairs in order; use the first present.
CURATED_PAIRS = [('Rps3', 'Rps5'), ('Rpl3', 'Rpl5'), ('Rps3', 'Rpl5'),
                 ('Eef1a1', 'Eef2'), ('Actb', 'Tubb5')]
sanity2_pair = None
sanity2_result = None
for (A, B) in CURATED_PAIRS:
    if A in adata.var_names and B in adata.var_names:
        corners = stratify_2x2(adata, A, B, 'outcome_marker')
        if corners is None: continue
        boot = bootstrap_epsilon(corners, n_boot=N_BOOT, seed=SEED)
        if boot is None: continue
        sanity2_pair = (A, B)
        sanity2_result = boot
        break

assert sanity2_pair is not None, 'no curated pair worked — investigate gene names'
A, B = sanity2_pair
print(f'curated pair tested: ({A}, {B})')
print(f'  ε  = {sanity2_result["epsilon"]:+.5f}')
print(f'  SE = {sanity2_result["se"]:.5f}')
print(f'  z  = {sanity2_result["z"]:+.2f}')
SANITY2_PASS = abs(sanity2_result['z']) > 3
print(f'Sanity #2 (sensitivity): {"PASS" if SANITY2_PASS else "FAIL"}')"""))


cells.append(md("""## 11. Sanity check 3 — outcome consistency

Pearson ρ between ε computed from marker outcome and ε from distance
outcome, across the 45 top-10 pairs. PASS if ρ > 0.7 (outcomes agree
on epistasis direction). FAIL → outcome choice is fragile."""))

cells.append(code(r"""eps_marker  = df_top['marker_eps'].values
eps_dist    = df_top['distance_eps'].values
mask = np.isfinite(eps_marker) & np.isfinite(eps_dist)
rho_consist, p_consist = pearsonr(eps_marker[mask], eps_dist[mask])
print(f'pairs with both ε computed: {mask.sum()}/{len(df_top)}')
print(f'Pearson ρ(marker, distance): {rho_consist:+.3f}  (p={p_consist:.4f})')
SANITY3_PASS = rho_consist > 0.7
print(f'Sanity #3 (outcome consistency): {"PASS" if SANITY3_PASS else "FAIL"}')

# Visualise
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(eps_marker[mask], eps_dist[mask], alpha=0.6)
ax.axhline(0, color='gray', lw=0.5)
ax.axvline(0, color='gray', lw=0.5)
ax.set_xlabel('ε (marker-score outcome)')
ax.set_ylabel('ε (PCA-distance outcome)')
ax.set_title(f'Outcome consistency on 45 top-10 pairs (ρ={rho_consist:.3f})')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'sanity3_outcome_consistency.png'), dpi=120)
plt.show()"""))


cells.append(md("""## 12. Sanity check 4 — sample size per quadrant

Across all 45 top-10 pairs, every corner of the 2×2 must have ≥ 30
cells. PASS if min over all corners ≥ 30. FAIL → tertile split is
too aggressive at this dataset scale."""))

cells.append(code(r"""corner_cols = [c for c in df_top.columns if c.startswith('marker_n')]
corner_mins = df_top[corner_cols].min(axis=1)
overall_min = int(corner_mins.min())
n_below_30  = int((corner_mins < 30).sum())
print(f'min cells/corner across all 45 pairs: {overall_min}')
print(f'pairs with any corner < 30 cells   : {n_below_30}/{len(df_top)}')
SANITY4_PASS = (overall_min >= 30) and (n_below_30 == 0)
print(f'Sanity #4 (sample size): {"PASS" if SANITY4_PASS else "FAIL"}')"""))


# ─────────────────────────────────────────────────────────────────────────────
# Verdict
# ─────────────────────────────────────────────────────────────────────────────

cells.append(md("""## 13. Calibration verdict — gate to pre-registration"""))

cells.append(code(r"""sanity_results = {
    'sanity1_calibration_p95': abs_z_p95,
    'sanity1_pass':            bool(SANITY1_PASS),
    'sanity2_curated_pair':    list(sanity2_pair),
    'sanity2_z':               sanity2_result['z'],
    'sanity2_pass':            bool(SANITY2_PASS),
    'sanity3_pearson_rho':     float(rho_consist),
    'sanity3_pass':            bool(SANITY3_PASS),
    'sanity4_min_corner_cells': overall_min,
    'sanity4_pass':            bool(SANITY4_PASS),
}
all_pass = all([SANITY1_PASS, SANITY2_PASS, SANITY3_PASS, SANITY4_PASS])

verdict = {
    'project':           'BioEpistasis',
    'phase':             'Pilot — calibration on Schiebinger day 8',
    'dataset':           'Schiebinger 2019 reprogramming MEF→iPSC',
    'timepoint':         'day 8 (mid-trajectory)',
    'n_cells':           int(adata.n_obs),
    'top10_hvg':         TOP10,
    'n_top10_pairs':     len(PAIRS_TOP10),
    'n_random_pairs':    n_random,
    'outcome_primary':   'iPSC-marker score (mean log-norm pluripotency markers)',
    'outcome_secondary': 'PCA distance to iPSC-marker-high centroid',
    'sign_convention':   'ε in loss space; ε>0 = synthetic-lethal/redundancy, ε<0 = suppression',
    'sanity':            sanity_results,
    'all_sanity_pass':   bool(all_pass),
    'next_step':         'Draft pre-registration v1 → lock → full Schiebinger trajectory'
                          if all_pass else
                          'Halt. Address failing sanity check before pre-reg drafting.',
}

out = os.path.join(OUT_DIR, 'pilot_calibration_verdict.json')
with open(out, 'w') as f:
    json.dump(verdict, f, indent=2)
print(json.dumps(verdict, indent=2))
print()
print('=' * 70)
print('PILOT CALIBRATION VERDICT:',
      'PASS — proceed to pre-reg drafting' if all_pass
      else 'FAIL — halt, address sanity issues before pre-reg')
print('=' * 70)
print(f'\\nSaved → {out}')"""))


nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                   "language_info": {"name": "python"},
                   "colab": {"provenance": []}},
      "nbformat": 4, "nbformat_minor": 5}


def main() -> None:
    with open(NB_PATH, "w") as f:
        json.dump(nb, f, indent=1)
    print(f"Wrote {NB_PATH}  ({len(cells)} cells)")


if __name__ == "__main__":
    main()
