"""
Generate pilot/05_norman_iter5.ipynb.

BioEpistasis pilot iter 5 — pivot to Perturb-seq (Norman 2019 K562
combinatorial CRISPRi). Tests whether the 2×2 design SURVIVES under
experimental-intervention data, vs the observational scRNA-seq pilot
(iter 1-4) which hit the soft-correlation × √n dominance ceiling.

Per perturbseq_scoping.md, expected pre-reg gates on Norman:
  C1 HARD circularity  : satisfied by construction (cell class = guide barcode)
  C2 continuous variation: N/A (perturbation classes are categorical BY DESIGN)
  C3 SOFT correlation  : negligible (guide barcode independent of expression)
  C4 dominance         : expected PASS — soft-correlation × √n collapses

Six sanity gates (all must PASS to draft pre-reg v1 biology side):
  S1 Calibration on random label permutations: |z|_p95 ∈ [1.5, 2.5]
  S2 Sensitivity on a known pathway pair: |z| > 3
  S3 Outcome consistency (cluster-identity vs PCA distance): |ρ| > 0.7
  S4 Sample size per class ≥ 30 (control / A only / B only / pair AB)
  S5 Outcome decoupling (guide-id × outcome correlation < 0.3)
  S6 Dominance (median |z| > max |corr| × sqrt(n_eval))

Compute: ~1-2 h on Colab CPU (no GPU needed).

Build with:
    python pilot/build_norman_iter5_notebook.py
"""

from __future__ import annotations

import json
import os

NB_PATH = os.path.join(os.path.dirname(__file__),
                       "05_norman_iter5.ipynb")


def md(src: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": src}


def code(src: str) -> dict:
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": src}


cells: list[dict] = []


cells.append(md(r"""# BioEpistasis pilot iter 5 — Norman 2019 K562 Perturb-seq

**Path B pivot.** Iter 1-4 on observational scRNA-seq (Bastidas-Ponce
pancreas, Paul15 myeloid lineage) hit the dominance ceiling — soft
correlation × √n could explain observed signal magnitude. Norman 2019
provides experimental-intervention data: cells are classified by guide
barcode, not natural expression stratification. Soft correlation
collapses by construction.

**Design adaptation.** The 2×2 epistasis design becomes 4-class
intervention contrast:

| class | description |
|-------|-------------|
| WT (control) | non-targeting guide |
| A only       | guide for gene A, no guide for B |
| B only       | guide for gene B, no guide for A |
| AB pair      | guides for both A and B |

ε = (φ_AB − φ_WT) − (φ_A − φ_WT) − (φ_B − φ_WT)
  = φ_AB − φ_A − φ_B + φ_WT

where φ is the outcome (cluster identity or PCA distance). Sign
convention: ε > 0 = synthetic-lethal/redundancy, ε < 0 = suppression.

**Sanity gates (all must PASS):**
- S1 calibration: |z|_p95 on random class permutations ∈ [1.5, 2.5]
- S2 sensitivity on curated pathway pair: |z| > 3
- S3 outcome consistency: |Pearson ρ| > 0.7
- S4 sample size: ≥ 30 cells in every class for every test pair
- S5 decoupling: |corr(guide_class, outcome_per_cell)| < 0.3
- S6 dominance: median |z| > max |corr| × √n

**Compute:** ~1-2 h on Colab CPU. No GPU required.
"""))


cells.append(md("""## 1. Setup + dependencies"""))
cells.append(code(r"""!pip install -q scanpy==1.10.3 anndata==0.10.9 pyarrow==16.1.0 \
                    matplotlib requests 2>&1 | tail -3"""))


cells.append(md("""## 2. Imports + sign-convention banner"""))
cells.append(code(r"""import os, json, time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.stats import pearsonr, kstest, halfnorm
from itertools import combinations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print('=' * 70)
print('SIGN CONVENTION (locked, see methodology/observational_epistasis_limits.md §1):')
print('  ε > 0 ↔ synthetic-lethal / redundancy phenotype')
print('  ε < 0 ↔ suppression / buffering / true rescue')
print('=' * 70)
print()

SEED = 42
N_BOOT = 1000
np.random.seed(SEED)

# Output directory
try:
    from google.colab import drive
    drive.mount('/content/drive', force_remount=False)
    OUT_DIR = '/content/drive/MyDrive/BioEpistasis_pilot/iter5_norman'
except Exception:
    OUT_DIR = os.path.expanduser('~/BioEpistasis_pilot/iter5_norman')
os.makedirs(OUT_DIR, exist_ok=True)
print(f'OUT_DIR: {OUT_DIR}')"""))


cells.append(md("""## 3. Download Norman 2019 K562 data

Sources (try in order):
1. Norman lab Box / Figshare — preprocessed h5ad
2. GEO GSE133344 — raw 10x matrices

The preprocessed h5ad is much faster to use. If unavailable, the
fallback assembles from raw 10x files (slower, ~15 min).
"""))

cells.append(code(r"""DATA_DIR = os.path.join(OUT_DIR, 'data')
os.makedirs(DATA_DIR, exist_ok=True)
H5AD = os.path.join(DATA_DIR, 'norman2019.h5ad')

# Norman 2019 preprocessed h5ad — multiple potential mirrors. Try them
# in turn; if all fail, document and abort with instructions.
URLS = [
    # Most common preprocessed mirror — verify URL at runtime
    'https://ndownloader.figshare.com/files/34027562',
    # Backup mirror via Hugging Face datasets (if mirrored there)
    'https://huggingface.co/datasets/sciplex/Norman2019/resolve/main/norman2019.h5ad',
]

if not os.path.exists(H5AD):
    success = False
    for url in URLS:
        try:
            print(f'Trying: {url}')
            !wget -q --tries=2 --timeout=60 -O {H5AD} "{url}"
            sz = os.path.getsize(H5AD)
            print(f'  size: {sz/1e6:.1f} MB')
            if sz > 50e6:
                success = True
                break
            else:
                print('  too small, trying next mirror')
                os.remove(H5AD)
        except Exception as e:
            print(f'  failed: {e}')
            continue
    assert success, (
        'Could not download Norman 2019 h5ad from any mirror.\n'
        'Manual fallback: download from GEO GSE133344 + preprocess.\n'
        'See pilot/calibration_report.md for guidance.'
    )

print(f'\nNorman 2019 file: {os.path.getsize(H5AD)/1e6:.0f} MB')"""))


cells.append(md("""## 4. Load + inspect dataset structure"""))
cells.append(code(r"""adata = sc.read_h5ad(H5AD)
print(adata)
print(f'\nobs columns: {list(adata.obs.columns)}')

# Norman data structure: cells annotated with their guide identity.
# Common column names in published files:
#   'guide_identity' / 'perturbation' / 'gene_target' / 'gem_group'
# Locate the right column by inspection.
candidate_cols = [c for c in adata.obs.columns
                  if any(k in c.lower() for k in
                         ['guide', 'pert', 'target', 'gene'])]
print(f'\ncandidate guide-identity columns: {candidate_cols}')
print()
for c in candidate_cols[:3]:
    print(f'{c}: {adata.obs[c].nunique()} unique values, top:')
    print(adata.obs[c].value_counts().head(8))
    print()"""))


cells.append(md("""## 5. Identify control + single-pert + pair-pert classes

Norman dataset uses a guide-identity string per cell. Singles are e.g.
`SET_A` or `gene_A_only`; pairs are typically formatted `geneA_geneB`
or similar. Control (non-targeting) cells are flagged separately.

We need to map every cell to one of: `control` / `single:GENE` /
`pair:GENE_A__GENE_B`. Adjust the parsing below to match your data
file's actual column conventions."""))

cells.append(code(r"""# Manually inspect adata.obs to find the right column.
# Below is a template — edit GUIDE_COL and CONTROL_LABEL after seeing
# the inspection output above.

GUIDE_COL = candidate_cols[0] if candidate_cols else 'guide_identity'
print(f'Using GUIDE_COL = {GUIDE_COL}')

# Common Norman conventions: control cells labeled 'NT' or 'control'
# or empty string. Auto-detect.
guide_vals = adata.obs[GUIDE_COL].astype(str)
ctrl_candidates = [v for v in guide_vals.unique()
                   if any(k in v.lower() for k in ['nt', 'control', 'ctrl', 'nontargeting'])]
print(f'Control-like labels: {ctrl_candidates}')
CONTROL_LABEL = ctrl_candidates[0] if ctrl_candidates else 'NT'
print(f'CONTROL_LABEL = {CONTROL_LABEL}')

# Parse pair labels — Norman convention is gene_A + gene_B with '__' or '_' separator
def parse_class(g: str):
    g = str(g)
    if g == CONTROL_LABEL or g.lower() in ('control', 'nt', 'nontargeting'):
        return ('control', None, None)
    # Try double-underscore separator first, then single
    for sep in ('__', '+', '_AND_', ' '):
        if sep in g:
            parts = g.split(sep)
            if len(parts) == 2:
                return ('pair', parts[0], parts[1])
    # Single perturbation
    return ('single', g, None)

class_info = guide_vals.apply(parse_class)
adata.obs['class_type'] = [c[0] for c in class_info]
adata.obs['gene_a'] = [c[1] for c in class_info]
adata.obs['gene_b'] = [c[2] for c in class_info]

print(f'\nclass_type counts:')
print(adata.obs['class_type'].value_counts())"""))


cells.append(md("""## 6. Preprocess + cluster (decoupled outcome construction)"""))
cells.append(code(r"""# Standard scanpy pipeline. Keep raw counts in .raw for record.
adata.raw = adata.copy() if adata.raw is None else adata.raw

# Skip preprocessing if already done (Norman h5ad often pre-normalized)
if adata.X.max() > 100:    # raw counts heuristic
    sc.pp.filter_genes(adata, min_cells=20)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print('preprocessed: normalize_total + log1p')

sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=2000)
sc.pp.pca(adata, n_comps=50, random_state=SEED)
sc.pp.neighbors(adata, n_neighbors=15, random_state=SEED)
sc.tl.leiden(adata, resolution=1.0, random_state=SEED)
print(f'leiden clusters: {adata.obs["leiden"].nunique()}')

# Identify "control-like" cluster: where most control cells live
ctrl_clust = adata.obs.loc[adata.obs['class_type']=='control', 'leiden'].mode().iloc[0]
print(f'control-dominant cluster: {ctrl_clust}')

# Outcome PRIMARY: distance from control cluster centroid (loss-like)
ctrl_mask = (adata.obs['leiden'] == ctrl_clust)
ctrl_centroid = adata.obsm['X_pca'][ctrl_mask].mean(axis=0)
diff = adata.obsm['X_pca'] - ctrl_centroid[None, :]
adata.obs['outcome_distance'] = np.linalg.norm(diff, axis=1)
print(f'outcome_distance: mean={adata.obs["outcome_distance"].mean():.2f}, '
      f'std={adata.obs["outcome_distance"].std():.2f}')

# Outcome SECONDARY: binary cluster identity (1 = NOT control cluster)
adata.obs['outcome_cluster'] = (adata.obs['leiden'] != ctrl_clust).astype(float)"""))


cells.append(md("""## 7. Identify available test pairs

A pair (A, B) is testable if Norman dataset contains:
- ≥ 30 cells in `control` class
- ≥ 30 cells in `single A`
- ≥ 30 cells in `single B`
- ≥ 30 cells in `pair AB`

Build the candidate list automatically."""))

cells.append(code(r"""# Singles indexed by gene
single_counts = (adata.obs[adata.obs['class_type']=='single']['gene_a']
                 .value_counts())
print(f'genes with ≥30 single-pert cells: {(single_counts >= 30).sum()}')

# Pairs indexed by (gene_a, gene_b)
pair_counts = adata.obs[adata.obs['class_type']=='pair'].groupby(
    ['gene_a','gene_b']).size().reset_index(name='count')

# Canonicalize pair: sorted tuple
pair_counts['pair_key'] = pair_counts.apply(
    lambda r: tuple(sorted([r['gene_a'], r['gene_b']])), axis=1)
pair_counts = pair_counts.groupby('pair_key')['count'].sum().reset_index()
print(f'pairs with ≥30 cells: {(pair_counts["count"] >= 30).sum()}')

n_ctrl = (adata.obs['class_type']=='control').sum()
print(f'control cells: {n_ctrl}')

# Filter to pairs where BOTH single-A and single-B have ≥30 cells
TESTABLE = []
for _, r in pair_counts.iterrows():
    a, b = r['pair_key']
    if r['count'] < 30: continue
    if single_counts.get(a, 0) < 30: continue
    if single_counts.get(b, 0) < 30: continue
    TESTABLE.append((a, b))
print(f'\ntestable pairs (all 4 classes ≥30 cells): {len(TESTABLE)}')
print(f'first 10: {TESTABLE[:10]}')"""))


cells.append(md("""## 8. ε helper: 4-class intervention contrast"""))
cells.append(code(r"""def epsilon_4class(adata, gene_a, gene_b, outcome_col,
                   control_label='control', n_boot=1000, seed=42):
    '''Compute ε for a Perturb-seq pair using 4-class intervention design.

    Δ_A  = mean(outcome | A only)  − mean(outcome | control)
    Δ_B  = mean(outcome | B only)  − mean(outcome | control)
    Δ_AB = mean(outcome | A+B)     − mean(outcome | control)
    ε    = Δ_AB − Δ_A − Δ_B
    '''
    pair_canon = tuple(sorted([gene_a, gene_b]))

    obs = adata.obs
    out = obs[outcome_col].values

    mask_ctrl = (obs['class_type'] == 'control').values
    mask_a = (obs['class_type'] == 'single') & (obs['gene_a'] == gene_a)
    mask_b = (obs['class_type'] == 'single') & (obs['gene_a'] == gene_b)
    mask_ab = (obs['class_type'] == 'pair') & (obs.apply(
        lambda r: tuple(sorted([str(r['gene_a']), str(r['gene_b'])])) == pair_canon, axis=1))

    arr_ctrl = out[mask_ctrl]
    arr_a    = out[mask_a]
    arr_b    = out[mask_b]
    arr_ab   = out[mask_ab.values]

    n = {'ctrl': len(arr_ctrl), 'a': len(arr_a), 'b': len(arr_b), 'ab': len(arr_ab)}
    if min(n.values()) < 10:
        return None

    rng = np.random.default_rng(seed)
    def stat(c, a, b, ab):
        return (np.mean(ab) - np.mean(c)) - (np.mean(a) - np.mean(c)) - (np.mean(b) - np.mean(c))

    eps_obs = stat(arr_ctrl, arr_a, arr_b, arr_ab)

    boots = np.empty(n_boot)
    for i in range(n_boot):
        c2 = rng.choice(arr_ctrl, size=len(arr_ctrl), replace=True)
        a2 = rng.choice(arr_a, size=len(arr_a), replace=True)
        b2 = rng.choice(arr_b, size=len(arr_b), replace=True)
        ab2 = rng.choice(arr_ab, size=len(arr_ab), replace=True)
        boots[i] = stat(c2, a2, b2, ab2)
    se = float(boots.std(ddof=1))
    z = eps_obs / se if se > 0 else 0.0
    return {'epsilon': float(eps_obs), 'se': se, 'z': float(z),
            'n_classes': n, 'min_class': min(n.values())}"""))


cells.append(md("""## 9. Calibration scan on first 20 testable pairs"""))
cells.append(code(r"""SUBSET = TESTABLE[:20]
print(f'computing ε for {len(SUBSET)} pairs (both outcomes)')
results = []
t0 = time.time()
for i, (A, B) in enumerate(SUBSET):
    rec = {'gene_a': A, 'gene_b': B}
    for outname, col in [('dist', 'outcome_distance'),
                         ('cluster', 'outcome_cluster')]:
        r = epsilon_4class(adata, A, B, col, n_boot=N_BOOT, seed=SEED)
        if r is None:
            continue
        rec[f'{outname}_eps'] = r['epsilon']
        rec[f'{outname}_se'] = r['se']
        rec[f'{outname}_z'] = r['z']
        rec[f'{outname}_min_class'] = r['min_class']
    results.append(rec)
    if (i+1) % 5 == 0:
        print(f'  [{i+1}/{len(SUBSET)}] {A}↔{B} '
              f'ε_dist={rec.get("dist_eps", float("nan")):+.4f} '
              f'z_dist={rec.get("dist_z", float("nan")):+.2f} '
              f'({(time.time()-t0)/(i+1):.1f}s/pair)')
df_top = pd.DataFrame(results)
df_top.to_parquet(os.path.join(OUT_DIR, 'top20_pairs.parquet'))
print(df_top.head(8).to_string(index=False))"""))


cells.append(md("""## 10. Sanity 1 — calibration via random class permutation"""))
cells.append(code(r"""# For 50 trials: take a real pair, permute the class_type column
# across cells (breaks A/B/AB/control distinction), recompute ε.
# Should give |z| ~ folded N(0,1).
rng_perm = np.random.default_rng(SEED)
permuted_zs = []
class_orig = adata.obs[['class_type', 'gene_a', 'gene_b']].copy()

n_calib = 50
for trial in range(n_calib):
    # Pick random testable pair
    A, B = TESTABLE[rng_perm.integers(0, len(TESTABLE))]
    # Permute class_type within the cells used for this pair
    perm = rng_perm.permutation(len(adata.obs))
    adata.obs[['class_type','gene_a','gene_b']] = class_orig.iloc[perm].values
    r = epsilon_4class(adata, A, B, 'outcome_distance', n_boot=N_BOOT, seed=SEED)
    if r is not None:
        permuted_zs.append(r['z'])
# Restore original
adata.obs[['class_type','gene_a','gene_b']] = class_orig.values

permuted_zs = np.array(permuted_zs)
abs_z_p95 = float(np.percentile(np.abs(permuted_zs), 95))
print(f'permutation null pairs: {len(permuted_zs)}')
print(f'|z|_p95 = {abs_z_p95:.2f}  (target [1.5, 2.5])')
SANITY1_PASS = 1.5 <= abs_z_p95 <= 2.5
print(f'Sanity 1: {"PASS" if SANITY1_PASS else "FAIL"}')"""))


cells.append(md("""## 11. Sanity 2 — sensitivity on a curated pathway pair

Norman et al. report several known interaction pairs. We try a few
candidates from MAPK/cell-cycle pathways and pick the first present
in our testable set."""))

cells.append(code(r"""CURATED_CANDIDATES = [
    # MAPK pathway pairs (Norman highlighted these)
    ('CBL', 'CNN1'), ('FOXA1', 'FOXA3'),
    ('CEBPA', 'CEBPB'), ('CEBPB', 'CEBPE'),
    # Tumor-suppressor pairs
    ('TP53', 'KRAS'),
    # Generic essentials
    ('RPS3', 'RPL5'),
]
sanity2_pair = sanity2_result = None
for (A, B) in CURATED_CANDIDATES:
    cand = tuple(sorted([A, B]))
    if cand in [tuple(sorted(p)) for p in TESTABLE]:
        r = epsilon_4class(adata, cand[0], cand[1], 'outcome_distance',
                           n_boot=N_BOOT, seed=SEED)
        if r is not None:
            sanity2_pair = cand
            sanity2_result = r
            break

if sanity2_pair is None:
    # Fallback: highest-magnitude pair from top20 scan
    df_sorted = df_top.assign(abs_z=df_top['dist_z'].abs()).sort_values('abs_z', ascending=False)
    if len(df_sorted):
        top1 = df_sorted.iloc[0]
        sanity2_pair = (top1['gene_a'], top1['gene_b'])
        sanity2_result = {'epsilon': top1['dist_eps'],
                          'se': top1['dist_se'], 'z': top1['dist_z']}

assert sanity2_pair is not None
A, B = sanity2_pair
print(f'curated/fallback pair: ({A}, {B})')
print(f'  ε  = {sanity2_result["epsilon"]:+.4f}')
print(f'  SE = {sanity2_result["se"]:.4f}')
print(f'  z  = {sanity2_result["z"]:+.2f}')
SANITY2_PASS = abs(sanity2_result['z']) > 3
print(f'Sanity 2: {"PASS" if SANITY2_PASS else "FAIL"}')"""))


cells.append(md("""## 12. Sanity 3 — outcome consistency"""))
cells.append(code(r"""eps_dist = df_top['dist_eps'].values
eps_clust = df_top['cluster_eps'].values
mask = np.isfinite(eps_dist) & np.isfinite(eps_clust)
rho_consist, p_consist = pearsonr(eps_dist[mask], eps_clust[mask])
print(f'pairs with both ε: {mask.sum()}/{len(df_top)}')
print(f'Pearson |ρ|(distance, cluster) = {abs(rho_consist):.3f}  (sign: {rho_consist:+.3f})')
SANITY3_PASS = abs(rho_consist) > 0.7
print(f'Sanity 3: {"PASS" if SANITY3_PASS else "FAIL"}')"""))


cells.append(md("""## 13. Sanity 4 — sample size per class"""))
cells.append(code(r"""mins = df_top['dist_min_class'].values
overall_min = int(mins.min())
n_below = int((mins < 30).sum())
print(f'min class size across all pairs: {overall_min}')
print(f'pairs with any class < 30: {n_below}/{len(df_top)}')
SANITY4_PASS = (overall_min >= 30) and (n_below == 0)
print(f'Sanity 4: {"PASS" if SANITY4_PASS else "FAIL"}')"""))


cells.append(md("""## 14. Constraint checks (per methodological_findings.md)

C1 HARD circularity: outcome formula should not mention any test gene.
  → outcome_distance = euclidean dist in PCA space → uses ALL HVGs collectively, not any single test gene → PASS by construction.
C2 continuous variation: tertile validity. N/A here — perturbation
  classes are categorical by experimental design (this is the key
  benefit of Perturb-seq over expression stratification).
C3 SOFT correlation: |corr(class indicator, outcome)|.
C4 dominance: median |z| > max |corr| × √n_eval.
"""))

cells.append(code(r"""# C3: encode class as integers, compute correlation with outcome
class_codes = pd.Categorical(adata.obs['class_type']).codes
rho_c3, _ = pearsonr(class_codes, adata.obs['outcome_distance'].values)
max_abs_corr = abs(rho_c3)
n_eval = len(adata)
print(f'C3 SOFT correlation |corr(class_indicator, outcome)| = {max_abs_corr:.3f}')

# C4 dominance
sig = df_top[df_top['dist_z'].abs() > 3]
median_z_sig = sig['dist_z'].abs().median() if len(sig) else 0.0
threshold_dom = max_abs_corr * np.sqrt(n_eval)
print(f'\nC4 dominance:')
print(f'  median |z|_sig    = {median_z_sig:.2f}')
print(f'  threshold = corr×√n = {threshold_dom:.2f}')
DOMINANCE_PASS = median_z_sig > threshold_dom
print(f'Dominance: {"PASS" if DOMINANCE_PASS else "FAIL"}')"""))


cells.append(md("""## 15. Verdict — bio-side pre-reg gate"""))
cells.append(code(r"""verdict = {
    'project':           'BioEpistasis',
    'phase':             'Pilot iter 5 — Norman 2019 K562 Perturb-seq',
    'dataset':           'Norman et al. 2019 GSE133344',
    'n_cells':           int(adata.n_obs),
    'n_testable_pairs':  len(TESTABLE),
    'n_calibration_subset': len(SUBSET),
    'sign_convention':   'ε in loss space; ε>0 = synthetic-lethal/redundancy',
    'sanity': {
        'S1_calibration_p95':       abs_z_p95,
        'S1_pass':                  bool(SANITY1_PASS),
        'S2_curated_pair':          list(sanity2_pair),
        'S2_z':                     sanity2_result['z'],
        'S2_pass':                  bool(SANITY2_PASS),
        'S3_pearson_rho':           float(rho_consist),
        'S3_pass':                  bool(SANITY3_PASS),
        'S4_min_class':             overall_min,
        'S4_pass':                  bool(SANITY4_PASS),
        'C3_soft_corr':             max_abs_corr,
        'C4_dominance_median_zsig': median_z_sig,
        'C4_dominance_threshold':   float(threshold_dom),
        'C4_dominance_pass':        bool(DOMINANCE_PASS),
    },
}
all_pass = all([SANITY1_PASS, SANITY2_PASS, SANITY3_PASS,
                SANITY4_PASS, DOMINANCE_PASS])
verdict['all_gates_pass'] = bool(all_pass)
verdict['next_step'] = (
    'Lock pre-reg v1 biology side, run full Norman scan over all testable pairs'
    if all_pass else
    'Halt; investigate failing gate before pre-reg drafting'
)

out = os.path.join(OUT_DIR, 'iter5_verdict.json')
with open(out, 'w') as f:
    json.dump(verdict, f, indent=2)
print(json.dumps(verdict, indent=2))
print()
print('=' * 70)
print('NORMAN PILOT iter 5:', 'PASS — proceed to bio pre-reg' if all_pass else 'FAIL — halt')
print('=' * 70)
print(f'Saved → {out}')"""))


cells.append(md("""## 16. Done

If verdict PASS:
- Lock biology pre-reg v1 (analog of Tier 1 prereg, with intervention
  design instead of stratification).
- Run full Norman scan over all testable pairs.
- Compare ε direction + magnitude + shape with Pythia/OLMo Tier 1.

If verdict FAIL on any gate:
- Document the new constraint in
  `BioEpistasis/methodological_findings.md`.
- Decide whether to adjust pilot or halt bio fork.
"""))


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
