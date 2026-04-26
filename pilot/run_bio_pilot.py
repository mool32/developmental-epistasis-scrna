"""
Local execution of the BioEpistasis pilot — pivot to Bastidas-Ponce
pancreas (Schiebinger wot URL is now 404; secondary dataset becomes
the pilot calibration source. Schiebinger remains the production
pre-reg target if/when data resurfaces, OR Bastidas-Ponce is promoted
to primary).

Bastidas-Ponce E15.5 pancreas:
- ~3700 cells, multiple endocrine lineages (alpha, beta, delta, epsilon)
- Pre-processed (normalized + log-transformed) by scvelo authors
- Endpoint analog of "iPSC" = beta cells
- Marker score: insulin + Pdx1 + Mafa + Nkx6-1
"""

import os
import time
import json
from itertools import combinations

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


print("=" * 70)
print("SIGN CONVENTION (locked):")
print("  ε > 0 ↔ synthetic-lethal / redundancy phenotype")
print("  ε < 0 ↔ suppression / buffering")
print("=" * 70)

SEED = 42
N_BOOT = 1000
np.random.seed(SEED)

OUT_DIR = os.path.expanduser("~/BioEpistasis_pilot")
os.makedirs(OUT_DIR, exist_ok=True)
print(f"OUT_DIR = {OUT_DIR}\n")


# ── Load Bastidas-Ponce pancreas (built into scvelo) ─────────────────────────
print("=== Load Bastidas-Ponce pancreas ===")
adata = scv.datasets.pancreas()
print(adata)
print(f"\ncluster sizes:")
print(adata.obs["clusters"].value_counts())

# Standard scanpy preprocessing pipeline (scvelo's filter_and_normalize
# has API incompatibility in our scvelo/scanpy versions).
print("\n=== Preprocess: filter genes + normalize_total + log1p ===")
sc.pp.filter_genes(adata, min_cells=20)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print(adata)


# ── HVG selection ────────────────────────────────────────────────────────────
print("\n=== HVG selection ===")
# scvelo already flagged HVGs; rank by dispersion via cell_ranger flavor on
# the now-normalized .X. Fall back to scvelo's existing flag if HVG redo fails.
try:
    sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=200)
    hvg_df = adata.var[adata.var["highly_variable"]].copy()
    hvg_df = hvg_df.sort_values("dispersions_norm", ascending=False)
except Exception as e:
    print(f"cell_ranger HVG failed ({e}); using scvelo's flag")
    hvg_df = adata.var[adata.var["highly_variable"]].copy()
    if "dispersions_norm" in hvg_df.columns:
        hvg_df = hvg_df.sort_values("dispersions_norm", ascending=False)
    elif "means" in hvg_df.columns:
        hvg_df = hvg_df.sort_values("means", ascending=False)

TOP10 = hvg_df.head(10).index.tolist()
print(f"top-10 HVGs: {TOP10}")


# ── Outcome 1 (primary): beta-cell marker score ──────────────────────────────
print("\n=== Outcome: beta-cell marker score (primary) ===")
BETA_CANDIDATES = [
    ["Ins1", "Ins2", "Pdx1", "Mafa", "Nkx6-1"],
    ["INS1", "INS2", "PDX1", "MAFA", "NKX6-1"],
]
present = []
for guess in BETA_CANDIDATES:
    present = [g for g in guess if g in adata.var_names]
    if len(present) >= 3:
        break
print(f"beta markers present: {present}")
if len(present) < 3:
    print(f"sample of var_names: {list(adata.var_names[:20])}")
    print(f"genes containing 'ins': {[g for g in adata.var_names if 'ns' in g.lower()][:10]}")
    raise RuntimeError("too few beta markers; check var_names format")

sc.tl.score_genes(adata, gene_list=present, score_name="beta_score", random_state=SEED)

# Outcome PRIMARY (decoupled from any single gene's expression):
# cluster identity. 1 = not beta cell, 0 = beta cell. "Loss-like" sign.
# This avoids the circularity that plagued marker-score outcome
# (when test genes are themselves in the marker set).
adata.obs["outcome_marker"] = (adata.obs["clusters"] != "Beta").astype(float)
print(f"outcome_marker (1 = not-beta): mean={adata.obs['outcome_marker'].mean():.3f}")

# ── Outcome 2 (secondary): PCA distance to beta centroid ─────────────────────
print("\n=== Outcome: PCA distance to beta cluster centroid (secondary) ===")
sc.pp.pca(adata, n_comps=50, random_state=SEED)
beta_mask = (adata.obs["clusters"] == "Beta").values
beta_centroid = adata.obsm["X_pca"][beta_mask].mean(axis=0)
diff = adata.obsm["X_pca"] - beta_centroid[None, :]
adata.obs["outcome_distance"] = np.linalg.norm(diff, axis=1)
print(f"outcome_distance mean={adata.obs['outcome_distance'].mean():.3f}")
print(f"beta cluster cells: {beta_mask.sum()}")


# ── Helpers ──────────────────────────────────────────────────────────────────
def stratify_2x2(adata, gene_a, gene_b, outcome_col, q_low=1/3, q_high=2/3):
    if gene_a not in adata.var_names or gene_b not in adata.var_names:
        return None
    expr_a = adata[:, gene_a].X
    expr_b = adata[:, gene_b].X
    if hasattr(expr_a, "toarray"):
        expr_a = expr_a.toarray()
    if hasattr(expr_b, "toarray"):
        expr_b = expr_b.toarray()
    expr_a = np.asarray(expr_a).flatten()
    expr_b = np.asarray(expr_b).flatten()

    a_low_thr = np.quantile(expr_a, q_low)
    a_high_thr = np.quantile(expr_a, q_high)
    b_low_thr = np.quantile(expr_b, q_low)
    b_high_thr = np.quantile(expr_b, q_high)

    out = adata.obs[outcome_col].values
    a_high = expr_a >= a_high_thr
    a_low = expr_a <= a_low_thr
    b_high = expr_b >= b_high_thr
    b_low = expr_b <= b_low_thr

    return {
        "HH": out[a_high & b_high],
        "HL": out[a_high & b_low],
        "LH": out[a_low & b_high],
        "LL": out[a_low & b_low],
    }


def bootstrap_epsilon(corners, n_boot=1000, seed=42):
    n = {k: len(v) for k, v in corners.items()}
    if min(n.values()) < 10:
        return None
    rng = np.random.default_rng(seed)

    # MEAN (not median) per corner. Median is unstable on sparse scRNA
    # data because bootstrap-resampled medians flip discretely between
    # mass points (e.g., 0 vs first nonzero), inflating SE artificially.
    def stat(c):
        return (
            (np.mean(c["LL"]) - np.mean(c["HH"]))
            - (np.mean(c["LH"]) - np.mean(c["HH"]))
            - (np.mean(c["HL"]) - np.mean(c["HH"]))
        )

    eps_obs = stat(corners)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        resampled = {k: rng.choice(v, size=len(v), replace=True) for k, v in corners.items()}
        boots[i] = stat(resampled)
    se = float(boots.std(ddof=1))
    z = eps_obs / se if se > 0 else 0.0
    return {"epsilon": float(eps_obs), "se": se, "z": float(z), "n_corners": n}


# ── Top-10 pairs ─────────────────────────────────────────────────────────────
print("\n=== Top-10 HVG pairs ===")
PAIRS_TOP10 = list(combinations(TOP10, 2))
print(f"top-10 pair count: {len(PAIRS_TOP10)}")

results_top10 = []
for (A, B) in PAIRS_TOP10:
    rec = {"gene_a": A, "gene_b": B}
    for outcome_name, col in [("marker", "outcome_marker"), ("distance", "outcome_distance")]:
        corners = stratify_2x2(adata, A, B, col)
        if corners is None:
            continue
        boot = bootstrap_epsilon(corners, n_boot=N_BOOT, seed=SEED)
        if boot is None:
            continue
        rec[f"{outcome_name}_eps"] = boot["epsilon"]
        rec[f"{outcome_name}_se"] = boot["se"]
        rec[f"{outcome_name}_z"] = boot["z"]
        rec[f"{outcome_name}_nHH"] = boot["n_corners"]["HH"]
        rec[f"{outcome_name}_nHL"] = boot["n_corners"]["HL"]
        rec[f"{outcome_name}_nLH"] = boot["n_corners"]["LH"]
        rec[f"{outcome_name}_nLL"] = boot["n_corners"]["LL"]
    results_top10.append(rec)
df_top = pd.DataFrame(results_top10)
df_top.to_parquet(os.path.join(OUT_DIR, "top10_pairs.parquet"))
print(df_top[["gene_a", "gene_b", "marker_eps", "marker_z",
              "distance_eps", "distance_z"]].to_string(index=False))


# ── Sanity 1: bootstrap calibration ──────────────────────────────────────────
print("\n=== Sanity 1 — bootstrap calibration ===")
rng_random = np.random.default_rng(SEED)
all_genes = adata.var_names.tolist()
n_random = 50
random_pairs, seen = [], set()
while len(random_pairs) < n_random:
    i, j = rng_random.integers(0, len(all_genes), size=2)
    if i == j:
        continue
    a, b = all_genes[int(i)], all_genes[int(j)]
    if a == b or (a, b) in seen or (b, a) in seen:
        continue
    seen.add((a, b))
    random_pairs.append((a, b))

random_zs = []
for (A, B) in random_pairs:
    corners = stratify_2x2(adata, A, B, "outcome_marker")
    if corners is None:
        continue
    boot = bootstrap_epsilon(corners, n_boot=N_BOOT, seed=SEED)
    if boot is None:
        continue
    random_zs.append(boot["z"])
random_zs = np.array(random_zs)
abs_z_p95 = float(np.percentile(np.abs(random_zs), 95))
print(f"random pairs measured: {len(random_zs)}")
print(f"|z| 95th percentile  : {abs_z_p95:.2f}  (target [1.5, 2.5])")
SANITY1_PASS = 1.5 <= abs_z_p95 <= 2.5
print(f"Sanity 1: {'PASS' if SANITY1_PASS else 'FAIL'}")


# ── Sanity 2: curated same-pathway pair ──────────────────────────────────────
print("\n=== Sanity 2 — curated same-pathway pair ===")
CURATED_PAIRS = [
    ("Ins1", "Ins2"),    # both insulin paralogs — strong same-pathway
    ("Rps3", "Rps5"),
    ("Rpl3", "Rpl5"),
    ("Eef1a1", "Eef2"),
    ("Actb", "Tubb5"),
]
sanity2_pair = sanity2_result = None
for (A, B) in CURATED_PAIRS:
    if A in adata.var_names and B in adata.var_names:
        corners = stratify_2x2(adata, A, B, "outcome_marker")
        if corners is None:
            continue
        boot = bootstrap_epsilon(corners, n_boot=N_BOOT, seed=SEED)
        if boot is None:
            continue
        sanity2_pair = (A, B)
        sanity2_result = boot
        break

if sanity2_pair is None:
    rib_or_actin = [g for g in adata.var_names if g.lower().startswith(("rps", "rpl", "act"))][:6]
    print(f"fallback candidates: {rib_or_actin}")
    for i, j in combinations(range(len(rib_or_actin)), 2):
        A, B = rib_or_actin[i], rib_or_actin[j]
        corners = stratify_2x2(adata, A, B, "outcome_marker")
        if corners is None:
            continue
        boot = bootstrap_epsilon(corners, n_boot=N_BOOT, seed=SEED)
        if boot is None:
            continue
        sanity2_pair = (A, B); sanity2_result = boot; break

assert sanity2_pair is not None, "no usable sanity-2 pair"
A, B = sanity2_pair
print(f"curated pair: ({A}, {B})")
print(f"  ε  = {sanity2_result['epsilon']:+.5f}")
print(f"  SE = {sanity2_result['se']:.5f}")
print(f"  z  = {sanity2_result['z']:+.2f}")
SANITY2_PASS = abs(sanity2_result["z"]) > 3
print(f"Sanity 2: {'PASS' if SANITY2_PASS else 'FAIL'}")


# ── Sanity 3: outcome consistency ────────────────────────────────────────────
print("\n=== Sanity 3 — outcome consistency ===")
eps_marker = df_top["marker_eps"].values
eps_dist = df_top["distance_eps"].values
mask = np.isfinite(eps_marker) & np.isfinite(eps_dist)
rho_consist, p_consist = pearsonr(eps_marker[mask], eps_dist[mask])
print(f"pairs with both ε: {mask.sum()}/{len(df_top)}")
print(f"Pearson ρ(marker, distance): {rho_consist:+.3f}  (p={p_consist:.4f})")
SANITY3_PASS = rho_consist > 0.7
print(f"Sanity 3: {'PASS' if SANITY3_PASS else 'FAIL'}")

fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(eps_marker[mask], eps_dist[mask], alpha=0.6)
ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
ax.set_xlabel("ε (marker outcome)")
ax.set_ylabel("ε (PCA distance outcome)")
ax.set_title(f"Outcome consistency on 45 top-10 pairs (ρ={rho_consist:.3f})")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "sanity3_outcome_consistency.png"), dpi=120)


# ── Sanity 4: sample size per quadrant ───────────────────────────────────────
print("\n=== Sanity 4 — sample size per quadrant ===")
corner_cols = [c for c in df_top.columns if c.startswith("marker_n")]
corner_mins = df_top[corner_cols].min(axis=1)
overall_min = int(corner_mins.min())
n_below_30 = int((corner_mins < 30).sum())
print(f"min cells/corner: {overall_min}")
print(f"pairs with corner < 30: {n_below_30}/{len(df_top)}")
SANITY4_PASS = (overall_min >= 30) and (n_below_30 == 0)
print(f"Sanity 4: {'PASS' if SANITY4_PASS else 'FAIL'}")


# ── Verdict ──────────────────────────────────────────────────────────────────
print("\n=== Verdict ===")
sanity_results = {
    "sanity1_calibration_p95": abs_z_p95,
    "sanity1_pass": bool(SANITY1_PASS),
    "sanity2_curated_pair": list(sanity2_pair),
    "sanity2_z": sanity2_result["z"],
    "sanity2_pass": bool(SANITY2_PASS),
    "sanity3_pearson_rho": float(rho_consist),
    "sanity3_pass": bool(SANITY3_PASS),
    "sanity4_min_corner_cells": overall_min,
    "sanity4_pass": bool(SANITY4_PASS),
}
all_pass = all([SANITY1_PASS, SANITY2_PASS, SANITY3_PASS, SANITY4_PASS])

verdict = {
    "project": "BioEpistasis",
    "phase": "Pilot — calibration on Bastidas-Ponce E15.5 pancreas (Schiebinger URL 404, pivot)",
    "dataset": "Bastidas-Ponce 2019 / scvelo built-in pancreas",
    "n_cells": int(adata.n_obs),
    "top10_hvg": TOP10,
    "n_top10_pairs": len(PAIRS_TOP10),
    "n_random_pairs": n_random,
    "outcome_primary": "beta-cell marker score",
    "outcome_secondary": "PCA distance to beta-marker-high centroid",
    "sign_convention": "ε in loss space; ε>0 = synthetic-lethal/redundancy",
    "sanity": sanity_results,
    "all_sanity_pass": bool(all_pass),
    "next_step": "Draft pre-reg v1 if PASS; address sanity if FAIL.",
}
out = os.path.join(OUT_DIR, "pilot_calibration_verdict.json")
with open(out, "w") as f:
    json.dump(verdict, f, indent=2)
print(json.dumps(verdict, indent=2))
print()
print("=" * 70)
print("PILOT VERDICT:", "PASS — proceed to pre-reg" if all_pass else "FAIL — halt")
print("=" * 70)
print(f"\nSaved → {out}")
