"""
BioEpistasis pilot Iteration 4 — synthetic null test on paul15.

Final gate before pre-reg lock. Tests whether the pipeline produces a
correct null distribution under "no real epistasis" — i.e., when one
gene's expression is permuted across cells (breaks any A-B
relationship while preserving marginal distributions).

Three pass criteria (all must PASS):
  N1. |mean(ε_null)| < 0.1 × |mean(ε_observed_top10)|
  N2. SE_null consistent with bootstrap SE estimate (median ratio
      within [0.5, 2.0])
  N3. |z_null| distribution close to folded N(0,1):
      one-sample KS test against half-normal CDF, p > 0.05.

Loads paul15 + reuses all preprocessing + outcome from Iter 3.
"""

import os
import json
from itertools import combinations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, kstest, halfnorm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


print("=" * 70)
print("ITER 4 — synthetic null test (final pre-reg gate)")
print("=" * 70)

SEED = 42
N_BOOT = 1000
N_NULL = 100      # synthetic null pairs
np.random.seed(SEED)

OUT_DIR = os.path.expanduser("~/BioEpistasis_pilot/iter4_synth_null")
os.makedirs(OUT_DIR, exist_ok=True)


# ── Load + preprocess (identical to Iter 3) ──────────────────────────────────
print("\n=== Load paul15 + preprocess ===")
adata = sc.datasets.paul15()
sc.pp.recipe_zheng17(adata, n_top_genes=200)

sc.tl.pca(adata, random_state=SEED)
sc.pp.neighbors(adata, n_neighbors=15, random_state=SEED)
sc.tl.diffmap(adata, random_state=SEED)
adata.uns["iroot"] = int(np.flatnonzero(adata.obs["paul15_clusters"] == "7MEP")[0])
sc.tl.dpt(adata)
adata.obs["outcome_pseudotime"] = adata.obs["dpt_pseudotime"].astype(float)

hvg_df = adata.var.copy().sort_values("std", ascending=False)
TOP30 = hvg_df.head(30).index.tolist()
print(f"top-30 HVG: {TOP30[:6]}...")
print(f"n_cells = {adata.n_obs}")


def get_expr(adata, g):
    x = adata[:, g].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).flatten()


def bootstrap_epsilon_2x2(expr_a, expr_b, outcome, n_boot=1000, seed=42):
    """Compute ε via 2×2 tertile stratification + bootstrap SE."""
    a_low_thr = np.quantile(expr_a, 1/3)
    a_high_thr = np.quantile(expr_a, 2/3)
    b_low_thr = np.quantile(expr_b, 1/3)
    b_high_thr = np.quantile(expr_b, 2/3)

    a_high = expr_a >= a_high_thr
    a_low = expr_a <= a_low_thr
    b_high = expr_b >= b_high_thr
    b_low = expr_b <= b_low_thr
    corners = {
        "HH": outcome[a_high & b_high],
        "HL": outcome[a_high & b_low],
        "LH": outcome[a_low & b_high],
        "LL": outcome[a_low & b_low],
    }
    n = {k: len(v) for k, v in corners.items()}
    if min(n.values()) < 10:
        return None

    rng = np.random.default_rng(seed)
    def stat(c):
        return ((np.mean(c["LL"]) - np.mean(c["HH"]))
                - (np.mean(c["LH"]) - np.mean(c["HH"]))
                - (np.mean(c["HL"]) - np.mean(c["HH"])))
    eps_obs = stat(corners)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        resampled = {k: rng.choice(v, size=len(v), replace=True) for k, v in corners.items()}
        boots[i] = stat(resampled)
    se = float(boots.std(ddof=1))
    z = eps_obs / se if se > 0 else 0.0
    return {"epsilon": float(eps_obs), "se": se, "z": float(z),
            "n_corners": n, "min_corner": min(n.values())}


# ── Compute observed ε for top-10 (= reference for N1 dominance test) ───────
print("\n=== Compute observed ε on top-10 (filtered, min_corner ≥ 20) ===")
TOP10 = TOP30[:10]
outcome = adata.obs["outcome_pseudotime"].values
observed = []
for (A, B) in combinations(TOP10, 2):
    eA = get_expr(adata, A)
    eB = get_expr(adata, B)
    r = bootstrap_epsilon_2x2(eA, eB, outcome, n_boot=N_BOOT, seed=SEED)
    if r is None:
        continue
    if r["min_corner"] < 20:        # S4 refined gate
        continue
    r.update({"gene_a": A, "gene_b": B})
    observed.append(r)
df_obs = pd.DataFrame(observed)
print(f"observed pairs surviving min_corner ≥ 20: {len(df_obs)}")
mean_eps_observed = float(df_obs["epsilon"].mean())
median_z_observed = float(df_obs["z"].abs().median())
sig_observed = df_obs[df_obs["z"].abs() > 3]
n_sig_observed = len(sig_observed)
median_z_sig = float(sig_observed["z"].abs().median()) if n_sig_observed else 0
print(f"mean ε observed       = {mean_eps_observed:+.5f}")
print(f"median |z| observed   = {median_z_observed:.2f}")
print(f"|z|>3 significant     = {n_sig_observed}/{len(df_obs)}")
print(f"median |z| of significant: {median_z_sig:.2f}")


# ── Generate synthetic null pairs ────────────────────────────────────────────
print("\n=== Generate synthetic null: permute one gene per pair ===")
rng = np.random.default_rng(SEED + 1)
all_genes = adata.var_names.tolist()
n_cells = adata.n_obs

null_records = []
attempts = 0
while len(null_records) < N_NULL and attempts < N_NULL * 5:
    attempts += 1
    # Pick A from top-30, B random gene (broader pool for variety)
    i = rng.integers(0, len(TOP30))
    j = rng.integers(0, len(all_genes))
    A = TOP30[int(i)]
    B = all_genes[int(j)]
    if A == B:
        continue
    eA = get_expr(adata, A)
    eB = get_expr(adata, B)
    # Permute B's expression across cells — kills A-B correlation,
    # preserves B's marginal distribution
    eB_shuf = rng.permutation(eB)
    r = bootstrap_epsilon_2x2(eA, eB_shuf, outcome, n_boot=N_BOOT, seed=SEED + len(null_records))
    if r is None:
        continue
    if r["min_corner"] < 20:
        continue
    r.update({"gene_a": A, "gene_b_shuffled": B})
    null_records.append(r)

df_null = pd.DataFrame(null_records)
print(f"synthetic null pairs: {len(df_null)} (target {N_NULL})")
print(f"  mean ε_null     = {df_null['epsilon'].mean():+.5f}")
print(f"  std ε_null      = {df_null['epsilon'].std():.5f}")
print(f"  median SE_null  = {df_null['se'].median():.5f}")
print(f"  median |z|_null = {df_null['z'].abs().median():.2f}")
print(f"  |z|_null > 3    = {(df_null['z'].abs() > 3).sum()}/{len(df_null)}")


# ── N1: mean ε_null close to 0 ───────────────────────────────────────────────
print("\n=== N1: |mean(ε_null)| < 0.1 × |mean(ε_observed)| ===")
mean_null = float(df_null["epsilon"].mean())
threshold_N1 = 0.1 * abs(mean_eps_observed)
print(f"|mean(ε_null)|    = {abs(mean_null):.5f}")
print(f"threshold (0.1×)  = {threshold_N1:.5f}")
N1_PASS = abs(mean_null) < threshold_N1
print(f"N1: {'PASS' if N1_PASS else 'FAIL'}")


# ── N2: SE_null consistent with observed SE ──────────────────────────────────
print("\n=== N2: SE_null consistent with bootstrap SE estimate ===")
median_se_null = float(df_null["se"].median())
median_se_observed = float(df_obs["se"].median())
ratio = median_se_null / median_se_observed if median_se_observed > 0 else float("inf")
print(f"median SE_null      = {median_se_null:.5f}")
print(f"median SE_observed  = {median_se_observed:.5f}")
print(f"ratio (null/obs)    = {ratio:.2f}  (target [0.5, 2.0])")
N2_PASS = 0.5 <= ratio <= 2.0
print(f"N2: {'PASS' if N2_PASS else 'FAIL'}")


# ── N3: |z_null| ~ folded N(0,1) via KS test ─────────────────────────────────
print("\n=== N3: |z_null| ~ folded N(0,1) via KS test ===")
abs_z_null = np.abs(df_null["z"].values)
ks_stat, ks_p = kstest(abs_z_null, halfnorm.cdf)
print(f"KS statistic D = {ks_stat:.3f}")
print(f"KS p-value     = {ks_p:.4f}  (target > 0.05)")
N3_PASS = ks_p > 0.05
print(f"N3: {'PASS' if N3_PASS else 'FAIL'}")


# ── Constraint-3 dominance criterion (per user spec) ─────────────────────────
print("\n=== Constraint 3 — dominance: median |z| > max |corr| × sqrt(n) ===")
# max |corr(outcome, gene_in_top10)| from Iter 3 was 0.476.
# Recompute here for closure.
max_abs_corr = 0.0
for g in TOP10:
    e = get_expr(adata, g)
    rho, _ = pearsonr(e, outcome)
    if abs(rho) > max_abs_corr:
        max_abs_corr = abs(rho)
n_eval = adata.n_obs
threshold_dom = max_abs_corr * np.sqrt(n_eval)
print(f"max |corr|        = {max_abs_corr:.3f}")
print(f"sqrt(n_eval)      = {np.sqrt(n_eval):.1f}")
print(f"threshold = corr×sqrt(n) = {threshold_dom:.2f}")
print(f"median |z| significant   = {median_z_sig:.2f}")
DOMINANCE_PASS = median_z_sig > threshold_dom
print(f"Dominance: {'PASS' if DOMINANCE_PASS else 'FAIL'}")


# ── Verdict ──────────────────────────────────────────────────────────────────
print("\n=== Synthetic null verdict ===")
all_pass = bool(N1_PASS and N2_PASS and N3_PASS and DOMINANCE_PASS)

verdict = {
    "project":         "BioEpistasis",
    "phase":           "Pilot Iter 4 — synthetic null test (final pre-reg gate)",
    "dataset":         "scanpy.datasets.paul15 (myeloid lineage)",
    "n_cells":         int(adata.n_obs),
    "n_top10_pairs_surviving_S4": len(df_obs),
    "n_null_pairs":    len(df_null),
    "observed": {
        "mean_epsilon":          mean_eps_observed,
        "median_abs_z":          median_z_observed,
        "n_significant_z_gt_3":  int(n_sig_observed),
        "median_abs_z_significant": median_z_sig,
    },
    "synthetic_null": {
        "mean_epsilon":     float(mean_null),
        "std_epsilon":      float(df_null['epsilon'].std()),
        "median_se":        median_se_null,
        "median_abs_z":     float(df_null['z'].abs().median()),
        "frac_abs_z_gt_3":  float((df_null['z'].abs() > 3).mean()),
    },
    "N1_mean_null_close_zero": bool(N1_PASS),
    "N2_se_consistent":        bool(N2_PASS),
    "N3_z_distribution_KS":    bool(N3_PASS),
    "N3_KS_stat":              float(ks_stat),
    "N3_KS_p":                 float(ks_p),
    "constraint3_dominance":   bool(DOMINANCE_PASS),
    "constraint3_max_corr":    max_abs_corr,
    "constraint3_threshold":   float(threshold_dom),
    "all_synthetic_null_pass": all_pass,
    "next_step": "Pre-reg v1 lock under all refined + synthetic-null gates"
                 if all_pass else
                 "Halt before pre-reg; investigate failing null gate",
}

# Save
df_obs.to_parquet(os.path.join(OUT_DIR, "observed_pairs.parquet"))
df_null.to_parquet(os.path.join(OUT_DIR, "synthetic_null_pairs.parquet"))
out = os.path.join(OUT_DIR, "iter4_synth_null_verdict.json")
with open(out, "w") as f:
    json.dump(verdict, f, indent=2)
print(json.dumps(verdict, indent=2))


# ── Plots ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# (a) ε distributions: observed vs null
ax = axes[0]
ax.hist(df_null["epsilon"], bins=20, alpha=0.6, label=f"null (n={len(df_null)})", color="C0")
ax.hist(df_obs["epsilon"], bins=15, alpha=0.6, label=f"observed (n={len(df_obs)})", color="C3")
ax.axvline(0, color="gray", lw=0.5)
ax.set_xlabel("ε"); ax.set_ylabel("count")
ax.set_title(f"ε: observed vs synthetic null")
ax.legend()

# (b) |z| distributions
ax = axes[1]
ax.hist(np.abs(df_null["z"]), bins=20, alpha=0.6, label="null", color="C0", density=True)
ax.hist(np.abs(df_obs["z"]), bins=15, alpha=0.6, label="observed", color="C3", density=True)
xx = np.linspace(0, max(np.abs(df_null["z"]).max(), 5), 200)
ax.plot(xx, halfnorm.pdf(xx), "k--", lw=1.5, label="folded N(0,1)")
ax.set_xlabel("|z|"); ax.set_ylabel("density")
ax.set_title(f"|z|: null vs N(0,1), KS p={ks_p:.3f}")
ax.legend()

# (c) SE: null vs observed
ax = axes[2]
ax.hist(df_null["se"], bins=20, alpha=0.6, label="null SE", color="C0")
ax.hist(df_obs["se"], bins=15, alpha=0.6, label="observed SE", color="C3")
ax.set_xlabel("bootstrap SE(ε)"); ax.set_ylabel("count")
ax.set_title(f"SE distribution (ratio null/obs = {ratio:.2f})")
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "iter4_synth_null_diagnostic.png"), dpi=120)

print()
print("=" * 70)
print("ITER 4 SYNTHETIC NULL VERDICT:",
      "PASS — proceed to pre-reg v1 lock" if all_pass else "FAIL — investigate")
print("=" * 70)
print(f"Saved → {out}")
