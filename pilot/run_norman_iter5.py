"""
Local execution of Norman 2019 pilot iter 5.
Mirrors developmental-epistasis-scrna/pilot/05_norman_iter5.ipynb but
runs as a script.

Compute: ~30-60 min CPU. Download: ~700 MB on first run.
"""

import os, json, time, subprocess
from itertools import combinations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


print("=" * 70)
print("SIGN CONVENTION (locked):")
print("  ε > 0 ↔ synthetic-lethal / redundancy phenotype")
print("  ε < 0 ↔ suppression / buffering")
print("=" * 70)
print()

SEED = 42
N_BOOT = 1000
np.random.seed(SEED)

OUT_DIR = os.path.expanduser("~/BioEpistasis_pilot/iter5_norman")
DATA_DIR = os.path.join(OUT_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)
print(f"OUT_DIR: {OUT_DIR}")


# ── Download Norman 2019 from scPerturb Zenodo ───────────────────────────────
URL = "https://zenodo.org/records/10044268/files/NormanWeissman2019_filtered.h5ad?download=1"
H5AD = os.path.join(DATA_DIR, "NormanWeissman2019_filtered.h5ad")

if not os.path.exists(H5AD) or os.path.getsize(H5AD) < 100e6:
    print(f"\nDownloading {URL}")
    print("Expected ~698 MB; takes 2-10 min on residential.")
    t0 = time.time()
    subprocess.check_call([
        "curl", "-L", "--retry", "3", "--max-time", "600",
        "-o", H5AD, URL,
    ])
    sz = os.path.getsize(H5AD)
    print(f"Downloaded: {sz/1e6:.0f} MB in {time.time()-t0:.0f}s")
    assert sz > 100e6, f"Download incomplete: {sz} bytes"
else:
    print(f"Cached: {os.path.getsize(H5AD)/1e6:.0f} MB")


# ── Load + inspect ────────────────────────────────────────────────────────────
print("\n=== Load AnnData ===")
adata = sc.read_h5ad(H5AD)
print(adata)
print(f"\nobs columns: {list(adata.obs.columns)}")

candidate_cols = [c for c in adata.obs.columns
                  if any(k in c.lower() for k in
                         ["guide", "pert", "target", "gene"])]
print(f"\nguide-id candidate columns: {candidate_cols}")
for c in candidate_cols[:3]:
    n_unique = adata.obs[c].nunique()
    print(f"\n{c}: {n_unique} unique")
    print(adata.obs[c].value_counts().head(8))


# ── Norman dataset uses scPerturb harmonized format ──────────────────────────
# `perturbation` column has clean labels (control / single / pair)
# `nperts` column has 0/1/2 (number of real perturbations per cell)
GUIDE_COL = "perturbation"
CONTROL_LABEL = "control"
print(f"\nUsing GUIDE_COL = {GUIDE_COL}")
print(f"CONTROL_LABEL = {CONTROL_LABEL}")
print(f"\n`nperts` distribution:")
print(adata.obs["nperts"].value_counts())


def parse_row(g: str, n: int):
    if g == CONTROL_LABEL or n == 0:
        return ("control", None, None)
    if n == 1:
        return ("single", g, None)
    if n == 2:
        # Pair label is GENE_A_GENE_B (Norman convention).
        parts = g.split("_")
        if len(parts) == 2:
            return ("pair", parts[0], parts[1])
        # Edge case (rarely): gene name with underscore — treat first as A, rest as B.
        return ("pair", parts[0], "_".join(parts[1:]))
    return ("unknown", None, None)


class_info = [parse_row(str(g), int(n))
              for g, n in zip(adata.obs[GUIDE_COL], adata.obs["nperts"])]
adata.obs["class_type"] = [c[0] for c in class_info]
adata.obs["gene_a"] = [c[1] for c in class_info]
adata.obs["gene_b"] = [c[2] for c in class_info]

print(f"\nclass_type counts:")
print(adata.obs["class_type"].value_counts())


# ── Preprocess if needed ─────────────────────────────────────────────────────
print("\n=== Preprocess ===")
if adata.X.max() > 100:
    print("Raw counts detected, normalizing")
    sc.pp.filter_genes(adata, min_cells=20)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
else:
    print("Already log-normalized")

sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=2000)
sc.pp.pca(adata, n_comps=50, random_state=SEED)
sc.pp.neighbors(adata, n_neighbors=15, random_state=SEED)
sc.tl.leiden(adata, resolution=1.0, random_state=SEED)
print(f"leiden clusters: {adata.obs['leiden'].nunique()}")

ctrl_mask = (adata.obs["class_type"] == "control")
ctrl_clust = adata.obs.loc[ctrl_mask, "leiden"].mode().iloc[0]
print(f"control-dominant cluster: {ctrl_clust}")

# Outcomes
ctrl_clust_mask = (adata.obs["leiden"] == ctrl_clust)
ctrl_centroid = adata.obsm["X_pca"][ctrl_clust_mask].mean(axis=0)
diff = adata.obsm["X_pca"] - ctrl_centroid[None, :]
adata.obs["outcome_distance"] = np.linalg.norm(diff, axis=1)
adata.obs["outcome_cluster"] = (adata.obs["leiden"] != ctrl_clust).astype(float)
print(f"outcome_distance mean={adata.obs['outcome_distance'].mean():.2f}")
print(f"outcome_cluster mean={adata.obs['outcome_cluster'].mean():.3f}")


# ── Identify testable pairs ──────────────────────────────────────────────────
print("\n=== Identify testable pairs ===")
single_counts = (adata.obs[adata.obs["class_type"] == "single"]["gene_a"]
                 .value_counts())
print(f"genes with ≥30 single-pert cells: {(single_counts >= 30).sum()}")

pair_counts = (adata.obs[adata.obs["class_type"] == "pair"]
               .groupby(["gene_a", "gene_b"]).size().reset_index(name="count"))
pair_counts["pair_key"] = pair_counts.apply(
    lambda r: tuple(sorted([r["gene_a"], r["gene_b"]])), axis=1)
pair_counts = pair_counts.groupby("pair_key")["count"].sum().reset_index()
print(f"pairs with ≥30 cells: {(pair_counts['count'] >= 30).sum()}")

n_ctrl = (adata.obs["class_type"] == "control").sum()
print(f"control cells: {n_ctrl}")

TESTABLE = []
for _, r in pair_counts.iterrows():
    a, b = r["pair_key"]
    if r["count"] < 30:
        continue
    if single_counts.get(a, 0) < 30:
        continue
    if single_counts.get(b, 0) < 30:
        continue
    TESTABLE.append((a, b))
print(f"\ntestable pairs (4 classes ≥30 cells): {len(TESTABLE)}")
print(f"first 10: {TESTABLE[:10]}")
assert len(TESTABLE) >= 5, "fewer than 5 testable pairs — design unviable"


# ── 4-class ε helper ─────────────────────────────────────────────────────────
def epsilon_4class(adata, gene_a, gene_b, outcome_col, n_boot=1000, seed=42):
    pair_canon = tuple(sorted([gene_a, gene_b]))
    obs = adata.obs
    out = obs[outcome_col].values

    mask_ctrl = (obs["class_type"] == "control").values
    mask_a = ((obs["class_type"] == "single") & (obs["gene_a"] == gene_a)).values
    mask_b = ((obs["class_type"] == "single") & (obs["gene_a"] == gene_b)).values
    mask_ab = ((obs["class_type"] == "pair") & obs.apply(
        lambda r: tuple(sorted([str(r["gene_a"]), str(r["gene_b"])])) == pair_canon,
        axis=1)).values

    arr_ctrl, arr_a, arr_b, arr_ab = (
        out[mask_ctrl], out[mask_a], out[mask_b], out[mask_ab])
    n = {"ctrl": len(arr_ctrl), "a": len(arr_a),
         "b": len(arr_b), "ab": len(arr_ab)}
    if min(n.values()) < 10:
        return None

    rng = np.random.default_rng(seed)
    def stat(c, a, b, ab):
        return ((np.mean(ab) - np.mean(c))
                - (np.mean(a) - np.mean(c))
                - (np.mean(b) - np.mean(c)))
    eps_obs = stat(arr_ctrl, arr_a, arr_b, arr_ab)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        c2 = rng.choice(arr_ctrl, len(arr_ctrl), replace=True)
        a2 = rng.choice(arr_a, len(arr_a), replace=True)
        b2 = rng.choice(arr_b, len(arr_b), replace=True)
        ab2 = rng.choice(arr_ab, len(arr_ab), replace=True)
        boots[i] = stat(c2, a2, b2, ab2)
    se = float(boots.std(ddof=1))
    z = eps_obs / se if se > 0 else 0.0
    return {"epsilon": float(eps_obs), "se": se, "z": float(z),
            "n_classes": n, "min_class": min(n.values())}


# ── Calibration scan on first 20 pairs ───────────────────────────────────────
SUBSET = TESTABLE[:20]
print(f"\n=== Calibration scan on {len(SUBSET)} pairs ===")
results = []
t0 = time.time()
for i, (A, B) in enumerate(SUBSET):
    rec = {"gene_a": A, "gene_b": B}
    for outname, col in [("dist", "outcome_distance"),
                         ("cluster", "outcome_cluster")]:
        r = epsilon_4class(adata, A, B, col, n_boot=N_BOOT, seed=SEED)
        if r is None:
            continue
        rec[f"{outname}_eps"] = r["epsilon"]
        rec[f"{outname}_se"] = r["se"]
        rec[f"{outname}_z"] = r["z"]
        rec[f"{outname}_min_class"] = r["min_class"]
    results.append(rec)
    if (i + 1) % 5 == 0:
        print(f"  [{i+1}/{len(SUBSET)}] {A}↔{B} "
              f"ε_dist={rec.get('dist_eps', float('nan')):+.4f} "
              f"z_dist={rec.get('dist_z', float('nan')):+.2f} "
              f"({(time.time()-t0)/(i+1):.1f}s/pair)")
df_top = pd.DataFrame(results)
df_top.to_parquet(os.path.join(OUT_DIR, "top20_pairs.parquet"))
print(df_top[["gene_a", "gene_b", "dist_eps", "dist_z",
              "cluster_eps", "cluster_z"]].to_string(index=False))


# ── S1 Calibration via class permutation ─────────────────────────────────────
print("\n=== Sanity 1 — calibration via class permutation ===")
class_orig = adata.obs[["class_type", "gene_a", "gene_b"]].copy()
permuted_zs = []
n_calib = 50
rng_perm = np.random.default_rng(SEED)
t0 = time.time()
for trial in range(n_calib):
    A, B = TESTABLE[rng_perm.integers(0, len(TESTABLE))]
    perm = rng_perm.permutation(len(adata.obs))
    adata.obs[["class_type", "gene_a", "gene_b"]] = class_orig.iloc[perm].values
    r = epsilon_4class(adata, A, B, "outcome_distance",
                       n_boot=N_BOOT, seed=SEED)
    if r is not None:
        permuted_zs.append(r["z"])
adata.obs[["class_type", "gene_a", "gene_b"]] = class_orig.values
print(f"  ({time.time()-t0:.0f}s)")

permuted_zs = np.array(permuted_zs)
abs_z_p95 = float(np.percentile(np.abs(permuted_zs), 95))
print(f"random-class permutation null pairs: {len(permuted_zs)}")
print(f"|z|_p95 = {abs_z_p95:.2f}  (target [1.5, 2.5])")
SANITY1_PASS = 1.5 <= abs_z_p95 <= 2.5
print(f"Sanity 1: {'PASS' if SANITY1_PASS else 'FAIL'}")


# ── S2 Curated pair sensitivity ──────────────────────────────────────────────
print("\n=== Sanity 2 — curated pathway pair ===")
CURATED = [
    ("CBL", "CNN1"), ("FOXA1", "FOXA3"),
    ("CEBPA", "CEBPB"), ("CEBPB", "CEBPE"),
    ("TP53", "KRAS"),
    ("RPS3", "RPL5"),
]
sanity2_pair = sanity2_result = None
testable_set = {tuple(sorted(p)) for p in TESTABLE}
for (A, B) in CURATED:
    cand = tuple(sorted([A, B]))
    if cand in testable_set:
        r = epsilon_4class(adata, cand[0], cand[1], "outcome_distance",
                           n_boot=N_BOOT, seed=SEED)
        if r is not None:
            sanity2_pair = cand
            sanity2_result = r
            break
if sanity2_pair is None:
    df_sorted = (df_top.assign(abs_z=df_top["dist_z"].abs())
                 .sort_values("abs_z", ascending=False))
    if len(df_sorted):
        top1 = df_sorted.iloc[0]
        sanity2_pair = (top1["gene_a"], top1["gene_b"])
        sanity2_result = {"epsilon": top1["dist_eps"],
                          "se": top1["dist_se"], "z": top1["dist_z"]}

assert sanity2_pair is not None
A, B = sanity2_pair
print(f"curated/fallback pair: ({A}, {B})")
print(f"  ε  = {sanity2_result['epsilon']:+.4f}")
print(f"  SE = {sanity2_result['se']:.4f}")
print(f"  z  = {sanity2_result['z']:+.2f}")
SANITY2_PASS = abs(sanity2_result["z"]) > 3
print(f"Sanity 2: {'PASS' if SANITY2_PASS else 'FAIL'}")


# ── S3 Outcome consistency ───────────────────────────────────────────────────
print("\n=== Sanity 3 — outcome consistency ===")
eps_dist = df_top["dist_eps"].values
eps_clust = df_top["cluster_eps"].values
mask = np.isfinite(eps_dist) & np.isfinite(eps_clust)
rho_consist, p_consist = pearsonr(eps_dist[mask], eps_clust[mask])
print(f"pairs with both ε: {mask.sum()}/{len(df_top)}")
print(f"|Pearson ρ|(distance, cluster) = {abs(rho_consist):.3f}  (sign: {rho_consist:+.3f})")
SANITY3_PASS = abs(rho_consist) > 0.7
print(f"Sanity 3: {'PASS' if SANITY3_PASS else 'FAIL'}")


# ── S4 Sample size ───────────────────────────────────────────────────────────
print("\n=== Sanity 4 — sample size per class ===")
mins = df_top["dist_min_class"].values
overall_min = int(mins.min())
n_below = int((mins < 30).sum())
print(f"min class size: {overall_min}, pairs <30: {n_below}/{len(df_top)}")
SANITY4_PASS = (overall_min >= 30) and (n_below == 0)
print(f"Sanity 4: {'PASS' if SANITY4_PASS else 'FAIL'}")


# ── C3 SOFT correlation, C4 dominance ────────────────────────────────────────
print("\n=== Constraint 3 + 4 ===")
class_codes = pd.Categorical(adata.obs["class_type"]).codes
rho_c3, _ = pearsonr(class_codes, adata.obs["outcome_distance"].values)
max_abs_corr = abs(rho_c3)
n_eval = len(adata)
print(f"C3 |corr(class_indicator, outcome)| = {max_abs_corr:.3f}")

sig = df_top[df_top["dist_z"].abs() > 3]
median_z_sig = sig["dist_z"].abs().median() if len(sig) else 0.0
threshold_dom = max_abs_corr * np.sqrt(n_eval)
print(f"C4 median |z|_sig = {median_z_sig:.2f}, threshold = {threshold_dom:.2f}")
DOMINANCE_PASS = median_z_sig > threshold_dom
print(f"Dominance: {'PASS' if DOMINANCE_PASS else 'FAIL'}")


# ── Verdict ──────────────────────────────────────────────────────────────────
print("\n=== Verdict ===")
all_pass = bool(SANITY1_PASS and SANITY2_PASS and SANITY3_PASS
                and SANITY4_PASS and DOMINANCE_PASS)
verdict = {
    "project":           "BioEpistasis",
    "phase":             "Pilot iter 5 — Norman 2019 K562 Perturb-seq",
    "dataset":           "Norman et al. 2019 GSE133344 (scPerturb mirror)",
    "n_cells":           int(adata.n_obs),
    "n_testable_pairs":  len(TESTABLE),
    "n_calibration":     len(SUBSET),
    "sign_convention":   "ε in loss space; ε>0 = synthetic-lethal/redundancy",
    "sanity": {
        "S1_calibration_p95":       abs_z_p95,
        "S1_pass":                  bool(SANITY1_PASS),
        "S2_curated_pair":          list(sanity2_pair),
        "S2_z":                     sanity2_result["z"],
        "S2_pass":                  bool(SANITY2_PASS),
        "S3_pearson_rho":           float(rho_consist),
        "S3_pass":                  bool(SANITY3_PASS),
        "S4_min_class":             overall_min,
        "S4_pass":                  bool(SANITY4_PASS),
        "C3_soft_corr":             max_abs_corr,
        "C4_dominance_median_zsig": median_z_sig,
        "C4_dominance_threshold":   float(threshold_dom),
        "C4_dominance_pass":        bool(DOMINANCE_PASS),
    },
    "all_gates_pass": all_pass,
    "next_step": ("Lock pre-reg v1 biology side, run full Norman scan"
                  if all_pass else
                  "Halt; investigate failing gate before pre-reg drafting"),
}
# Cast np floats to native python types for JSON serialization
def _to_native(x):
    if isinstance(x, dict):
        return {k: _to_native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_to_native(v) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x

verdict = _to_native(verdict)
out = os.path.join(OUT_DIR, "iter5_verdict.json")
with open(out, "w") as f:
    json.dump(verdict, f, indent=2)
print(json.dumps(verdict, indent=2))
print()
print("=" * 70)
print("NORMAN PILOT iter 5:",
      "PASS — proceed to bio pre-reg" if all_pass else "FAIL — halt")
print("=" * 70)
print(f"Saved → {out}")
