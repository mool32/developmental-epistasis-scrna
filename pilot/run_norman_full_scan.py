"""
Full Norman 2019 K562 Perturb-seq scan over all 131 testable pairs.

Pre-registered as v1 (commit 241d22f, tag bio_prereg_v1_locked) in
mool32/developmental-epistasis-scrna.

LOCKED design:
- Single primary outcome: PCA distance to control-cluster centroid (50 PCs).
- 4-class intervention contrast (control / A only / B only / AB pair).
- Per-pair bootstrap n=1000, seed=42.
- Primary statistic: frac_synth = #{|z|>3 AND ε>0} / #{|z|>3}.
- Permutation null: 10,000 iterations, sign-flip per pair, seed 20260504.
- 4-tier decision: PASS / PARTIAL / WEAK / FAIL_REVERSED.
- 5 pre-flight gates (G1-G5) — all must hold or HOLD primary verdict.

Compute: ~10-15 min CPU.

Output:
  data/analysis/norman/norman_pairs.parquet
  data/analysis/norman/bio_norman_verdict.json
  data/analysis/norman/bio_norman_headline.png
"""

import os, json, time, subprocess
from itertools import combinations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, kstest, halfnorm, mannwhitneyu, norm, t as student_t, laplace
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


print("=" * 70)
print("Norman 2019 K562 Perturb-seq — FULL SCAN per LOCKED pre-reg v1")
print("Pre-reg commit: 241d22f, tag: bio_prereg_v1_locked")
print()
print("SIGN CONVENTION (locked):")
print("  ε > 0 ↔ synthetic-lethal / redundancy phenotype")
print("  ε < 0 ↔ suppression / buffering")
print("=" * 70)
print()

SEED = 42
PERMUTATION_SEED = 20260504
N_BOOT = 1000
N_PERM = 10000
PRE_REG_COMMIT = "241d22f"
PRE_REG_TAG = "bio_prereg_v1_locked"

OUT_DIR = os.path.expanduser("~/BioEpistasis_pilot/full_scan_norman")
DATA_DIR = os.path.join(OUT_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)


# ── Load Norman h5ad ─────────────────────────────────────────────────────────
H5AD = os.path.expanduser("~/BioEpistasis_pilot/iter5_norman/data/NormanWeissman2019_filtered.h5ad")
if not os.path.exists(H5AD):
    URL = "https://zenodo.org/records/10044268/files/NormanWeissman2019_filtered.h5ad?download=1"
    print(f"Downloading {URL}")
    os.makedirs(os.path.dirname(H5AD), exist_ok=True)
    subprocess.check_call(["curl", "-L", "-o", H5AD, URL])

print(f"Loading {os.path.getsize(H5AD)/1e6:.0f} MB h5ad ...")
adata = sc.read_h5ad(H5AD)
print(adata)


# ── Parse classes (per pre-reg §1, scPerturb format) ─────────────────────────
def parse_row(g, n):
    if g == "control" or n == 0:
        return ("control", None, None)
    if n == 1:
        return ("single", g, None)
    if n == 2:
        parts = g.split("_")
        if len(parts) == 2:
            return ("pair", parts[0], parts[1])
        return ("pair", parts[0], "_".join(parts[1:]))
    return ("unknown", None, None)


class_info = [parse_row(str(g), int(n))
              for g, n in zip(adata.obs["perturbation"], adata.obs["nperts"])]
adata.obs["class_type"] = [c[0] for c in class_info]
adata.obs["gene_a"] = [c[1] for c in class_info]
adata.obs["gene_b"] = [c[2] for c in class_info]


# ── Preprocess + cluster (per pre-reg §3.1) ──────────────────────────────────
print("\n=== Preprocess + cluster (locked design) ===")
if adata.X.max() > 100:
    sc.pp.filter_genes(adata, min_cells=20)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=2000)
sc.pp.pca(adata, n_comps=50, random_state=SEED)
sc.pp.neighbors(adata, n_neighbors=15, random_state=SEED)
sc.tl.leiden(adata, resolution=1.0, random_state=SEED)

ctrl_clust = adata.obs.loc[adata.obs["class_type"] == "control", "leiden"].mode().iloc[0]
ctrl_clust_mask = (adata.obs["leiden"] == ctrl_clust)
ctrl_centroid = adata.obsm["X_pca"][ctrl_clust_mask].mean(axis=0)
diff = adata.obsm["X_pca"] - ctrl_centroid[None, :]
adata.obs["outcome_distance"] = np.linalg.norm(diff, axis=1)
print(f"Control-dominant cluster: {ctrl_clust}, "
      f"outcome_distance: mean={adata.obs['outcome_distance'].mean():.2f}")


# ── Identify all testable pairs (≥30 cells in each of 4 classes) ─────────────
print("\n=== Identify testable pairs (G3 sample size gate) ===")
single_counts = adata.obs[adata.obs["class_type"] == "single"]["gene_a"].value_counts()
pair_counts = (adata.obs[adata.obs["class_type"] == "pair"]
               .groupby(["gene_a", "gene_b"]).size().reset_index(name="count"))
pair_counts["pair_key"] = pair_counts.apply(
    lambda r: tuple(sorted([r["gene_a"], r["gene_b"]])), axis=1)
pair_counts = pair_counts.groupby("pair_key")["count"].sum().reset_index()

n_ctrl = (adata.obs["class_type"] == "control").sum()
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
print(f"Total testable pairs: {len(TESTABLE)} (control n={n_ctrl})")


# ── ε helper ─────────────────────────────────────────────────────────────────
def epsilon_4class(adata, gene_a, gene_b, n_boot=N_BOOT, seed=SEED):
    pair_canon = tuple(sorted([gene_a, gene_b]))
    obs = adata.obs
    out = obs["outcome_distance"].values
    mask_ctrl = (obs["class_type"] == "control").values
    mask_a = ((obs["class_type"] == "single") & (obs["gene_a"] == gene_a)).values
    mask_b = ((obs["class_type"] == "single") & (obs["gene_a"] == gene_b)).values
    mask_ab = ((obs["class_type"] == "pair") & obs.apply(
        lambda r: tuple(sorted([str(r["gene_a"]), str(r["gene_b"])])) == pair_canon,
        axis=1)).values
    arr_ctrl, arr_a, arr_b, arr_ab = out[mask_ctrl], out[mask_a], out[mask_b], out[mask_ab]
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
            "n_corners": n, "min_class": min(n.values()),
            "delta_a": float(np.mean(arr_a) - np.mean(arr_ctrl)),
            "delta_b": float(np.mean(arr_b) - np.mean(arr_ctrl)),
            "delta_ab": float(np.mean(arr_ab) - np.mean(arr_ctrl))}


# ── Full scan ────────────────────────────────────────────────────────────────
print(f"\n=== Full scan over {len(TESTABLE)} pairs ===")
results = []
t0 = time.time()
for i, (A, B) in enumerate(TESTABLE):
    r = epsilon_4class(adata, A, B)
    if r is None:
        continue
    rec = {"gene_a": A, "gene_b": B,
           "epsilon": r["epsilon"], "se": r["se"], "z": r["z"],
           "delta_a": r["delta_a"], "delta_b": r["delta_b"],
           "delta_ab": r["delta_ab"],
           "n_ctrl": r["n_corners"]["ctrl"], "n_a": r["n_corners"]["a"],
           "n_b": r["n_corners"]["b"], "n_ab": r["n_corners"]["ab"]}
    results.append(rec)
    if (i + 1) % 25 == 0 or (i + 1) == len(TESTABLE):
        rate = (time.time() - t0) / (i + 1)
        eta = rate * (len(TESTABLE) - i - 1)
        print(f"  [{i+1}/{len(TESTABLE)}] {A}↔{B} "
              f"ε={r['epsilon']:+.3f} z={r['z']:+.2f} "
              f"({rate:.1f}s/pair, ETA {eta:.0f}s)")
df = pd.DataFrame(results)
df.to_parquet(os.path.join(OUT_DIR, "norman_pairs.parquet"))
print(f"\nFull scan complete: {len(df)} pairs in {(time.time()-t0)/60:.1f} min")


# ── Pre-flight gates (G1-G5) ─────────────────────────────────────────────────
print("\n=== Pre-flight gates per pre-reg §5 ===")

# G1 calibration via class permutation
print("G1 calibration ...", end=" ")
class_orig = adata.obs[["class_type", "gene_a", "gene_b"]].copy()
permuted_zs = []
rng_perm = np.random.default_rng(SEED)
for trial in range(50):
    A, B = TESTABLE[rng_perm.integers(0, len(TESTABLE))]
    perm = rng_perm.permutation(len(adata.obs))
    adata.obs[["class_type", "gene_a", "gene_b"]] = class_orig.iloc[perm].values
    r = epsilon_4class(adata, A, B)
    if r is not None:
        permuted_zs.append(r["z"])
adata.obs[["class_type", "gene_a", "gene_b"]] = class_orig.values
abs_z_p95 = float(np.percentile(np.abs(permuted_zs), 95))
G1 = 1.5 <= abs_z_p95 <= 2.5
print(f"|z|_p95 = {abs_z_p95:.2f}  {'PASS' if G1 else 'FAIL'}")

# G2 sensitivity: at least one curated pair |z|>3
print("G2 sensitivity ...", end=" ")
CURATED = [("CBL", "CNN1"), ("CEBPA", "CEBPB"), ("CEBPB", "CEBPE")]
g2_results = []
for (A, B) in CURATED:
    cand = tuple(sorted([A, B]))
    matches = df[(df["gene_a"] == cand[0]) & (df["gene_b"] == cand[1])]
    if len(matches):
        z_val = matches.iloc[0]["z"]
        g2_results.append((cand, z_val))
G2 = any(abs(z) > 3 for _, z in g2_results)
print(f"max |z| on curated = {max(abs(z) for _, z in g2_results) if g2_results else 0:.2f}  "
      f"{'PASS' if G2 else 'FAIL'}")

# G3 sample size
print("G3 sample size ...", end=" ")
min_class = int(df[["n_ctrl", "n_a", "n_b", "n_ab"]].min().min())
G3 = min_class >= 30
print(f"min class = {min_class}  {'PASS' if G3 else 'FAIL'}")

# G4 SOFT correlation
print("G4 SOFT correlation ...", end=" ")
class_codes = pd.Categorical(adata.obs["class_type"]).codes
rho_g4, _ = pearsonr(class_codes, adata.obs["outcome_distance"].values)
G4 = abs(rho_g4) < 0.05
print(f"|corr| = {abs(rho_g4):.4f}  {'PASS' if G4 else 'FAIL'}")

# G5 dominance
print("G5 dominance ...", end=" ")
sig = df[df["z"].abs() > 3]
median_z_sig = float(sig["z"].abs().median()) if len(sig) else 0.0
threshold_dom = abs(rho_g4) * np.sqrt(len(adata))
G5 = median_z_sig > threshold_dom
print(f"median |z|_sig = {median_z_sig:.2f} vs threshold {threshold_dom:.4f}  "
      f"{'PASS' if G5 else 'FAIL'}")

ALL_GATES = G1 and G2 and G3 and G4 and G5
print(f"\nAll 5 gates: {'PASS' if ALL_GATES else 'HOLD'}")


# ── Primary test ─────────────────────────────────────────────────────────────
print("\n=== Primary test per pre-reg §4 ===")
sig = df[df["z"].abs() > 3].copy()
n_sig = len(sig)
n_synth = int((sig["epsilon"] > 0).sum())
frac_synth = n_synth / max(1, n_sig)
print(f"Significant pairs (|z|>3): {n_sig}/{len(df)}")
print(f"frac_synth (#{{|z|>3 AND ε>0}} / #{{|z|>3}}) = {frac_synth:.3f}")

# Permutation null
print(f"Computing permutation null ({N_PERM} iter, seed {PERMUTATION_SEED}) ...")
rng = np.random.default_rng(PERMUTATION_SEED)
perm_fracs = np.empty(N_PERM)
all_eps = sig["epsilon"].values  # only operate on |z|>3 subset
all_z = sig["z"].abs().values  # not needed for sign-flip null
for k in range(N_PERM):
    flips = rng.choice([-1, 1], size=len(all_eps))
    perm_eps = all_eps * flips
    perm_fracs[k] = (perm_eps > 0).sum() / max(1, len(perm_eps))
p_value = float(np.mean(np.abs(perm_fracs - 0.5) >= np.abs(frac_synth - 0.5)))
print(f"Permutation p (two-sided, vs symmetric null at 0.5): {p_value:.4f}")

# Decision tier
if frac_synth > 0.55 and p_value < 0.01:
    TIER = "PASS"
elif 0.50 < frac_synth <= 0.55 and p_value < 0.05:
    TIER = "PARTIAL"
elif 0.45 <= frac_synth <= 0.50 and p_value < 0.10:
    TIER = "WEAK"
elif frac_synth < 0.45:
    TIER = "FAIL_REVERSED"
else:
    TIER = "INDETERMINATE"
print(f"\nPRIMARY VERDICT: {TIER}")


# ── Secondaries ──────────────────────────────────────────────────────────────
print("\n=== Secondaries per pre-reg §6 ===")

# 6.1 Magnitude
median_abs_eps_sig = float(sig["epsilon"].abs().median()) if len(sig) else 0.0
median_abs_eps_all = float(df["epsilon"].abs().median())
print(f"6.1 median |ε|_significant = {median_abs_eps_sig:.4f}")
print(f"    median |ε|_all          = {median_abs_eps_all:.4f}")

# 6.2 Distribution shape AIC
abs_z = df["z"].abs().values
def aic_norm(x):
    mu, sd = x.mean(), x.std(ddof=1)
    return 2*2 - 2*norm.logpdf(x, mu, sd).sum()
def aic_lap(x):
    loc, sc_ = laplace.fit(x); return 2*2 - 2*laplace.logpdf(x, loc, sc_).sum()
def aic_t(x):
    df_, loc, sc_ = student_t.fit(x); return 2*3 - 2*student_t.logpdf(x, df_, loc, sc_).sum()
aics = {"gaussian": aic_norm(abs_z), "laplace": aic_lap(abs_z),
        "student_t": aic_t(abs_z)}
best_shape = min(aics, key=aics.get)
print(f"6.2 AIC shape: {aics}, best = {best_shape}")

# 6.3 Pathway coherence — skipped (no pathway annotations available locally)
print("6.3 pathway coherence: SKIPPED (no KEGG/GO mapping in pilot scope)")

# 6.4 Cross-substrate comparison
print("6.4 Cross-substrate frac(ε>0) of significant pairs:")
print(f"    Pythia 410M Tier 1:  0.78")
print(f"    OLMo 1B Phase 4:     0.57")
print(f"    Yeast Costanzo 2010: ~0.60-0.70")
print(f"    Norman 2019 K562:    {frac_synth:.3f}  ← THIS WORK")

# 6.5 Top-K identification
top30 = df.assign(abs_z=df["z"].abs()).sort_values("abs_z", ascending=False).head(30)
print(f"6.5 Top-30 by |z| computed (saved to verdict)")


# ── Verdict JSON ─────────────────────────────────────────────────────────────
verdict = {
    "pre_registration_tag":  PRE_REG_TAG,
    "pre_reg_commit":        PRE_REG_COMMIT,
    "model":                 "Norman 2019 K562 CRISPRi",
    "dataset":               "scPerturb mirror, NormanWeissman2019_filtered.h5ad",
    "n_cells":               int(adata.n_obs),
    "n_testable_pairs":      len(TESTABLE),
    "n_pairs_scanned":       len(df),
    "sign_convention":       "ε in loss space; ε>0 = synthetic-lethal/redundancy",
    "outcome":               "PCA distance to control-cluster centroid (50 PCs)",
    "n_boot":                N_BOOT,
    "n_permutations":        N_PERM,
    "permutation_seed":      PERMUTATION_SEED,

    "gates": {
        "G1_calibration_p95":           abs_z_p95,
        "G1_pass":                      bool(G1),
        "G2_max_curated_z":             max(abs(z) for _, z in g2_results) if g2_results else 0.0,
        "G2_pass":                      bool(G2),
        "G3_min_class":                 min_class,
        "G3_pass":                      bool(G3),
        "G4_soft_corr":                 abs(rho_g4),
        "G4_pass":                      bool(G4),
        "G5_dominance_median_zsig":     median_z_sig,
        "G5_dominance_threshold":       float(threshold_dom),
        "G5_pass":                      bool(G5),
        "all_gates_pass":               bool(ALL_GATES),
    },

    "primary": {
        "n_significant_pairs":  int(n_sig),
        "n_synth":              int(n_synth),
        "frac_synth":           float(frac_synth),
        "permutation_p":        float(p_value),
        "tier":                 TIER,
    },

    "secondaries": {
        "median_abs_eps_significant": median_abs_eps_sig,
        "median_abs_eps_all":         median_abs_eps_all,
        "aic_shape":                  {k: float(v) for k, v in aics.items()},
        "best_shape":                 best_shape,
        "cross_substrate_frac_synth": {
            "pythia_410m_tier1":      0.78,
            "olmo_1b_phase4":         0.57,
            "yeast_costanzo2010":     "0.60-0.70 (literature range)",
            "norman_2019_k562":       float(frac_synth),
        },
        "top30_by_abs_z":             top30[["gene_a","gene_b","epsilon","se","z"]].to_dict("records"),
    },
}

# Cast np types
def _to_native(x):
    if isinstance(x, dict): return {k: _to_native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [_to_native(v) for v in x]
    if isinstance(x, (np.floating, np.integer)): return x.item()
    if isinstance(x, np.ndarray): return x.tolist()
    return x

verdict = _to_native(verdict)
out_json = os.path.join(OUT_DIR, "bio_norman_verdict.json")
with open(out_json, "w") as f:
    json.dump(verdict, f, indent=2)
print(f"\nVerdict saved → {out_json}")


# ── Headline figure ──────────────────────────────────────────────────────────
print("\n=== Headline figure ===")
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Panel A: ε distribution of significant pairs
ax = axes[0]
ax.hist(sig["epsilon"], bins=25, edgecolor="black", alpha=0.75, color="C2")
ax.axvline(0, color="red", linestyle="--", lw=1.5)
ax.set_xlabel(r"$\varepsilon$ (significant pairs only, $|z|>3$)")
ax.set_ylabel("count")
ax.set_title(f"Norman 2019 K562 — {n_sig} significant pairs\n"
             f"frac(ε>0) = {frac_synth:.3f} → {TIER}")

# Panel B: cross-substrate comparison
ax = axes[1]
substrates = ["Pythia 410M\n(Tier 1)", "OLMo 1B\n(Phase 4)",
              "Yeast\n(Costanzo 2010)", "Norman K562\n(THIS WORK)"]
fracs = [0.78, 0.57, 0.65, frac_synth]  # yeast midpoint of 0.60-0.70
colors = ["#d62728", "#ff7f0e", "#9467bd", "#2ca02c"]
bars = ax.bar(substrates, fracs, color=colors,
              edgecolor="black", linewidth=0.8)
ax.axhline(0.5, color="gray", linestyle=":", alpha=0.6, label="symmetric null")
ax.axhline(0.55, color="orange", linestyle=":", alpha=0.6, label="PASS threshold (0.55)")
ax.set_ylabel(r"frac($\varepsilon > 0$) of significant pairs")
ax.set_ylim(0, 1)
ax.set_title("Cross-substrate synthetic-lethal regime")
for bar, f_ in zip(bars, fracs):
    ax.text(bar.get_x() + bar.get_width()/2, f_ + 0.02,
            f"{f_:.2f}", ha="center", fontsize=10, fontweight="bold")
ax.legend(loc="upper left", fontsize=9)

fig.suptitle(f"Bio pre-reg v1 verdict — Norman 2019 full scan ({len(df)} pairs)",
             fontsize=13, y=1.00)
plt.tight_layout()
out_png = os.path.join(OUT_DIR, "bio_norman_headline.png")
plt.savefig(out_png, dpi=170, bbox_inches="tight", facecolor="white")
print(f"Figure saved → {out_png}")

print()
print("=" * 70)
print(f"NORMAN FULL SCAN COMPLETE  →  TIER: {TIER}")
print(f"frac_synth = {frac_synth:.3f}  (predicted >0.55 for PASS)")
print(f"perm p     = {p_value:.4f}    (target <0.01 for PASS)")
print(f"all gates  = {'PASS' if ALL_GATES else 'HOLD'}")
print("=" * 70)
