"""
13: UMAP + GMM clustering on all 191 longitudinal subjects.

Prior clustering scripts (11, 12, 12c) were limited to the 57 subjects with ≥3
observations and ≥2 points within 200 days — an interpolation requirement.
This script avoids that bottleneck by using summary features computed directly
from each subject's data, enabling clustering on all 191 longitudinal subjects.

Feature vector (6 dimensions):
  1. log10(eu_first)    — starting titer (log10 EU)
  2. log10(eu_second)   — second-visit titer (log10 EU)
  3. log2(t_start+1)    — timing of first observation (days, log-scaled)
  4. log2_fc            — fold change first→second visit (log2 EU)
  5. p_resp             — posterior P(responder) from Step 4 mixture model
  6. is_typhi           — serovar indicator (1=Typhi, 0=Paratyphi)

Workflow:
  A. Build and standardize feature matrix for all 191 subjects
  B. Select GMM components k via BIC (k=1..8)
  C. Fit final GMM; assign cluster labels
  D. UMAP 2D embedding (colored by cluster, then by individual features)
  E. Cluster characterization plots

Run from repo root:
  reanalysis/.venv/bin/python3 \\
    metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/mixture_eda/13_umap_gmm_clustering.py
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
from importlib import import_module

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(".")
SCRIPT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/mixture_eda")
DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = SCRIPT_DIR

sys.path.insert(0, str(SCRIPT_DIR))

# ---------------------------------------------------------------------------
# Load Step 4 mixture model results (re-run 08's fit to get p_resp)
# ---------------------------------------------------------------------------

def load_features():
    """Build feature matrix for all 191 longitudinal subjects.

    Imports fold-change data from 08 and re-fits Step 4 to get p_resp.
    Merges with the raw longitudinal data for serovar and second-visit EU.
    """
    m08 = import_module("08_fold_change_mixture")

    # Fold-change dataframe: 191 subjects, first-two-point convention
    df = m08.load_fold_change_data()

    # Step 4: Teunis power-law mixture — refit to get p_noise posteriors
    # Replicates the main() call sequence: Step 1 → Step 4 (uses p1 as init)
    print("Fitting Step 4 mixture model...")
    x = df["log2_fc"].values
    t0 = df["t_start"].values
    t1 = df["t_end"].values
    p1 = m08.fit_two_gaussian(x)
    p4 = m08.fit_teunis_mixture(x, t0, t1, p1)
    df["p_resp"] = 1 - p4["p_noise"]

    # Pull eu_second from the raw long data
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    second_obs = {}
    serovar_map = {}
    for sid, grp in long[long["has_longitudinal"]].groupby("index_id"):
        grp_sorted = grp.sort_values("days_since_fever_onset")
        if len(grp_sorted) >= 2:
            second_obs[sid] = grp_sorted.iloc[1]["vi_igg_eu"]
            serovar_map[sid] = grp_sorted.iloc[0]["serovar"]

    df["eu_second"] = df["index_id"].map(second_obs)
    df["serovar"] = df["index_id"].map(serovar_map)
    df["is_typhi"] = (df["serovar"] == "typhi").astype(float)

    # Feature matrix
    feats = pd.DataFrame({
        "log10_eu_first":  np.log10(df["eu_start"]),
        "log10_eu_second": np.log10(df["eu_second"]),
        "log2_t_start":    np.log2(df["t_start"] + 1),
        "log2_fc":         df["log2_fc"],
        "p_resp":          df["p_resp"],
        "is_typhi":        df["is_typhi"],
    })

    return df, feats


def standardize(feats):
    """Z-score each feature column."""
    mu = feats.mean()
    sd = feats.std(ddof=0)
    sd = sd.where(sd > 0, 1.0)  # avoid /0
    return (feats - mu) / sd, mu, sd


# ---------------------------------------------------------------------------
# GMM BIC selection
# ---------------------------------------------------------------------------

def select_gmm_k(X_std, k_max=8, n_init=20, random_state=42):
    """Fit GMM for k=1..k_max, return BIC array and fitted models."""
    from sklearn.mixture import GaussianMixture
    bics = []
    models = []
    for k in range(1, k_max + 1):
        gm = GaussianMixture(n_components=k, n_init=n_init,
                             covariance_type="full", random_state=random_state)
        gm.fit(X_std)
        bics.append(gm.bic(X_std))
        models.append(gm)
    return np.array(bics), models


def plot_bic(bics, k_best, out_path):
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ks = np.arange(1, len(bics) + 1)
    ax.plot(ks, bics, "o-", color="#2166ac")
    ax.axvline(k_best, color="#d6604d", ls="--", lw=1.5, label=f"k={k_best} selected")
    ax.set_xlabel("Number of GMM components (k)")
    ax.set_ylabel("BIC")
    ax.set_title("GMM model selection by BIC")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {out_path}")


# ---------------------------------------------------------------------------
# UMAP embedding
# ---------------------------------------------------------------------------

def run_umap(X_std, random_state=42):
    import umap
    reducer = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1,
                        random_state=random_state)
    embedding = reducer.fit_transform(X_std)
    return embedding


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

CLUSTER_CMAP = plt.cm.get_cmap("tab10")


def _scatter_base(ax, emb, c, cmap, vmin=None, vmax=None, s=35, alpha=0.8,
                  label=None, norm=None):
    sc = ax.scatter(emb[:, 0], emb[:, 1], c=c, cmap=cmap,
                    vmin=vmin, vmax=vmax, s=s, alpha=alpha, norm=norm)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    return sc


def plot_umap_overview(emb, labels, df, feats, out_path):
    """4-panel UMAP: cluster label, P(resp), log10(EU start), serovar."""
    k = labels.max() + 1
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel 1: GMM cluster
    ax = axes[0, 0]
    colors = [CLUSTER_CMAP(i / max(k - 1, 1)) for i in labels]
    for cl in range(k):
        mask = labels == cl
        ax.scatter(emb[mask, 0], emb[mask, 1], s=35, alpha=0.8,
                   color=CLUSTER_CMAP(cl / max(k - 1, 1)), label=f"Cluster {cl + 1}")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title("GMM cluster assignment")
    ax.legend(markerscale=1.2, fontsize=8)

    # Panel 2: P(responder)
    ax = axes[0, 1]
    sc = _scatter_base(ax, emb, df["p_resp"].values, "RdYlBu_r", vmin=0, vmax=1)
    ax.set_title("P(responder) from Step 4 mixture")
    plt.colorbar(sc, ax=ax, label="P(resp)")

    # Panel 3: log10(EU first)
    ax = axes[1, 0]
    sc = _scatter_base(ax, emb, feats["log10_eu_first"].values, "viridis")
    ax.set_title("log₁₀(EU first visit)")
    plt.colorbar(sc, ax=ax, label="log₁₀(EU)")

    # Panel 4: serovar
    ax = axes[1, 1]
    typhi_mask = df["is_typhi"].values == 1
    ax.scatter(emb[typhi_mask, 0], emb[typhi_mask, 1], s=35, alpha=0.8,
               color="#d73027", label="S. Typhi")
    ax.scatter(emb[~typhi_mask, 0], emb[~typhi_mask, 1], s=35, alpha=0.8,
               color="#4575b4", label="S. Paratyphi")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title("Serovar")
    ax.legend()

    fig.suptitle("UMAP of 191 longitudinal subjects\n(6-feature GMM clustering)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {out_path}")


def plot_umap_features(emb, feats, df, out_path):
    """6-panel UMAP, one per raw feature."""
    feat_info = [
        ("log10_eu_first",  "log₁₀(EU first visit)",  "viridis",  None, None),
        ("log10_eu_second", "log₁₀(EU second visit)", "viridis",  None, None),
        ("log2_t_start",    "log₂(t_start + 1 days)", "plasma",   None, None),
        ("log2_fc",         "log₂(fold change)",       "coolwarm", -3,   3),
        ("p_resp",          "P(responder)",             "RdYlBu_r", 0,   1),
        ("is_typhi",        "Is Typhi (1=yes)",         "bwr",      0,   1),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.ravel()

    for ax, (col, title, cmap, vmin, vmax) in zip(axes, feat_info):
        vals = feats[col].values
        sc = ax.scatter(emb[:, 0], emb[:, 1], c=vals, cmap=cmap,
                        vmin=vmin, vmax=vmax, s=30, alpha=0.8)
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.set_title(title)
        plt.colorbar(sc, ax=ax)

    fig.suptitle("UMAP feature panels — 191 longitudinal subjects",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {out_path}")


def plot_cluster_summary(df, labels, feats, out_path):
    """Box/strip plots of key features by cluster."""
    k = labels.max() + 1
    df = df.copy()
    df["cluster"] = labels + 1  # 1-indexed labels

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.ravel()

    panels = [
        ("log10_eu_first",  feats["log10_eu_first"].values,  "log₁₀(EU first visit)"),
        ("log2_fc",         feats["log2_fc"].values,          "log₂(fold change)"),
        ("p_resp",          feats["p_resp"].values,            "P(responder)"),
        ("log10_eu_second", feats["log10_eu_second"].values,  "log₁₀(EU second visit)"),
        ("log2_t_start",    feats["log2_t_start"].values,     "log₂(t_start+1 days)"),
        ("is_typhi",        feats["is_typhi"].values,          "Fraction Typhi"),
    ]

    for ax, (col, vals, ylabel) in zip(axes, panels):
        df[col] = vals
        cl_labels = sorted(df["cluster"].unique())

        for i, cl in enumerate(cl_labels):
            sub = df[df["cluster"] == cl][col].values
            color = CLUSTER_CMAP(i / max(k - 1, 1))

            # Jitter strip
            jx = np.random.default_rng(i).uniform(-0.2, 0.2, size=len(sub))
            ax.scatter(np.full(len(sub), cl) + jx, sub,
                       color=color, s=10, alpha=0.4)
            # Box
            bp = ax.boxplot(sub, positions=[cl], widths=0.4, patch_artist=True,
                            medianprops=dict(color="black", lw=2),
                            whiskerprops=dict(color=color),
                            capprops=dict(color=color),
                            flierprops=dict(visible=False))
            bp["boxes"][0].set_facecolor((*mcolors.to_rgb(color), 0.3))
            bp["boxes"][0].set_edgecolor(color)

        ax.set_xlabel("Cluster")
        ax.set_ylabel(ylabel)
        ax.set_xticks(cl_labels)

    # Cluster size as text on first axis
    axes[0].set_title("log₁₀(EU first visit) by cluster")
    for cl in sorted(df["cluster"].unique()):
        n = (df["cluster"] == cl).sum()
        y_pos = df[df["cluster"] == cl]["log10_eu_first"].max()
        axes[0].text(cl, y_pos + 0.05, f"n={n}", ha="center", va="bottom",
                     fontsize=8, color="black")

    fig.suptitle("Cluster feature distributions — 191 longitudinal subjects",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {out_path}")


def plot_cluster_trajectories(df_long, df_feats, labels, out_path):
    """Trajectory spaghetti plots by cluster, using all available visits."""
    from statsmodels.nonparametric.smoothers_lowess import lowess

    # Merge cluster labels back onto full longitudinal data
    df_feats = df_feats.copy()
    df_feats["cluster"] = labels + 1  # 1-indexed

    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    merged = long[long["has_longitudinal"]].merge(
        df_feats[["index_id", "cluster"]], on="index_id", how="inner"
    )

    k = labels.max() + 1
    n_cols = min(k, 4)
    n_rows = (k + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 4 * n_rows),
                             squeeze=False)
    axes_flat = axes.ravel()

    for i, cl in enumerate(sorted(merged["cluster"].unique())):
        ax = axes_flat[i]
        sub = merged[merged["cluster"] == cl]
        color = CLUSTER_CMAP(i / max(k - 1, 1))

        n_subj = sub["index_id"].nunique()
        typhi_frac = (sub.drop_duplicates("index_id")["serovar"] == "typhi").mean()
        p_resp_median = df_feats[df_feats["cluster"] == cl]["p_resp"].median()

        # Individual trajectories
        for sid, grp in sub.groupby("index_id"):
            grp = grp.sort_values("days_since_fever_onset")
            ax.plot(grp["days_since_fever_onset"], np.log10(grp["vi_igg_eu"]),
                    color=color, alpha=0.25, lw=0.8)

        # LOESS overlay (all data points in cluster)
        t_all = sub["days_since_fever_onset"].values
        y_all = np.log10(sub["vi_igg_eu"].values)
        if len(t_all) >= 6:
            order = np.argsort(t_all)
            smoothed = lowess(y_all[order], t_all[order], frac=0.5, return_sorted=True)
            ax.plot(smoothed[:, 0], smoothed[:, 1], color="black", lw=2, alpha=0.85,
                    label="LOESS")

        ax.set_xlim(-5, 400)
        ax.set_ylim(1, 4.5)
        ax.set_xlabel("Days since fever onset")
        ax.set_ylabel("log₁₀(Vi IgG EU)")
        ax.set_title(
            f"Cluster {cl}  (n={n_subj})\n"
            f"Typhi {100*typhi_frac:.0f}%  |  P(resp)={p_resp_median:.2f}"
        )

    # Hide unused axes
    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle("Vi IgG trajectories by GMM cluster\n(all available visits)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("13: UMAP + GMM clustering — 191 longitudinal subjects")
    print("=" * 60)

    # ---- Feature matrix --------------------------------------------------
    print("\nBuilding feature matrix...")
    df, feats = load_features()
    print(f"  {len(df)} subjects, {feats.shape[1]} features")
    print(feats.describe().to_string())

    X_std, feat_mu, feat_sd = standardize(feats)
    X_arr = X_std.values

    # ---- GMM BIC selection -----------------------------------------------
    print("\nSelecting k via BIC...")
    bics, models = select_gmm_k(X_arr, k_max=8, n_init=30)
    k_best = int(np.argmin(bics)) + 1
    print(f"  BICs: {[f'{b:.1f}' for b in bics]}")
    print(f"  Best k = {k_best}")
    plot_bic(bics, k_best, OUT_DIR / "13_gmm_bic.png")

    # Fit final GMM
    from sklearn.mixture import GaussianMixture
    gm_best = models[k_best - 1]
    labels = gm_best.predict(X_arr)
    probs = gm_best.predict_proba(X_arr)

    # Sort clusters by median P(resp) descending
    cluster_presp = [df["p_resp"].values[labels == k].mean() for k in range(k_best)]
    order = np.argsort(cluster_presp)[::-1]
    remap = {old: new for new, old in enumerate(order)}
    labels = np.array([remap[l] for l in labels])

    print(f"\nCluster summary (sorted by mean P(resp)):")
    for cl in range(k_best):
        mask = labels == cl
        n = mask.sum()
        typhi_frac = feats["is_typhi"].values[mask].mean()
        p_resp = df["p_resp"].values[mask].mean()
        log2fc = np.median(feats["log2_fc"].values[mask])
        log10eu = np.median(feats["log10_eu_first"].values[mask])
        print(f"  Cluster {cl+1}: n={n}, P(resp)={p_resp:.2f}, "
              f"Typhi={100*typhi_frac:.0f}%, "
              f"median log2FC={log2fc:.2f}, median log10EU={log10eu:.2f}")

    # ---- UMAP embedding --------------------------------------------------
    print("\nRunning UMAP...")
    emb = run_umap(X_arr)

    # ---- Plots -----------------------------------------------------------
    print("\nGenerating figures...")
    df["cluster"] = labels + 1  # 1-indexed for plotting

    plot_umap_overview(emb, labels, df, feats, OUT_DIR / "13_umap_overview.png")
    plot_umap_features(emb, feats, df, OUT_DIR / "13_umap_features.png")
    plot_cluster_summary(df, labels, feats, OUT_DIR / "13_cluster_summary.png")
    plot_cluster_trajectories(df, df[["index_id"]].assign(p_resp=df["p_resp"],
                                                           cluster=labels + 1),
                              labels, OUT_DIR / "13_cluster_trajectories.png")

    print("\nDone.")


if __name__ == "__main__":
    main()
