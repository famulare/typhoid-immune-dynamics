"""
12: Enriched trajectory clustering — shape + covariates.

Extends script 11's shape-only clustering by adding:
  - P(resp) from Step 1 mixture model (fold-change signal)
  - Serovar indicator (typhi=1, paratyphi=0)

Feature vector per subject (all standardized to unit variance):
  - Shape vector: mean-centered, L2-normed interpolated log10(EU)
  - P(resp): scalar
  - is_typhi: binary indicator

Distance: Euclidean on the concatenated standardized features.
Clustering: Ward hierarchical.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist
from statsmodels.nonparametric.smoothers_lowess import lowess as sm_lowess
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from importlib import import_module
m08 = import_module("08_fold_change_mixture")

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/mixture_eda")

T_GRID = np.array([5, 15, 30, 60, 90, 120, 150, 180])
OUTLIER_IDS = {"indx387"}


def main():
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    raw = pd.read_csv(DATA_DIR / "TypoidCaseData_github_09.30.21.csv")
    serovar_map = raw.set_index("index_id")["bldculres"]

    # Fit working model for P(resp)
    df_fc = m08.load_fold_change_data()
    x = df_fc["log2_fc"].values
    p1 = m08.fit_two_gaussian(x)
    p_resp_map = dict(zip(df_fc["index_id"], 1 - p1["p_noise"]))

    # Select subjects: ≥3 obs, ≥2 within 200d, not outlier
    longi = long[long["has_longitudinal"]]
    nobs = longi.groupby("index_id").size()
    multi_ids = set(nobs[nobs >= 3].index) - OUTLIER_IDS
    good_ids = set()
    for sid in multi_ids:
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        if (grp["days_since_fever_onset"] <= 200).sum() >= 2:
            good_ids.add(sid)

    # Build feature matrix: shape + P(resp) + serovar
    feature_rows = {}
    for sid in sorted(good_ids):
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        t = grp["days_since_fever_onset"].values
        y = np.log10(np.maximum(grp["vi_igg_eu"].values, 1))

        t_min, t_max = t.min(), t.max()
        grid_mask = (T_GRID >= t_min) & (T_GRID <= t_max)
        if grid_mask.sum() < 3:
            continue

        t_use = T_GRID[grid_mask]
        y_interp = np.interp(t_use, t, y)
        y_centered = y_interp - y_interp.mean()
        norm = np.linalg.norm(y_centered)
        if norm < 1e-10:
            continue
        y_normed = y_centered / norm

        # Pad to full grid length (NaN for uncovered points)
        shape_full = np.full(len(T_GRID), np.nan)
        shape_full[grid_mask] = y_normed

        p_resp = p_resp_map.get(sid, np.nan)
        is_typhi = 1.0 if serovar_map.get(sid) == "typhi" else 0.0

        if np.isnan(p_resp):
            continue

        feature_rows[sid] = {
            "shape": shape_full,
            "grid_mask": grid_mask,
            "p_resp": p_resp,
            "is_typhi": is_typhi,
            "serovar": serovar_map.get(sid, "unknown"),
            "n_obs": len(grp),
            "subj_mean_log10": y.mean(),
        }

    sids = sorted(feature_rows.keys())
    n = len(sids)
    print(f"Subjects with complete features: {n}")

    # Build feature matrix — handle different grid coverages
    # Use only grid points covered by ALL subjects for the shape component,
    # OR compute pairwise distances on overlapping points
    # Simpler approach: fill NaN with 0 (neutral after centering+norming)
    # and standardize each feature dimension
    shape_matrix = np.array([feature_rows[s]["shape"] for s in sids])
    shape_matrix = np.nan_to_num(shape_matrix, nan=0.0)

    p_resp_vec = np.array([feature_rows[s]["p_resp"] for s in sids])
    typhi_vec = np.array([feature_rows[s]["is_typhi"] for s in sids])

    # Standardize each feature group to unit variance
    # Shape: already L2-normed per subject, but standardize across subjects per grid point
    shape_std = shape_matrix / (shape_matrix.std(axis=0, keepdims=True) + 1e-10)

    # P(resp): standardize
    p_resp_std = (p_resp_vec - p_resp_vec.mean()) / (p_resp_vec.std() + 1e-10)

    # Serovar: standardize
    typhi_std = (typhi_vec - typhi_vec.mean()) / (typhi_vec.std() + 1e-10)

    # Concatenate: shape (8 dims) + P(resp) (1 dim) + serovar (1 dim)
    # Weight: give shape, P(resp), and serovar roughly equal total influence
    # Shape has 8 dims, so downweight each by sqrt(8) to equalize with the scalars
    shape_weight = 1.0 / np.sqrt(len(T_GRID))
    features = np.column_stack([
        shape_std * shape_weight,
        p_resp_std,
        typhi_std,
    ])

    print(f"Feature matrix: {features.shape} (shape×{len(T_GRID)} + P(resp) + serovar)")

    # Hierarchical clustering
    dist_condensed = pdist(features, metric="euclidean")
    Z = linkage(dist_condensed, method="ward")

    for k in [2, 3, 4]:
        labels = fcluster(Z, k, criterion="maxclust")
        sizes = [f"C{c}:{(labels==c).sum()}" for c in range(1, k+1)]
        print(f"  k={k}: {', '.join(sizes)}")

    k_primary = 3
    cluster_labels = fcluster(Z, k_primary, criterion="maxclust")

    # Order clusters by mean P(resp)
    cluster_presp = {}
    for c in range(1, k_primary + 1):
        c_sids = [sids[i] for i in range(n) if cluster_labels[i] == c]
        cluster_presp[c] = np.mean([feature_rows[s]["p_resp"] for s in c_sids])
    rank_order = sorted(cluster_presp, key=lambda c: cluster_presp[c])
    relabel = {old: new + 1 for new, old in enumerate(rank_order)}
    cluster_labels_ord = np.array([relabel[c] for c in cluster_labels])
    cluster_map = dict(zip(sids, cluster_labels_ord))

    cluster_colors = {1: "#2166ac", 2: "#f4a582", 3: "#b2182b"}

    print(f"\nCluster summary (k={k_primary}):")
    for c in range(1, k_primary + 1):
        c_sids = [s for s in sids if cluster_map[s] == c]
        presps = [feature_rows[s]["p_resp"] for s in c_sids]
        serovars = [feature_rows[s]["serovar"] for s in c_sids]
        n_t = serovars.count("typhi")
        n_p = serovars.count("paratyphi")
        print(f"  C{c}: n={len(c_sids)}, mean P(resp)={np.mean(presps):.3f}, "
              f"typhi={n_t}, para={n_p}")

    # =====================================================================
    # Figure 1: Summary 2×2 (dendrogram, mean-centered, absolute, P(resp))
    # =====================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # Dendrogram
    ax = axes[0, 0]
    from scipy.cluster.hierarchy import set_link_color_palette
    set_link_color_palette([cluster_colors[c] for c in range(1, k_primary + 1)])
    dendrogram(Z, ax=ax, color_threshold=Z[-(k_primary-1), 2],
               above_threshold_color="gray", no_labels=True)
    ax.set_xlabel("Subjects"); ax.set_ylabel("Ward distance")
    ax.set_title(f"Hierarchical clustering (k={k_primary})\nFeatures: shape + P(resp) + serovar")
    set_link_color_palette(None)

    # Mean-centered trajectories
    ax = axes[0, 1]
    for sid in sids:
        c = cluster_map[sid]
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        y = np.log10(np.maximum(grp["vi_igg_eu"].values, 1))
        ax.plot(grp["days_since_fever_onset"].values, y - y.mean(),
                color=cluster_colors[c], linewidth=0.8, alpha=0.5)
    ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlim(-10, 400); ax.grid(True, alpha=0.3)
    ax.set_xlabel("Days since fever onset")
    ax.set_ylabel("log10(EU) − subject mean")
    ax.set_title("Mean-centered trajectories by cluster")
    legend_els = [Line2D([0], [0], color=cluster_colors[c], lw=2,
                         label=f"C{c}") for c in range(1, k_primary + 1)]
    ax.legend(handles=legend_els, fontsize=9)

    # Absolute trajectories
    ax = axes[1, 0]
    for sid in sids:
        c = cluster_map[sid]
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                color=cluster_colors[c], linewidth=0.8, alpha=0.5)
    ax.set_yscale("log"); ax.set_ylim(10, 5000); ax.set_xlim(-10, 400)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("Days since fever onset"); ax.set_ylabel("Vi IgG (EU)")
    ax.set_title("Absolute trajectories by cluster")

    # P(resp) + serovar by cluster
    ax = axes[1, 1]
    for c in range(1, k_primary + 1):
        c_sids = [s for s in sids if cluster_map[s] == c]
        presps = [feature_rows[s]["p_resp"] for s in c_sids]
        serovars = [feature_rows[s]["serovar"] for s in c_sids]
        jitter = np.random.RandomState(c).uniform(-0.15, 0.15, len(presps))
        for p, j, sv in zip(presps, jitter, serovars):
            marker = "o" if sv == "typhi" else "s"
            ax.scatter(c + j, p, color=cluster_colors[c], marker=marker,
                       s=30, alpha=0.7, edgecolors="black", linewidths=0.3)
        mean_p = np.mean(presps)
        ax.plot([c - 0.25, c + 0.25], [mean_p, mean_p], color="black",
                linewidth=2.5, zorder=10)
        ax.text(c + 0.28, mean_p, f"{mean_p:.2f}", fontsize=9, va="center",
                fontweight="bold")
    ax.set_xticks(range(1, k_primary + 1))
    ax.set_xticklabels([f"C{c}" for c in range(1, k_primary + 1)])
    ax.set_ylabel("P(responder)"); ax.axhline(0.5, color="gray", linestyle="--", lw=0.8)
    ax.set_title("P(resp) by cluster (— = mean)\n○=Typhi, □=Paratyphi")
    ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle(f"Enriched clustering: shape + P(resp) + serovar (n={n}, k={k_primary})",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "12_enriched_clustering.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("\nSaved 12_enriched_clustering.png")

    # =====================================================================
    # Figure 2: 4×k faceted trajectories (typhi/para × absolute/centered)
    # Colored by P(resp), black LOESS
    # =====================================================================
    presp_cmap = cm.RdYlBu_r
    presp_norm = plt.Normalize(vmin=0.2, vmax=1.0)

    row_defs = [
        ("S. Typhi — absolute", "typhi", "absolute"),
        ("S. Typhi — centered", "typhi", "centered"),
        ("S. Paratyphi A — absolute", "paratyphi", "absolute"),
        ("S. Paratyphi A — centered", "paratyphi", "centered"),
    ]

    fig, axes = plt.subplots(4, k_primary, figsize=(6 * k_primary + 1, 18), sharex=True)

    for col, c in enumerate(range(1, k_primary + 1)):
        c_sids = [s for s in sids if cluster_map[s] == c]

        for row_idx, (row_label, sv_filter, mode) in enumerate(row_defs):
            ax = axes[row_idx, col]
            sv_sids = [s for s in c_sids if feature_rows[s]["serovar"] == sv_filter]

            all_t, all_y = [], []
            for sid in sv_sids:
                grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
                t = grp["days_since_fever_onset"].values
                y_log = np.log10(np.maximum(grp["vi_igg_eu"].values, 1))
                p_r = feature_rows[sid]["p_resp"]
                color = presp_cmap(presp_norm(p_r))

                if mode == "absolute":
                    y_plot = grp["vi_igg_eu"].values
                    y_loess = y_log
                else:
                    y_plot = y_log - feature_rows[sid]["subj_mean_log10"]
                    y_loess = y_plot

                all_t.extend(t); all_y.extend(y_loess)
                ax.plot(t, y_plot, color=color, linewidth=0.8, alpha=0.6)
                ax.scatter(t, y_plot, color=color, s=15, alpha=0.6,
                           edgecolors="black", linewidths=0.2)

            all_t_arr = np.array(all_t)
            all_y_arr = np.array(all_y)
            if len(all_t_arr) >= 6 and len(sv_sids) >= 3:
                smooth = sm_lowess(all_y_arr, all_t_arr, frac=0.6, return_sorted=True)

                # Bootstrap CI by resampling subjects
                t_eval = smooth[:, 0]
                rng = np.random.RandomState(42 + col * 10 + row_idx)
                boot_curves = []
                for _ in range(200):
                    boot_sids = rng.choice(sv_sids, size=len(sv_sids), replace=True)
                    bt, by = [], []
                    for bsid in boot_sids:
                        bgrp = longi[longi["index_id"] == bsid].sort_values("days_since_fever_onset")
                        bt_vals = bgrp["days_since_fever_onset"].values
                        by_log = np.log10(np.maximum(bgrp["vi_igg_eu"].values, 1))
                        if mode == "centered":
                            by_log = by_log - feature_rows[bsid]["subj_mean_log10"]
                        bt.extend(bt_vals); by.extend(by_log)
                    if len(bt) >= 6:
                        bs = sm_lowess(np.array(by), np.array(bt), frac=0.6, return_sorted=True)
                        boot_interp = np.interp(t_eval, bs[:, 0], bs[:, 1])
                        boot_curves.append(boot_interp)
                if len(boot_curves) >= 10:
                    boot_arr = np.array(boot_curves)
                    lo = np.percentile(boot_arr, 2.5, axis=0)
                    hi = np.percentile(boot_arr, 97.5, axis=0)
                    if mode == "absolute":
                        ax.fill_between(t_eval, 10**lo, 10**hi, color="black", alpha=0.2, zorder=9)
                        ax.plot(smooth[:, 0], 10**smooth[:, 1], "k-", lw=3, alpha=0.9, zorder=10)
                    else:
                        ax.fill_between(t_eval, lo, hi, color="black", alpha=0.2, zorder=9)
                        ax.plot(smooth[:, 0], smooth[:, 1], "k-", lw=3, alpha=0.9, zorder=10)
                else:
                    if mode == "absolute":
                        ax.plot(smooth[:, 0], 10**smooth[:, 1], "k-", lw=3, alpha=0.9, zorder=10)
                    else:
                        ax.plot(smooth[:, 0], smooth[:, 1], "k-", lw=3, alpha=0.9, zorder=10)

            ax.set_xlim(-10, 400); ax.grid(True, alpha=0.3)
            if mode == "absolute":
                ax.set_yscale("log"); ax.set_ylim(10, 5000)
            else:
                ax.axhline(0, color="gray", linestyle=":", lw=0.8)

            if row_idx == 0:
                n_t = sum(1 for s in c_sids if feature_rows[s]["serovar"] == "typhi")
                n_p = len(c_sids) - n_t
                ax.set_title(f"C{c} (n={len(c_sids)}: {n_t}T, {n_p}P)\n"
                             f"{sv_filter}: n={len(sv_sids)}", fontsize=10)
            elif row_idx == 2:
                ax.set_title(f"{sv_filter}: n={len(sv_sids)}", fontsize=10)
            if col == 0:
                ax.set_ylabel(f"{row_label}")
            if row_idx == 3:
                ax.set_xlabel("Days since fever onset")

    sm_cb = cm.ScalarMappable(norm=presp_norm, cmap=presp_cmap)
    sm_cb.set_array([])
    fig.colorbar(sm_cb, ax=axes, location="right", shrink=0.4, pad=0.06, label="P(responder)")

    fig.suptitle(f"Enriched clusters × serovar (n={n}, k={k_primary})\n"
                 f"Colored by P(resp), black = LOESS",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "12_enriched_cluster_trajectories.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 12_enriched_cluster_trajectories.png")


if __name__ == "__main__":
    main()
