"""
12c: Clustering on raw (non-centered) log10 titer traces + covariates.

Unlike 12 (which mean-centered to cluster on shape only), this uses the
raw interpolated log10(EU) values — so titer LEVEL contributes to the
clustering alongside trajectory shape, P(resp), and serovar.

Feature vector per subject (all standardized):
  - Raw log10(EU) on common time grid (not centered, not L2-normed)
  - P(resp): scalar
  - is_typhi: binary indicator

Distance: Euclidean on standardized features.
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
m08 = __import__("08_fold_change_mixture")

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/mixture_eda")

T_GRID = np.array([5, 15, 30, 60, 90, 120, 150, 180])
OUTLIER_IDS = {"indx387"}


def main():
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    raw = pd.read_csv(DATA_DIR / "TypoidCaseData_github_09.30.21.csv")
    serovar_map = raw.set_index("index_id")["bldculres"]

    df_fc = m08.load_fold_change_data()
    p1 = m08.fit_two_gaussian(df_fc["log2_fc"].values)
    p_resp_map = dict(zip(df_fc["index_id"], 1 - p1["p_noise"]))

    longi = long[long["has_longitudinal"]]
    nobs = longi.groupby("index_id").size()
    multi_ids = set(nobs[nobs >= 3].index) - OUTLIER_IDS
    good_ids = set()
    for sid in multi_ids:
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        if (grp["days_since_fever_onset"] <= 200).sum() >= 2:
            good_ids.add(sid)

    # Build features: RAW log10(EU) on grid (not centered) + P(resp) + serovar
    feature_rows = {}
    for sid in sorted(good_ids):
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        t = grp["days_since_fever_onset"].values
        y = np.log10(np.maximum(grp["vi_igg_eu"].values, 1))

        t_min, t_max = t.min(), t.max()
        grid_mask = (T_GRID >= t_min) & (T_GRID <= t_max)
        if grid_mask.sum() < 3:
            continue

        y_interp = np.interp(T_GRID[grid_mask], t, y)

        # Pad uncovered grid points with NaN → fill with 0 after standardization
        raw_full = np.full(len(T_GRID), np.nan)
        raw_full[grid_mask] = y_interp

        p_resp = p_resp_map.get(sid, np.nan)
        if np.isnan(p_resp):
            continue
        is_typhi = 1.0 if serovar_map.get(sid) == "typhi" else 0.0

        feature_rows[sid] = {
            "raw_trace": raw_full,
            "grid_mask": grid_mask,
            "p_resp": p_resp,
            "is_typhi": is_typhi,
            "serovar": serovar_map.get(sid, "unknown"),
            "subj_mean_log10": y.mean(),
        }

    sids = sorted(feature_rows.keys())
    n = len(sids)
    print(f"Subjects: {n}")

    raw_matrix = np.array([feature_rows[s]["raw_trace"] for s in sids])
    raw_matrix = np.nan_to_num(raw_matrix, nan=0.0)
    p_resp_vec = np.array([feature_rows[s]["p_resp"] for s in sids])
    typhi_vec = np.array([feature_rows[s]["is_typhi"] for s in sids])

    # Standardize each feature group
    raw_std = (raw_matrix - raw_matrix.mean(axis=0, keepdims=True)) / (
        raw_matrix.std(axis=0, keepdims=True) + 1e-10)
    p_resp_std = (p_resp_vec - p_resp_vec.mean()) / (p_resp_vec.std() + 1e-10)
    typhi_std = (typhi_vec - typhi_vec.mean()) / (typhi_vec.std() + 1e-10)

    # Weight: raw trace (8d) downweighted per dim to equalize with scalars
    trace_weight = 1.0 / np.sqrt(len(T_GRID))
    features = np.column_stack([raw_std * trace_weight, p_resp_std, typhi_std])
    print(f"Features: {features.shape}")

    dist_condensed = pdist(features, metric="euclidean")
    Z = linkage(dist_condensed, method="ward")

    for k in [2, 3, 4, 5, 6, 8]:
        labels = fcluster(Z, k, criterion="maxclust")
        parts = []
        for c in range(1, k + 1):
            ci = [i for i in range(n) if labels[i] == c]
            c_sids = [sids[i] for i in ci]
            n_t = sum(1 for s in c_sids if feature_rows[s]["serovar"] == "typhi")
            mp = np.mean([feature_rows[s]["p_resp"] for s in c_sids])
            parts.append(f"{len(ci)}({n_t}T,{len(ci)-n_t}P,P̄={mp:.2f})")
        print(f"  k={k}: {' | '.join(parts)}")

    k_primary = 4
    cluster_labels = fcluster(Z, k_primary, criterion="maxclust")

    # Order by (serovar, then mean titer level)
    cluster_info = {}
    for c in range(1, k_primary + 1):
        ci = [i for i in range(n) if cluster_labels[i] == c]
        c_sids = [sids[i] for i in ci]
        n_t = sum(1 for s in c_sids if feature_rows[s]["serovar"] == "typhi")
        mp = np.mean([feature_rows[s]["p_resp"] for s in c_sids])
        mean_eu = np.mean([10**feature_rows[s]["subj_mean_log10"] for s in c_sids])
        cluster_info[c] = {"n": len(ci), "n_t": n_t, "n_p": len(ci) - n_t,
                           "mean_presp": mp, "mean_eu": mean_eu,
                           "typhi_maj": n_t > len(ci) - n_t}

    def sort_key(c):
        info = cluster_info[c]
        return (0 if info["typhi_maj"] else 1, info["mean_presp"])

    ordered = sorted(cluster_info.keys(), key=sort_key)
    relabel = {old: new + 1 for new, old in enumerate(ordered)}
    labels_ord = np.array([relabel[c] for c in cluster_labels])
    cluster_map = dict(zip(sids, labels_ord))

    cluster_colors = {1: "#2166ac", 2: "#92c5de", 3: "#f4a582", 4: "#b2182b"}
    if k_primary > 4:
        cmap = plt.cm.tab10
        cluster_colors = {c: cmap(c / k_primary) for c in range(1, k_primary + 1)}

    print(f"\nCluster summary (k={k_primary}, raw titers + P(resp) + serovar):")
    for c in range(1, k_primary + 1):
        info = cluster_info[ordered[c - 1]]
        c_sids = [s for s in sids if cluster_map[s] == c]
        print(f"  C{c}: n={len(c_sids)}, {info['n_t']}T/{info['n_p']}P, "
              f"mean P(resp)={info['mean_presp']:.2f}, mean EU={info['mean_eu']:.0f}")

    # =====================================================================
    # Figure 1: Summary 2×2
    # =====================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # Dendrogram
    ax = axes[0, 0]
    from scipy.cluster.hierarchy import set_link_color_palette
    set_link_color_palette([cluster_colors[c] for c in range(1, k_primary + 1)])
    dendrogram(Z, ax=ax, color_threshold=Z[-(k_primary - 1), 2],
               above_threshold_color="gray", no_labels=True)
    ax.set_xlabel("Subjects"); ax.set_ylabel("Ward distance")
    ax.set_title(f"Hierarchical clustering (k={k_primary})\nFeatures: raw log10(EU) + P(resp) + serovar")
    set_link_color_palette(None)

    # Absolute trajectories
    ax = axes[0, 1]
    for sid in sids:
        c = cluster_map[sid]
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                color=cluster_colors[c], linewidth=0.8, alpha=0.5)
    ax.set_yscale("log"); ax.set_ylim(10, 5000); ax.set_xlim(-10, 400)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("Days since fever onset"); ax.set_ylabel("Vi IgG (EU)")
    ax.set_title("Absolute trajectories by cluster")
    legend_els = [Line2D([0], [0], color=cluster_colors[c], lw=2,
                         label=f"C{c}") for c in range(1, k_primary + 1)]
    ax.legend(handles=legend_els, fontsize=9)

    # Mean-centered
    ax = axes[1, 0]
    for sid in sids:
        c = cluster_map[sid]
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        y = np.log10(np.maximum(grp["vi_igg_eu"].values, 1))
        ax.plot(grp["days_since_fever_onset"].values, y - y.mean(),
                color=cluster_colors[c], linewidth=0.8, alpha=0.5)
    ax.axhline(0, color="gray", linestyle=":", lw=0.8)
    ax.set_xlim(-10, 400); ax.grid(True, alpha=0.3)
    ax.set_xlabel("Days since fever onset")
    ax.set_ylabel("log10(EU) − subject mean")
    ax.set_title("Mean-centered trajectories by cluster")

    # P(resp) by cluster
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
        ax.plot([c - 0.25, c + 0.25], [mean_p, mean_p], "k-", lw=2.5, zorder=10)
        ax.text(c + 0.28, mean_p, f"{mean_p:.2f}", fontsize=9, va="center", fontweight="bold")
    ax.set_xticks(range(1, k_primary + 1))
    ax.set_xticklabels([f"C{c}" for c in range(1, k_primary + 1)])
    ax.set_ylabel("P(responder)"); ax.axhline(0.5, color="gray", linestyle="--", lw=0.8)
    ax.set_title("P(resp) by cluster (— = mean)\n○=Typhi, □=Paratyphi")
    ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle(f"Raw-titer clustering (n={n}, k={k_primary})\nFeatures: log10(EU) trace + P(resp) + serovar",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "12c_raw_titer_clustering.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("\nSaved 12c_raw_titer_clustering.png")

    # =====================================================================
    # Figure 2: 2×k faceted (absolute + centered), colored by P(resp), LOESS+CI
    # =====================================================================
    presp_cmap = cm.RdYlBu_r
    presp_norm = plt.Normalize(vmin=0.2, vmax=1.0)

    fig, axes = plt.subplots(2, k_primary, figsize=(5 * k_primary, 10), sharex=True)

    for col, c in enumerate(range(1, k_primary + 1)):
        c_sids = [s for s in sids if cluster_map[s] == c]
        info = cluster_info[ordered[col]]

        for row_idx, mode in enumerate(["absolute", "centered"]):
            ax = axes[row_idx, col]
            all_t, all_y = [], []

            for sid in c_sids:
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
                ax.scatter(t, y_plot, color=color, s=12, alpha=0.6,
                           edgecolors="black", linewidths=0.2)

            # LOESS + bootstrap CI
            all_t_arr, all_y_arr = np.array(all_t), np.array(all_y)
            if len(all_t_arr) >= 6 and len(c_sids) >= 3:
                smooth = sm_lowess(all_y_arr, all_t_arr, frac=0.6, return_sorted=True)
                t_eval = smooth[:, 0]
                rng = np.random.RandomState(42 + col * 10 + row_idx)
                boot_curves = []
                for _ in range(200):
                    bs = rng.choice(c_sids, size=len(c_sids), replace=True)
                    bt, by = [], []
                    for bsid in bs:
                        bg = longi[longi["index_id"] == bsid].sort_values("days_since_fever_onset")
                        btv = bg["days_since_fever_onset"].values
                        byv = np.log10(np.maximum(bg["vi_igg_eu"].values, 1))
                        if mode == "centered":
                            byv = byv - feature_rows[bsid]["subj_mean_log10"]
                        bt.extend(btv); by.extend(byv)
                    if len(bt) >= 6:
                        bsm = sm_lowess(np.array(by), np.array(bt), frac=0.6, return_sorted=True)
                        boot_curves.append(np.interp(t_eval, bsm[:, 0], bsm[:, 1]))
                if len(boot_curves) >= 10:
                    ba = np.array(boot_curves)
                    lo, hi = np.percentile(ba, 2.5, axis=0), np.percentile(ba, 97.5, axis=0)
                    if mode == "absolute":
                        ax.fill_between(t_eval, 10**lo, 10**hi, color="black", alpha=0.15)
                        ax.plot(smooth[:, 0], 10**smooth[:, 1], "k-", lw=2.5, zorder=10)
                    else:
                        ax.fill_between(t_eval, lo, hi, color="black", alpha=0.15)
                        ax.plot(smooth[:, 0], smooth[:, 1], "k-", lw=2.5, zorder=10)

            ax.set_xlim(-10, 400); ax.grid(True, alpha=0.3)
            if mode == "absolute":
                ax.set_yscale("log"); ax.set_ylim(10, 5000)
            else:
                ax.axhline(0, color="gray", linestyle=":", lw=0.8)

            if row_idx == 0:
                ax.set_title(f"C{c}: n={len(c_sids)} ({info['n_t']}T, {info['n_p']}P)\n"
                             f"P̄={info['mean_presp']:.2f}, mean EU={info['mean_eu']:.0f}",
                             fontsize=10)
            if col == 0:
                label = "Vi IgG (EU)" if mode == "absolute" else "log10 − subj mean"
                ax.set_ylabel(label)
            if row_idx == 1:
                ax.set_xlabel("Days since fever onset")

    fig.suptitle(f"Raw-titer clusters (n={n}, k={k_primary})\n"
                 f"Colored by P(resp), black = LOESS ± 95% CI",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.93, 0.95])
    cbar_ax = fig.add_axes([0.94, 0.25, 0.015, 0.5])
    sm_cb = cm.ScalarMappable(norm=presp_norm, cmap=presp_cmap)
    sm_cb.set_array([])
    fig.colorbar(sm_cb, cax=cbar_ax, label="P(resp)")
    fig.savefig(OUT_DIR / "12c_raw_titer_trajectories.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 12c_raw_titer_trajectories.png")


if __name__ == "__main__":
    main()
