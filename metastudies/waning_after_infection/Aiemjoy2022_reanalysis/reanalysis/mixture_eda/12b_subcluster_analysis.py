"""
12b: Sub-cluster analysis — finer cuts of the enriched clustering.

Cuts the Ward dendrogram from script 12 at k=8, which gives:
  - Typhi-declining (C1 at k=3) → 2 subclusters
  - Paratyphi (C2 at k=3) → 3 subclusters (the "sneaky" split)
  - Typhi-rising (C3 at k=3) → 3 subclusters

Shows 4×8 faceted trajectories (serovar × absolute/centered)
colored by P(resp) with LOESS + bootstrap CI.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm
from scipy.cluster.hierarchy import linkage, fcluster
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

    # Build features (same as script 12)
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
        y_centered = y_interp - y_interp.mean()
        norm = np.linalg.norm(y_centered)
        if norm < 1e-10:
            continue
        y_normed = y_centered / norm
        shape_full = np.full(len(T_GRID), 0.0)
        shape_full[grid_mask] = y_normed
        p_resp = p_resp_map.get(sid, np.nan)
        if np.isnan(p_resp):
            continue
        is_typhi = 1.0 if serovar_map.get(sid) == "typhi" else 0.0
        feature_rows[sid] = {
            "shape": shape_full, "p_resp": p_resp, "is_typhi": is_typhi,
            "serovar": serovar_map.get(sid, "unknown"),
            "subj_mean_log10": y.mean(),
        }

    sids = sorted(feature_rows.keys())
    n = len(sids)

    shape_matrix = np.array([feature_rows[s]["shape"] for s in sids])
    p_resp_vec = np.array([feature_rows[s]["p_resp"] for s in sids])
    typhi_vec = np.array([feature_rows[s]["is_typhi"] for s in sids])

    shape_std = shape_matrix / (shape_matrix.std(axis=0, keepdims=True) + 1e-10)
    p_resp_std = (p_resp_vec - p_resp_vec.mean()) / (p_resp_vec.std() + 1e-10)
    typhi_std = (typhi_vec - typhi_vec.mean()) / (typhi_vec.std() + 1e-10)
    sw = 1.0 / np.sqrt(len(T_GRID))
    features = np.column_stack([shape_std * sw, p_resp_std, typhi_std])

    dist_condensed = pdist(features, metric="euclidean")
    Z = linkage(dist_condensed, method="ward")

    # Cut at k=8 subclusters
    k = 8
    raw_labels = fcluster(Z, k, criterion="maxclust")

    # Order by (serovar majority, then mean P(resp))
    cluster_info = {}
    for c in range(1, k + 1):
        ci = [i for i in range(n) if raw_labels[i] == c]
        c_sids = [sids[i] for i in ci]
        n_t = sum(1 for s in c_sids if feature_rows[s]["serovar"] == "typhi")
        n_p = len(c_sids) - n_t
        mp = np.mean([feature_rows[s]["p_resp"] for s in c_sids])
        is_typhi_majority = n_t > n_p
        cluster_info[c] = {"n": len(ci), "n_t": n_t, "n_p": n_p,
                           "mean_presp": mp, "typhi_maj": is_typhi_majority}

    # Sort: typhi-low-P first, then paratyphi, then typhi-high-P
    def sort_key(c):
        info = cluster_info[c]
        if info["typhi_maj"] and info["mean_presp"] < 0.6:
            return (0, info["mean_presp"])
        elif not info["typhi_maj"]:
            return (1, info["mean_presp"])
        else:
            return (2, info["mean_presp"])

    ordered = sorted(cluster_info.keys(), key=sort_key)
    relabel = {old: new + 1 for new, old in enumerate(ordered)}
    labels = np.array([relabel[c] for c in raw_labels])
    cluster_map = dict(zip(sids, labels))

    # Print summary
    print(f"Sub-cluster analysis: n={n}, k={k}")
    parent_names = {}
    for c in range(1, k + 1):
        c_sids = [s for s in sids if cluster_map[s] == c]
        info_orig = cluster_info[ordered[c - 1]]
        n_t, n_p = info_orig["n_t"], info_orig["n_p"]
        mp = info_orig["mean_presp"]
        sv = "Typhi" if info_orig["typhi_maj"] else "Paratyphi"
        parent = "declining" if sv == "Typhi" and mp < 0.6 else (
            "Paratyphi" if sv == "Paratyphi" else "rising")
        parent_names[c] = parent
        print(f"  SC{c} ({parent}): n={len(c_sids)}, {n_t}T/{n_p}P, mean P(resp)={mp:.2f}")

    # =====================================================================
    # Figure: 2 rows × k cols — absolute + centered, colored by P(resp)
    # =====================================================================
    presp_cmap = cm.RdYlBu_r
    presp_norm = plt.Normalize(vmin=0.2, vmax=1.0)

    fig, axes = plt.subplots(2, k, figsize=(3.5 * k, 8), sharex=True)

    for col, c in enumerate(range(1, k + 1)):
        c_sids = [s for s in sids if cluster_map[s] == c]

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
                info = cluster_info[ordered[col]]
                sv_label = f"{info['n_t']}T,{info['n_p']}P"
                ax.set_title(f"SC{c} ({parent_names[c]})\n"
                             f"n={len(c_sids)} ({sv_label})\n"
                             f"P̄={info['mean_presp']:.2f}", fontsize=8)
            if col == 0:
                label = "Vi IgG (EU)" if mode == "absolute" else "log10 − subj mean"
                ax.set_ylabel(label)
            if row_idx == 1:
                ax.set_xlabel("Days")

    sm_cb = cm.ScalarMappable(norm=presp_norm, cmap=presp_cmap)
    sm_cb.set_array([])
    fig.colorbar(sm_cb, ax=axes, location="right", shrink=0.5, pad=0.04, label="P(resp)")

    fig.suptitle(f"Sub-cluster analysis (k={k}, n={n})\n"
                 f"Ordered: Typhi-declining → Paratyphi → Typhi-rising",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "12b_subclusters.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved 12b_subclusters.png")


if __name__ == "__main__":
    main()
