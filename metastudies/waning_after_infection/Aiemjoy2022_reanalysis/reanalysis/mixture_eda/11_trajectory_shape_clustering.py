"""
11: Unsupervised trajectory shape clustering.

For subjects with ≥3 points and ≥2 within 200 days:
  1. Linear-interpolate log10(EU) onto a common time grid
  2. Mean-center (remove level) + L2-normalize (unit shape vectors)
  3. Pairwise cosine similarity (= correlation for mean-centered data)
  4. Hierarchical clustering (average linkage — Ward requires Euclidean
     geometry and is invalid on a precomputed cosine-distance matrix)
  5. Visualize: dendrogram, clustered trajectories, cluster vs P(resp)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from importlib import import_module
m08 = import_module("08_fold_change_mixture")

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis/mixture_eda")

# Common time grid for interpolation (days since fever onset)
T_GRID = np.array([5, 15, 30, 60, 90, 120, 150, 180])


def main():
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    raw = pd.read_csv(DATA_DIR / "TypoidCaseData_github_09.30.21.csv")

    # Load fold-change data and fit working model for P(resp)
    df_fc = m08.load_fold_change_data()
    x = df_fc["log2_fc"].values
    p1 = m08.fit_two_gaussian(x)
    p_resp_map = dict(zip(df_fc["index_id"], 1 - p1["p_noise"]))

    # Select subjects: ≥3 observations, ≥2 within 200 days
    longi = long[long["has_longitudinal"]]
    nobs = longi.groupby("index_id").size()
    multi_ids = set(nobs[nobs >= 3].index)
    good_ids = set()
    for sid in multi_ids:
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        if (grp["days_since_fever_onset"] <= 200).sum() >= 2:
            good_ids.add(sid)

    print(f"Subjects with ≥3 pts, ≥2 by day 200: {len(good_ids)}")

    # Interpolate onto common grid
    shape_vectors = {}
    subject_info = {}
    serovar_map = raw.set_index("index_id")["bldculres"]

    for sid in sorted(good_ids):
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        t = grp["days_since_fever_onset"].values
        y = np.log10(np.maximum(grp["vi_igg_eu"].values, 1))

        # Interpolate — only within observed range (no extrapolation)
        t_min, t_max = t.min(), t.max()
        grid_mask = (T_GRID >= t_min) & (T_GRID <= t_max)
        if grid_mask.sum() < 3:
            continue  # need at least 3 grid points covered

        t_use = T_GRID[grid_mask]
        y_interp = np.interp(t_use, t, y)

        # Mean-center
        y_centered = y_interp - y_interp.mean()

        # L2 normalize
        norm = np.linalg.norm(y_centered)
        if norm < 1e-10:
            continue  # flat trajectory, skip
        y_normed = y_centered / norm

        # Store with the grid indices used
        shape_vectors[sid] = {"t_grid": t_use, "y_normed": y_normed,
                              "y_centered": y_centered, "grid_mask": grid_mask}
        subject_info[sid] = {
            "serovar": serovar_map.get(sid, "unknown"),
            "p_resp": p_resp_map.get(sid, np.nan),
            "n_obs": len(grp),
        }

    sids = sorted(shape_vectors.keys())
    n = len(sids)
    print(f"Subjects with ≥3 grid points covered: {n}")

    # Pairwise cosine similarity on the FULL grid
    # For subjects covering different grid subsets, compute similarity
    # only on overlapping grid points
    sim_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            mask_i = shape_vectors[sids[i]]["grid_mask"]
            mask_j = shape_vectors[sids[j]]["grid_mask"]
            overlap = mask_i & mask_j
            if overlap.sum() < 3:
                sim_matrix[i, j] = 0
                continue
            # Re-interpolate both onto overlap grid
            t_overlap = T_GRID[overlap]
            yi = np.interp(t_overlap, shape_vectors[sids[i]]["t_grid"],
                           shape_vectors[sids[i]]["y_normed"])
            yj = np.interp(t_overlap, shape_vectors[sids[j]]["t_grid"],
                           shape_vectors[sids[j]]["y_normed"])
            # Re-normalize on overlap
            ni, nj = np.linalg.norm(yi), np.linalg.norm(yj)
            if ni < 1e-10 or nj < 1e-10:
                sim_matrix[i, j] = 0
            else:
                sim_matrix[i, j] = np.dot(yi / ni, yj / nj)

    # Convert to distance: d = 1 - sim (cosine distance)
    dist_matrix = 1 - sim_matrix
    np.fill_diagonal(dist_matrix, 0)
    # Ensure symmetry and non-negative
    dist_matrix = (dist_matrix + dist_matrix.T) / 2
    dist_matrix = np.maximum(dist_matrix, 0)

    # Hierarchical clustering on precomputed cosine distances.
    # Use average (UPGMA) linkage — Ward minimizes Euclidean SSE and is
    # not well-defined on arbitrary distance matrices.
    dist_condensed = squareform(dist_matrix)
    Z = linkage(dist_condensed, method="average")

    # Cut at k clusters
    for k in [2, 3, 4]:
        labels = fcluster(Z, k, criterion="maxclust")
        sizes = [f"C{c}:{(labels==c).sum()}" for c in range(1, k+1)]
        print(f"  k={k}: {', '.join(sizes)}")

    # Use k=3 as primary (interpretable: rise-fall, flat, decline)
    k_primary = 3
    cluster_labels = fcluster(Z, k_primary, criterion="maxclust")
    cluster_map = dict(zip(sids, cluster_labels))

    # Order clusters by mean P(resp) for interpretability
    cluster_presp = {}
    for c in range(1, k_primary + 1):
        c_sids = [s for s in sids if cluster_map[s] == c]
        cluster_presp[c] = np.mean([subject_info[s]["p_resp"] for s in c_sids
                                    if not np.isnan(subject_info[s]["p_resp"])])
    rank_order = sorted(cluster_presp, key=lambda c: cluster_presp[c])
    relabel = {old: new + 1 for new, old in enumerate(rank_order)}
    cluster_labels_ordered = np.array([relabel[c] for c in cluster_labels])
    cluster_map_ordered = dict(zip(sids, cluster_labels_ordered))

    cluster_colors = {1: "#2166ac", 2: "#f4a582", 3: "#b2182b"}
    cluster_names = {1: "Declining", 2: "Mixed/flat", 3: "Rising"}

    # Print cluster summary
    print(f"\nCluster summary (k={k_primary}, ordered by mean P(resp)):")
    for c in range(1, k_primary + 1):
        c_sids = [s for s in sids if cluster_map_ordered[s] == c]
        presps = [subject_info[s]["p_resp"] for s in c_sids if not np.isnan(subject_info[s]["p_resp"])]
        serovars = [subject_info[s]["serovar"] for s in c_sids]
        n_typhi = serovars.count("typhi")
        n_para = serovars.count("paratyphi")
        print(f"  C{c} ({cluster_names[c]}): n={len(c_sids)}, "
              f"mean P(resp)={np.mean(presps):.3f}, "
              f"typhi={n_typhi}, para={n_para}")

    # =====================================================================
    # Figure: 2×2
    # =====================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # Panel A: Dendrogram
    ax = axes[0, 0]
    # Color the dendrogram by cluster
    from scipy.cluster.hierarchy import set_link_color_palette
    set_link_color_palette([cluster_colors[c] for c in range(1, k_primary + 1)])
    dendrogram(Z, ax=ax, color_threshold=Z[-(k_primary-1), 2],
               above_threshold_color="gray", no_labels=True)
    ax.set_xlabel("Subjects")
    ax.set_ylabel("Ward distance")
    ax.set_title(f"Hierarchical clustering (k={k_primary})")
    set_link_color_palette(None)

    # Panel B: Mean-centered trajectories by cluster (one panel, overlaid)
    ax = axes[0, 1]
    for sid in sids:
        c = cluster_map_ordered[sid]
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        t = grp["days_since_fever_onset"].values
        y = np.log10(np.maximum(grp["vi_igg_eu"].values, 1))
        y_centered = y - y.mean()
        ax.plot(t, y_centered, color=cluster_colors[c], linewidth=0.8, alpha=0.5)
        ax.scatter(t, y_centered, color=cluster_colors[c], s=12, alpha=0.5,
                   edgecolors="none")
    ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlim(-10, 400); ax.grid(True, alpha=0.3)
    ax.set_xlabel("Days since fever onset")
    ax.set_ylabel("log10(EU) − subject mean")
    ax.set_title("Mean-centered trajectories by cluster")
    legend_els = [Line2D([0], [0], color=cluster_colors[c], lw=2,
                         label=f"C{c}: {cluster_names[c]}") for c in range(1, k_primary + 1)]
    ax.legend(handles=legend_els, fontsize=9)

    # Panel C: Absolute trajectories faceted by cluster
    ax = axes[1, 0]
    for sid in sids:
        c = cluster_map_ordered[sid]
        grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
        ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                color=cluster_colors[c], linewidth=0.8, alpha=0.5)
        ax.scatter(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                   color=cluster_colors[c], s=12, alpha=0.5, edgecolors="none")
    ax.set_yscale("log"); ax.set_ylim(10, 5000); ax.set_xlim(-10, 400)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("Days since fever onset")
    ax.set_ylabel("Vi IgG (EU)")
    ax.set_title("Absolute trajectories by cluster")

    # Panel D: P(resp) by cluster — individual points + mean
    ax = axes[1, 1]
    for c in range(1, k_primary + 1):
        c_sids = [s for s in sids if cluster_map_ordered[s] == c]
        presps = [subject_info[s]["p_resp"] for s in c_sids
                  if not np.isnan(subject_info[s]["p_resp"])]
        jitter = np.random.RandomState(c).uniform(-0.15, 0.15, len(presps))
        serovars = [subject_info[s]["serovar"] for s in c_sids
                    if not np.isnan(subject_info[s]["p_resp"])]
        for p, j, sv in zip(presps, jitter, serovars):
            marker = "o" if sv == "typhi" else "s"
            ax.scatter(c + j, p, color=cluster_colors[c], marker=marker,
                       s=30, alpha=0.7, edgecolors="black", linewidths=0.3)
        # Mean P(resp) per cluster
        mean_p = np.mean(presps)
        ax.plot([c - 0.25, c + 0.25], [mean_p, mean_p], color="black",
                linewidth=2.5, zorder=10)
        ax.text(c + 0.28, mean_p, f"{mean_p:.2f}", fontsize=9, va="center",
                fontweight="bold")
    ax.set_xticks(range(1, k_primary + 1))
    ax.set_xticklabels([f"C{c}\n{cluster_names[c]}" for c in range(1, k_primary + 1)])
    ax.set_ylabel("P(responder) — Step 1 mixture")
    ax.set_title("P(resp) by shape cluster (— = mean)\n(○=Typhi, □=Paratyphi)")
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
    ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle(f"Trajectory shape clustering (n={n}, k={k_primary})\n"
                 f"Linear interp → mean-center → L2 norm → cosine distance → Ward",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "11_shape_clustering.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("\nSaved 11_shape_clustering.png")


    # =====================================================================
    # Figure 2: 4×3 faceted trajectories by cluster × serovar
    # Rows: typhi absolute, typhi centered, paratyphi absolute, para centered
    # Colored by P(resp). LOESS per panel. Filter outlier (indx387).
    # =====================================================================
    from statsmodels.nonparametric.smoothers_lowess import lowess as sm_lowess
    import matplotlib.cm as cm

    # P(resp) colormap
    presp_cmap = cm.RdYlBu_r
    presp_norm = plt.Normalize(vmin=0.2, vmax=1.0)

    # Filter obvious outlier (indx387: collapses to ~0 EU)
    outlier_ids = {"indx387"}

    row_defs = [
        ("S. Typhi — absolute", "typhi", "absolute"),
        ("S. Typhi — mean-centered", "typhi", "centered"),
        ("S. Paratyphi A — absolute", "paratyphi", "absolute"),
        ("S. Paratyphi A — mean-centered", "paratyphi", "centered"),
    ]

    fig, axes = plt.subplots(4, k_primary, figsize=(6 * k_primary, 18),
                             sharex=True)

    # Precompute subject means
    subj_means = {}
    for sid in sids:
        grp = longi[longi["index_id"] == sid]
        subj_means[sid] = np.log10(np.maximum(grp["vi_igg_eu"].values, 1)).mean()

    for col, c in enumerate(range(1, k_primary + 1)):
        c_sids = [s for s in sids if cluster_map_ordered[s] == c and s not in outlier_ids]

        for row_idx, (row_label, sv_filter, mode) in enumerate(row_defs):
            ax = axes[row_idx, col]
            sv_sids = [s for s in c_sids if subject_info[s]["serovar"] == sv_filter]

            # Collect for LOESS
            all_t, all_y = [], []

            for sid in sv_sids:
                grp = longi[longi["index_id"] == sid].sort_values("days_since_fever_onset")
                t = grp["days_since_fever_onset"].values
                y_log = np.log10(np.maximum(grp["vi_igg_eu"].values, 1))
                p_r = subject_info[sid]["p_resp"]
                color = presp_cmap(presp_norm(p_r))

                if mode == "absolute":
                    y_plot = grp["vi_igg_eu"].values
                    y_for_loess = y_log
                else:
                    y_plot = y_log - subj_means[sid]
                    y_for_loess = y_plot

                all_t.extend(t)
                all_y.extend(y_for_loess)

                if mode == "absolute":
                    ax.plot(t, y_plot, color=color, linewidth=0.8, alpha=0.6)
                    ax.scatter(t, y_plot, color=color, s=15, alpha=0.6,
                               edgecolors="black", linewidths=0.2)
                else:
                    ax.plot(t, y_plot, color=color, linewidth=0.8, alpha=0.6)
                    ax.scatter(t, y_plot, color=color, s=15, alpha=0.6,
                               edgecolors="black", linewidths=0.2)

            # LOESS
            all_t = np.array(all_t)
            all_y = np.array(all_y)
            if len(all_t) >= 6:
                smooth = sm_lowess(all_y, all_t, frac=0.6, return_sorted=True)
                if mode == "absolute":
                    ax.plot(smooth[:, 0], 10**smooth[:, 1], color="black",
                            linewidth=3, alpha=0.9, zorder=10)
                else:
                    ax.plot(smooth[:, 0], smooth[:, 1], color="black",
                            linewidth=3, alpha=0.9, zorder=10)

            ax.set_xlim(-10, 400); ax.grid(True, alpha=0.3)
            if mode == "absolute":
                ax.set_yscale("log"); ax.set_ylim(10, 5000)
            else:
                ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)

            if row_idx == 0:
                n_sv = len(sv_sids)
                ax.set_title(f"C{c}: {cluster_names[c]} (n={len(c_sids)})\n"
                             f"{sv_filter}: n={n_sv}", fontsize=10)
            elif row_idx == 2:
                n_sv = len(sv_sids)
                ax.set_title(f"{sv_filter}: n={n_sv}", fontsize=10)

            if col == 0:
                if mode == "absolute":
                    ax.set_ylabel(f"{row_label}\nVi IgG (EU)")
                else:
                    ax.set_ylabel(f"{row_label}\nlog10 − subj mean")
            if row_idx == 3:
                ax.set_xlabel("Days since fever onset")

    # Colorbar
    sm_cb = cm.ScalarMappable(norm=presp_norm, cmap=presp_cmap)
    sm_cb.set_array([])
    cbar = fig.colorbar(sm_cb, ax=axes, location="right", shrink=0.5, pad=0.02)
    cbar.set_label("P(responder)")

    fig.suptitle(f"Trajectory shape clusters × serovar (n={n}, k={k_primary})\n"
                 f"Colored by P(resp), black = LOESS",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "11_cluster_trajectories.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 11_cluster_trajectories.png")

    # =====================================================================
    # Figure 3: Cluster membership vs P(resp)
    # Stacked bar / proportion plot showing cluster composition at each
    # P(resp) level, and P(resp) distribution colored by cluster
    # =====================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: P(resp) histogram colored by cluster
    ax = axes[0]
    presp_by_cluster = {}
    for c in range(1, k_primary + 1):
        c_sids = [s for s in sids if cluster_map_ordered[s] == c]
        presp_by_cluster[c] = [subject_info[s]["p_resp"] for s in c_sids
                               if not np.isnan(subject_info[s]["p_resp"])]
    bins = np.linspace(0, 1, 16)
    bottom = np.zeros(len(bins) - 1)
    for c in range(1, k_primary + 1):
        counts, _ = np.histogram(presp_by_cluster[c], bins=bins)
        ax.bar(bins[:-1] + (bins[1] - bins[0]) / 2, counts, width=bins[1] - bins[0],
               bottom=bottom, color=cluster_colors[c], alpha=0.8,
               edgecolor="black", linewidth=0.3, label=f"C{c}: {cluster_names[c]}")
        bottom += counts
    ax.set_xlabel("P(responder)")
    ax.set_ylabel("Count")
    ax.set_title("P(resp) distribution by shape cluster")
    ax.legend(fontsize=9)
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=0.8)

    # Panel B: Cluster proportion vs P(resp) bin
    ax = axes[1]
    bin_centers = bins[:-1] + (bins[1] - bins[0]) / 2
    totals = np.zeros(len(bins) - 1)
    counts_by_c = {}
    for c in range(1, k_primary + 1):
        counts_by_c[c], _ = np.histogram(presp_by_cluster[c], bins=bins)
        totals += counts_by_c[c]
    bottom = np.zeros(len(bins) - 1)
    for c in range(1, k_primary + 1):
        fracs = np.where(totals > 0, counts_by_c[c] / totals, 0)
        ax.bar(bin_centers, fracs, width=bins[1] - bins[0], bottom=bottom,
               color=cluster_colors[c], alpha=0.8, edgecolor="black", linewidth=0.3,
               label=f"C{c}: {cluster_names[c]}")
        bottom += fracs
    ax.set_xlabel("P(responder)")
    ax.set_ylabel("Cluster proportion")
    ax.set_title("Shape cluster composition vs P(resp)")
    ax.legend(fontsize=9)
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=0.8)
    ax.set_ylim(0, 1)

    fig.suptitle(f"Shape cluster membership vs mixture P(resp) (n={n})",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "11_cluster_vs_presp.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved 11_cluster_vs_presp.png")


if __name__ == "__main__":
    main()
