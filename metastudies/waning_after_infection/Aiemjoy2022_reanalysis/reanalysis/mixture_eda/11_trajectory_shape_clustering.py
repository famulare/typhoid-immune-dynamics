"""
11: Unsupervised trajectory shape clustering.

For subjects with ≥3 points and ≥2 within 200 days:
  1. Linear-interpolate log10(EU) onto a common time grid
  2. Mean-center (remove level) + L2-normalize (unit shape vectors)
  3. Pairwise cosine similarity (= correlation for mean-centered data)
  4. Hierarchical clustering (Ward's method)
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

    # Hierarchical clustering (Ward's method)
    dist_condensed = squareform(dist_matrix)
    Z = linkage(dist_condensed, method="ward")

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
    cluster_names = {1: "Low P(resp)", 2: "Medium", 3: "High P(resp)"}

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

    # Panel D: P(resp) by cluster + serovar breakdown
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
    ax.set_xticks(range(1, k_primary + 1))
    ax.set_xticklabels([f"C{c}\n{cluster_names[c]}" for c in range(1, k_primary + 1)])
    ax.set_ylabel("P(responder) — Step 1 mixture")
    ax.set_title("P(resp) distribution by shape cluster\n(○=Typhi, □=Paratyphi)")
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
    ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle(f"Trajectory shape clustering (n={n}, k={k_primary})\n"
                 f"Linear interp → mean-center → L2 norm → cosine distance → Ward",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "11_shape_clustering.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("\nSaved 11_shape_clustering.png")


if __name__ == "__main__":
    main()
