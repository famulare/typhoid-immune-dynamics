"""
Reanalysis of Aiemjoy 2022 Figure S11 Panel A:
Individual longitudinal Vi IgG trajectories after natural typhoid infection.

Reads extracted trajectory data and produces exploratory analyses
addressing the question: is there a detectable Vi IgG waning signal?
"""

import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

EXTRACTION_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/fig_s11_extraction")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/fig_s11_reanalysis")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Load extracted data
# =============================================================================

def load_trajectories():
    """Load Panel A trajectories from extracted CSV."""
    with open(EXTRACTION_DIR / "panel_A_trajectories.csv") as f:
        rows = list(csv.DictReader(f))

    # Group by trajectory_id
    trajs = {}
    for r in rows:
        tid = int(r["trajectory_id"])
        trajs.setdefault(tid, {
            "trajectory_type": r["trajectory_type"],
            "color_hex": r["color_hex"],
            "points": [],
        })
        trajs[tid]["points"].append((
            float(r["days_since_fever_onset"]),
            float(r["elisa_units"]),
        ))

    # Sort points within each trajectory by time
    for t in trajs.values():
        t["points"].sort()

    return trajs


# =============================================================================
# Analysis functions
# =============================================================================

def compute_slopes(trajs):
    """Compute log10(EU) slope per day for each individual trajectory."""
    slopes = []
    for tid, t in trajs.items():
        if t["trajectory_type"] != "individual":
            continue
        pts = t["points"]
        if len(pts) < 2:
            continue
        t0, eu0 = pts[0]
        t1, eu1 = pts[-1]
        if t1 <= t0 or eu0 <= 0 or eu1 <= 0:
            continue
        dt = t1 - t0
        log_slope = (np.log10(eu1) - np.log10(eu0)) / dt
        slopes.append({
            "trajectory_id": tid,
            "t_start": t0,
            "t_end": t1,
            "eu_start": eu0,
            "eu_end": eu1,
            "duration_days": dt,
            "log10_slope_per_day": log_slope,
            "fold_change": eu1 / eu0,
            "n_points": len(pts),
        })
    return slopes


# =============================================================================
# Plots
# =============================================================================

def plot_panel_a_reproduction(trajs):
    """Reproduce Panel A as closely as possible to the original."""
    fig, ax = plt.subplots(figsize=(8, 5))

    for tid, t in trajs.items():
        pts = t["points"]
        days = [p[0] for p in pts]
        eus = [p[1] for p in pts]
        if t["trajectory_type"] == "overall":
            ax.plot(days, eus, color=t["color_hex"], linewidth=2.5,
                    label="Overall (model fit)", zorder=10)
        elif t["trajectory_type"] == "individual":
            ax.plot(days, eus, color=t["color_hex"], linewidth=0.8,
                    alpha=0.5, zorder=1)

    ax.set_xlabel("Days since fever onset")
    ax.set_ylabel("Anti-Vi IgG (ELISA units)")
    ax.set_yscale("log")
    ax.set_xlim(0, 580)
    ax.set_ylim(50, 3000)
    ax.set_title("Aiemjoy 2022 Fig S11A: Individual Vi IgG trajectories\n(extracted from PDF vector data)")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "panel_a_reproduction.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved panel_a_reproduction.png")


def plot_slope_analysis(slopes):
    """Analyze the distribution of trajectory slopes."""
    log_slopes = [s["log10_slope_per_day"] for s in slopes]
    durations = [s["duration_days"] for s in slopes]
    fold_changes = [s["fold_change"] for s in slopes]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Histogram of log-slopes
    ax = axes[0, 0]
    ax.hist(log_slopes, bins=30, edgecolor="black", alpha=0.7, color="#9133be")
    ax.axvline(0, color="red", linestyle="--", linewidth=1.5, label="Zero slope (no waning)")
    ax.axvline(np.median(log_slopes), color="blue", linestyle="-", linewidth=1.5,
               label=f"Median = {np.median(log_slopes):.6f}")
    n_decline = sum(1 for s in log_slopes if s < 0)
    n_rise = sum(1 for s in log_slopes if s > 0)
    ax.set_xlabel("Log10(EU) slope per day")
    ax.set_ylabel("Count")
    ax.set_title(f"Slope distribution: {n_decline} declining, {n_rise} rising")
    ax.legend(fontsize=8)

    # 2. Fold-change distribution
    ax = axes[0, 1]
    ax.hist(fold_changes, bins=30, edgecolor="black", alpha=0.7, color="#9133be")
    ax.axvline(1.0, color="red", linestyle="--", linewidth=1.5, label="No change")
    ax.axvline(np.median(fold_changes), color="blue", linestyle="-", linewidth=1.5,
               label=f"Median = {np.median(fold_changes):.2f}")
    ax.set_xlabel("Fold change (end / start)")
    ax.set_ylabel("Count")
    ax.set_title("Fold change distribution")
    ax.legend(fontsize=8)

    # 3. Slope vs duration
    ax = axes[1, 0]
    colors = ["#d62728" if s < 0 else "#2ca02c" for s in log_slopes]
    ax.scatter(durations, log_slopes, c=colors, s=20, alpha=0.6)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.5)
    ax.set_xlabel("Observation duration (days)")
    ax.set_ylabel("Log10(EU) slope per day")
    ax.set_title("Slope vs. observation duration\n(red=declining, green=rising)")
    ax.grid(True, alpha=0.3)

    # 4. Start EU vs slope
    ax = axes[1, 1]
    start_eus = [s["eu_start"] for s in slopes]
    ax.scatter(start_eus, log_slopes, c=colors, s=20, alpha=0.6)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("Starting ELISA units")
    ax.set_ylabel("Log10(EU) slope per day")
    ax.set_title("Slope vs. starting antibody level")
    ax.grid(True, alpha=0.3)

    fig.suptitle("Aiemjoy 2022 Vi IgG waning signal analysis\n(extracted from Fig S11A)",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "slope_analysis.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved slope_analysis.png")


def plot_spaghetti_by_duration(trajs):
    """Spaghetti plot colored by observation duration to highlight
    long-duration trajectories that are most informative about waning."""
    fig, ax = plt.subplots(figsize=(8, 5))

    # Collect individual trajectories with duration info
    indiv = []
    for tid, t in trajs.items():
        if t["trajectory_type"] != "individual":
            continue
        pts = t["points"]
        if len(pts) < 2:
            continue
        duration = pts[-1][0] - pts[0][0]
        indiv.append((duration, pts))

    # Sort so long-duration trajectories plot on top
    indiv.sort()

    # Colormap by duration
    max_dur = max(d for d, _ in indiv)
    cmap = plt.cm.viridis

    for dur, pts in indiv:
        days = [p[0] for p in pts]
        eus = [p[1] for p in pts]
        color = cmap(dur / max_dur)
        ax.plot(days, eus, color=color, linewidth=1.0 + dur / max_dur,
                alpha=0.6)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, max_dur))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, label="Observation duration (days)")

    ax.set_xlabel("Days since fever onset")
    ax.set_ylabel("Anti-Vi IgG (ELISA units)")
    ax.set_yscale("log")
    ax.set_xlim(0, 580)
    ax.set_ylim(50, 3000)
    ax.set_title("Vi IgG trajectories colored by observation duration\n(longer = more informative for waning)")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "trajectories_by_duration.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved trajectories_by_duration.png")


def plot_facet_by_npoints(trajs):
    """Faceted spaghetti plot: columns = time point count, rows = rising/declining."""
    # Group individual trajectories by point count and direction
    groups = {}  # (npts, direction) -> [trajectories]
    for tid, t in trajs.items():
        if t["trajectory_type"] != "individual":
            continue
        n = len(t["points"])
        pts = t["points"]
        if pts[-1][1] < pts[0][1]:
            direction = "declining"
        else:
            direction = "rising"  # includes flat
        groups.setdefault((n, direction), []).append(t)

    npts_values = sorted(set(n for n, _ in groups.keys()))
    directions = ["declining", "rising"]
    dir_colors = {"declining": "#d62728", "rising": "#2ca02c"}
    dir_labels = {"declining": "Net declining", "rising": "Net rising / flat"}

    n_cols = len(npts_values)
    fig, axes = plt.subplots(2, n_cols, figsize=(5 * n_cols, 8),
                             sharex=True, sharey=True)

    for col, npts in enumerate(npts_values):
        for row, direction in enumerate(directions):
            ax = axes[row, col]
            group = groups.get((npts, direction), [])
            color = dir_colors[direction]

            for t in group:
                pts = t["points"]
                days = [p[0] for p in pts]
                eus = [p[1] for p in pts]
                ax.plot(days, eus, color=color, linewidth=1.2, alpha=0.5,
                        marker="o", markersize=3)

            ax.set_yscale("log")
            ax.set_xlim(-10, 580)
            ax.set_ylim(50, 3000)
            ax.grid(True, alpha=0.3)

            # Column titles on top row
            if row == 0:
                total_col = len(groups.get((npts, "declining"), [])) + \
                            len(groups.get((npts, "rising"), []))
                ax.set_title(f"{npts} time points (n={total_col})", fontsize=11)

            # Row labels on left column
            if col == 0:
                ax.set_ylabel(f"{dir_labels[direction]}\n(n={len(group)})\n\n"
                              f"Anti-Vi IgG (EU)", fontsize=9)
            else:
                # Still annotate count
                ax.text(0.02, 0.95, f"n={len(group)}", transform=ax.transAxes,
                        fontsize=9, va="top", fontweight="bold", color=color)

            # X labels on bottom row
            if row == 1:
                ax.set_xlabel("Days since fever onset")

    fig.suptitle("Vi IgG trajectories: time points × direction\n"
                 "(all post-infection; no pre-infection baseline observed)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "trajectories_faceted_by_npoints.png", dpi=150,
                bbox_inches="tight")
    plt.close()
    print(f"Saved trajectories_faceted_by_npoints.png")


# =============================================================================
# Summary statistics
# =============================================================================

def print_summary(trajs, slopes):
    """Print summary statistics."""
    indiv = {k: v for k, v in trajs.items() if v["trajectory_type"] == "individual"}
    overall = {k: v for k, v in trajs.items() if v["trajectory_type"] == "overall"}

    print("=" * 60)
    print("PANEL A REANALYSIS SUMMARY")
    print("=" * 60)
    print(f"Total trajectories: {len(indiv)} individual + {len(overall)} overall")

    n_pts = [len(t["points"]) for t in indiv.values()]
    print(f"Points per trajectory: {min(n_pts)}-{max(n_pts)} "
          f"(median {np.median(n_pts):.0f})")

    durations = [t["points"][-1][0] - t["points"][0][0] for t in indiv.values()
                 if len(t["points"]) >= 2]
    print(f"Observation durations: median {np.median(durations):.0f} days, "
          f"range {min(durations):.0f}-{max(durations):.0f}")

    all_eus = [eu for t in indiv.values() for _, eu in t["points"]]
    print(f"ELISA unit range: {min(all_eus):.0f}-{max(all_eus):.0f}")

    print(f"\n--- Slope analysis (n={len(slopes)}) ---")
    log_slopes = [s["log10_slope_per_day"] for s in slopes]
    n_decline = sum(1 for s in log_slopes if s < 0)
    n_rise = sum(1 for s in log_slopes if s > 0)
    print(f"Declining: {n_decline}, Rising: {n_rise}")
    print(f"Median log10 slope/day: {np.median(log_slopes):.6f}")
    print(f"Mean log10 slope/day:   {np.mean(log_slopes):.6f}")
    print(f"IQR: [{np.percentile(log_slopes, 25):.6f}, "
          f"{np.percentile(log_slopes, 75):.6f}]")

    fold_changes = [s["fold_change"] for s in slopes]
    print(f"Median fold change: {np.median(fold_changes):.3f}")
    print(f"IQR fold change: [{np.percentile(fold_changes, 25):.3f}, "
          f"{np.percentile(fold_changes, 75):.3f}]")

    # Long-duration subset (>200 days)
    long_slopes = [s for s in slopes if s["duration_days"] > 200]
    if long_slopes:
        ls = [s["log10_slope_per_day"] for s in long_slopes]
        n_dec = sum(1 for s in ls if s < 0)
        n_ris = sum(1 for s in ls if s > 0)
        print(f"\n--- Long-duration trajectories (>200 days, n={len(long_slopes)}) ---")
        print(f"Declining: {n_dec}, Rising: {n_ris}")
        print(f"Median log10 slope/day: {np.median(ls):.6f}")
        print(f"Median fold change: {np.median([s['fold_change'] for s in long_slopes]):.3f}")


# =============================================================================
# Main
# =============================================================================

def main():
    trajs = load_trajectories()
    slopes = compute_slopes(trajs)

    print_summary(trajs, slopes)

    print("\nGenerating plots...")
    plot_panel_a_reproduction(trajs)
    plot_slope_analysis(slopes)
    plot_spaghetti_by_duration(trajs)
    plot_facet_by_npoints(trajs)

    print(f"\nAll outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
