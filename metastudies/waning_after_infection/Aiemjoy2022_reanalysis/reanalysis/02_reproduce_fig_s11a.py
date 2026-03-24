"""
02: Reproduce Aiemjoy 2022 Figure S11A from CSV data.

Verify that the CSV data matches the published figure.
Also reconstruct the overall model curve from Table S1 parameters.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/reanalysis")


def teunis_decay(t, y1, t1, alpha, shape):
    """Teunis 2016 power-function decay.

    y(t > t1) = [y1^(1-shape) - (1-shape) * alpha * (t - t1)]^(1/(1-shape))

    From the JAGS model (v9na.model.jags):
      y1^(1-shape) - (1-shape) * alpha * (t - t1)  must stay positive.
    """
    tau = t - t1
    inner = y1 ** (1 - shape) - (1 - shape) * alpha * tau
    inner = np.maximum(inner, 1e-10)  # avoid negative/zero
    return inner ** (1 / (1 - shape))


def teunis_rise(t, y0, y1, t1):
    """Exponential rise phase: y(t <= t1) = y0 * exp(mu * t), mu = log(y1/y0)/t1"""
    mu = np.log(y1 / y0) / t1
    return y0 * np.exp(mu * t)


def full_trajectory(t_arr, y0, y1, t1, alpha, shape):
    """Full two-phase trajectory."""
    y = np.empty_like(t_arr, dtype=float)
    rise = t_arr <= t1
    decay = t_arr > t1
    y[rise] = teunis_rise(t_arr[rise], y0, y1, t1)
    y[decay] = teunis_decay(t_arr[decay], y1, t1, alpha, shape)
    return y


# Table S1 overall Vi IgG parameters (median)
VI_PARAMS_OVERALL = {
    "y0": 126.72,
    "y1": 542.7,
    "t1": 2.92,
    "alpha": 0.0,  # decay rate = 0
    "shape": 1.25,
}


def main():
    long = pd.read_csv(DATA_DIR / "vi_igg_long.csv")
    print(f"Loaded {len(long)} Vi IgG observations, "
          f"{long['index_id'].nunique()} subjects")

    # Filter to longitudinal
    longitudinal = long[long["has_longitudinal"]].copy()
    n_long = longitudinal["index_id"].nunique()

    # Separate by serovar
    typhi = longitudinal[longitudinal["serovar"] == "typhi"]
    para = longitudinal[longitudinal["serovar"] == "paratyphi"]
    n_typhi = typhi["index_id"].nunique()
    n_para = para["index_id"].nunique()

    print(f"Longitudinal: {n_long} subjects ({n_typhi} typhi, {n_para} paratyphi)")

    # Figure: 3 versions for comparison
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    # Panel 1: All longitudinal (matches Fig S11A most closely)
    ax = axes[0]
    for sid, grp in typhi.groupby("index_id"):
        grp = grp.sort_values("days_since_fever_onset")
        ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                color="#9133be", linewidth=0.8, alpha=0.4, marker=".", markersize=2)
    # Overall model curve
    t_model = np.linspace(0.1, 600, 500)
    y_model = full_trajectory(t_model, **VI_PARAMS_OVERALL)
    ax.plot(t_model, y_model, color="#4B0082", linewidth=2.5, label="Overall (Table S1)")
    ax.set_xlabel("Days since fever onset")
    ax.set_ylabel("Anti-Vi IgG (ELISA units)")
    ax.set_yscale("log")
    ax.set_xlim(0, 580)
    ax.set_ylim(10, 5000)
    ax.set_title(f"S. Typhi only (n={n_typhi})\n(should match Fig S11A)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 2: S. Paratyphi
    ax = axes[1]
    for sid, grp in para.groupby("index_id"):
        grp = grp.sort_values("days_since_fever_onset")
        ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                color="#2ca02c", linewidth=0.8, alpha=0.5, marker=".", markersize=2)
    ax.plot(t_model, y_model, color="#4B0082", linewidth=2.5, alpha=0.3,
            linestyle="--", label="S. Typhi overall")
    ax.set_xlabel("Days since fever onset")
    ax.set_yscale("log")
    ax.set_xlim(0, 580)
    ax.set_title(f"S. Paratyphi only (n={n_para})")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 3: All serovars, extended time range
    ax = axes[2]
    for sid, grp in longitudinal.groupby("index_id"):
        grp = grp.sort_values("days_since_fever_onset")
        color = "#9133be" if grp["serovar"].iloc[0] == "typhi" else "#2ca02c"
        ax.plot(grp["days_since_fever_onset"], grp["vi_igg_eu"],
                color=color, linewidth=0.8, alpha=0.4, marker=".", markersize=2)
    t_ext = np.linspace(0.1, 1200, 500)
    y_ext = full_trajectory(t_ext, **VI_PARAMS_OVERALL)
    ax.plot(t_ext, y_ext, color="#4B0082", linewidth=2.5, label="Overall model")
    ax.set_xlabel("Days since fever onset")
    ax.set_yscale("log")
    ax.set_xlim(0, 1200)
    ax.set_title(f"All serovars, full time range (n={n_long})")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Reproducing Aiemjoy 2022 Fig S11A from CSV data",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "02_fig_s11a_reproduction.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved 02_fig_s11a_reproduction.png")

    # Trajectory count comparison
    print(f"\n--- Trajectory count investigation ---")
    print(f"CSV S. Typhi longitudinal: {n_typhi}")
    print(f"PDF extraction: 141 trajectories")
    print(f"Discrepancy: {n_typhi - 141} subjects")

    # Check if time-range filtering explains it
    typhi_in_range = typhi[typhi["days_since_fever_onset"] <= 550]
    n_typhi_range = typhi_in_range["index_id"].nunique()
    print(f"S. Typhi with all obs <= 550 days: {n_typhi_range}")

    # Check if there are subjects with obs starting very late
    first_obs = typhi.groupby("index_id")["days_since_fever_onset"].min()
    print(f"S. Typhi first obs > 200 days: {(first_obs > 200).sum()}")
    print(f"S. Typhi first obs > 400 days: {(first_obs > 400).sum()}")


if __name__ == "__main__":
    main()
