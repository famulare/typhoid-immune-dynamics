"""
01: Prepare Vi IgG data from Aiemjoy 2022 wide-format CSV into clean long format.
"""

import pandas as pd
from pathlib import Path

DATA_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/data")


def main():
    df = pd.read_csv(DATA_DIR / "TypoidCaseData_github_09.30.21.csv")
    print(f"Loaded {len(df)} subjects")
    print(f"Countries: {df['Country'].value_counts().to_dict()}")
    print(f"Serovars: {df['bldculres'].value_counts().to_dict()}")

    # Melt Vi IgG visits to long format
    rows = []
    for _, subj in df.iterrows():
        for visit in range(1, 8):
            vi_col = f"Vi_IgG_visit{visit}"
            time_col = f"TimeInDays_visit{visit}"
            vi_val = subj[vi_col]
            time_val = subj[time_col]
            if pd.notna(vi_val) and pd.notna(time_val):
                rows.append({
                    "index_id": subj["index_id"],
                    "age": subj["age"],
                    "country": subj["Country"],
                    "serovar": subj["bldculres"],
                    "hospitalized": subj["Hospitalized"],
                    "reinf_obs": subj["reinf_obs"],
                    "visit_number": visit,
                    "days_since_fever_onset": time_val,
                    "vi_igg_eu": vi_val,
                })

    long = pd.DataFrame(rows)
    long = long.sort_values(["index_id", "days_since_fever_onset"]).reset_index(drop=True)

    # Add per-subject summary columns
    obs_counts = long.groupby("index_id").size().rename("n_vi_obs")
    long = long.merge(obs_counts, on="index_id")
    long["has_longitudinal"] = long["n_vi_obs"] >= 2

    # Summary
    n_subj = long["index_id"].nunique()
    n_obs = len(long)
    n_long = long.loc[long["has_longitudinal"], "index_id"].nunique()
    n_single = n_subj - n_long

    print(f"\n{'='*60}")
    print(f"Vi IgG LONG-FORMAT SUMMARY")
    print(f"{'='*60}")
    print(f"Total subjects with Vi IgG: {n_subj}")
    print(f"Total observations: {n_obs}")
    print(f"Longitudinal (>=2 obs): {n_long} subjects")
    print(f"Single observation: {n_single} subjects")
    print(f"Country: {long.groupby('country')['index_id'].nunique().to_dict()}")

    print(f"\nBy serovar:")
    for sv, g in long.groupby("serovar"):
        ns = g["index_id"].nunique()
        nl = g.loc[g["has_longitudinal"], "index_id"].nunique()
        print(f"  {sv}: {ns} subjects ({nl} longitudinal), {len(g)} obs")

    print(f"\nObservations per subject:")
    print(obs_counts.value_counts().sort_index().to_string())

    print(f"\nBy reinf_obs (longitudinal only):")
    long_subj = long.loc[long["has_longitudinal"]].drop_duplicates("index_id")
    print(long_subj["reinf_obs"].value_counts().to_string())

    print(f"\nTime range: {long['days_since_fever_onset'].min():.0f} - "
          f"{long['days_since_fever_onset'].max():.0f} days")
    print(f"Vi IgG range: {long['vi_igg_eu'].min():.2f} - "
          f"{long['vi_igg_eu'].max():.1f} EU")
    n_zero = (long["vi_igg_eu"] == 0).sum()
    print(f"Zero values: {n_zero}")

    print(f"\nAge distribution (all Vi subjects):")
    age_bins = pd.cut(long.drop_duplicates("index_id")["age"],
                      bins=[0, 5, 15, 100], labels=["<5", "5-15", "16+"])
    print(age_bins.value_counts().sort_index().to_string())

    # Save
    out_path = DATA_DIR / "vi_igg_long.csv"
    long.to_csv(out_path, index=False)
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    main()
