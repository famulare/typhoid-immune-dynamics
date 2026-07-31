<!-- GENERATED FILE — do not hand-edit. Regenerate:
     Rscript -e 'source("tier_specs.R"); write_tier_ladder_md()' -->
# Tier ladder — generated from tier_specs.R

Generated 2026-07-31 from `tier_specs.R` + `dose_response_data.csv` +
`../analysis_data/darton_individual_endpoints.csv` at git `01df289`.

Counts are **computed** by calling `build_stan_data()` for each spec and counting
the rows it returns; the generator hard-errors when a computed count disagrees
with the spec's declared `expect`. If a number here is wrong, the registry or the
data is wrong — fix it there and regenerate. **No other document in this repo
states a tier observation count.**

| key | tier_col | individualize_darton | N_obs | grouped | individual | stage | status |
|---|---|---|---|---|---|---|---|
| t1-grouped   | tier1_active | FALSE        | 25           | 25           |  0           | phi          | runnable     |
| t1-indiv     | tier1_active | TRUE         | 80           | 24           | 56           | phi          | runnable     |
| t2-grouped   | tier2_active | FALSE        | 31           | 31           |  0           | phi-eta      | blocked      |
| t2-indiv     | tier2_active | TRUE         | 85           | 29           | 56           | phi-eta      | blocked      |

## Likelihood groups per tier

| key | groups |
|---|---|
| t1-grouped                                 | hornick_cond=1 md_fev=11 md_inf=6 ox_fev=7 |
| t1-indiv                                                                       | hornick_cond=1 md_fev=11 md_inf=6 ox_fev=6 ox_fevginf_indiv=26 ox_inf_indiv=30 |
| t2-grouped                                          | hornick_cond=1 md_fev=11 md_inf=6 ox_fev=7 ox_inf=6 |
| t2-indiv                                                                                | hornick_cond=1 md_fev=11 md_inf=6 ox_fev=6 ox_fevginf_indiv=26 ox_inf=5 ox_inf_indiv=30 |

## Run directories

| key | run_dir |
|---|---|
| t1-grouped              | results/t1-grouped__phi |
| t1-indiv              | results/t1-indiv__phi |
| t2-grouped                  | results/t2-grouped__phi-eta |
| t2-indiv                  | results/t2-indiv__phi-eta |

Run dirs are `<data config>__<model stage>`. The stage token is derived from the
`.stan`'s `parameters{}` and asserted against the registry, so a directory name
cannot drift from the model that produced it. Prior-predictive companions take a
`-prior` suffix. See `TIER_LOCK.md` for the naming rule and the blocked rungs.

