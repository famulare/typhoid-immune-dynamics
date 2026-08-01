<!-- GENERATED FILE — do not hand-edit. Regenerate:
     Rscript -e 'source("tier_specs.R"); write_tier_ladder_md()' -->
# Tier ladder — generated from tier_specs.R

Generated 2026-07-31 from `tier_specs.R` + `dose_response_data.csv` +
`../analysis_data/darton_individual_endpoints.csv` at git `e31dcf5`.

Counts are **computed** by calling `build_stan_data()` for each spec and counting
the rows it returns; the generator hard-errors when a computed count disagrees
with the spec's declared `expect`. If a number here is wrong, the registry or the
data is wrong — fix it there and regenerate. **No other document in this repo
states a tier observation count.**

| key | tier_col | individualize_darton | N_obs | grouped | individual | stage | status |
|---|---|---|---|---|---|---|---|
| t1-grouped   | tier1_active | FALSE        |  25          | 25           |   0          | phi-rho      | runnable     |
| t1-indiv     | tier1_active | TRUE         |  80          | 24           |  56          | phi-rho      | runnable     |
| t1-indiv-vax | tier1_active | TRUE         | 177          | 24           | 153          | phi-rho-vax  | runnable     |
| t2-grouped   | tier2_active | FALSE        |  31          | 31           |   0          | phi-rho-eta  | blocked      |
| t2-indiv        | tier2_active    | TRUE            |  85             | 29              |  56             | phi-rho-eta-psi | runnable        |
| t2-indiv-vax        | tier2_active        | TRUE                | 182                 | 29                  | 153                 | phi-rho-eta-psi-vax | runnable            |

## Likelihood groups per tier

| key | groups |
|---|---|
| t1-grouped                                 | hornick_cond=1 md_fev=11 md_inf=6 ox_fev=7 |
| t1-indiv                                                                       | hornick_cond=1 md_fev=11 md_inf=6 ox_fev=6 ox_fevginf_indiv=26 ox_inf_indiv=30 |
| t1-indiv-vax                                                                   | hornick_cond=1 md_fev=11 md_inf=6 ox_fev=6 ox_fevginf_indiv=63 ox_inf_indiv=90 |
| t2-grouped                                          | hornick_cond=1 md_fev=11 md_inf=6 ox_fev=7 ox_inf=6 |
| t2-indiv                                                                                | hornick_cond=1 md_fev=11 md_inf=6 ox_fev=6 ox_fevginf_indiv=26 ox_inf=5 ox_inf_indiv=30 |
| t2-indiv-vax                                                                            | hornick_cond=1 md_fev=11 md_inf=6 ox_fev=6 ox_fevginf_indiv=63 ox_inf=5 ox_inf_indiv=90 |

## Run directories

| key | run_dir |
|---|---|
| t1-grouped                  | results/t1-grouped__phi-rho |
| t1-indiv                  | results/t1-indiv__phi-rho |
| t1-indiv-vax                      | results/t1-indiv-vax__phi-rho-vax |
| t2-grouped                      | results/t2-grouped__phi-rho-eta |
| t2-indiv                          | results/t2-indiv__phi-rho-eta-psi |
| t2-indiv-vax                              | results/t2-indiv-vax__phi-rho-eta-psi-vax |

## Naming rule

A configuration is named by its DATA switches and nothing else:
`<row set>-<Darton representation>[-vax]`. `indiv` means the Darton placebo arm
enters as per-subject n=1 rows (the issue-#15 cascade: `ox_inf_indiv` +
`ox_fevginf_indiv`) rather than as a grouped binomial. (`cascade` is deliberately
not used as a configuration name -- in this repo it already denotes the
P_inf x P_fev|inf factorization.) The optional third switch, `-vax`
(`include_vaccine_arms`), ALSO individualizes Darton's M01ZH09 and Ty21a arms
(+vaccine-terms, tier1.5_plan.md) -- it requires `indiv` (there is no grouped
representation of those arms in any active tier).

Run directories are `<config>__<model stage>`, with `-prior` appended for the prior
predictive. The stage token lists the parameter increments that are ACTIVE; it is
**derived** from the `.stan`'s `parameters{}` and **asserted** against the registry,
so a directory cannot misdescribe the model that produced it, and a new stage is a
new directory rather than an overwrite. An increment counts only when the `.stan`
declares its parameters AND the data reach them -- `eta_lo`/`kappa` sit in exactly
one likelihood branch (group 2), so declared-but-unreached is not a stage. Within a
stage, runs are distinguished by `run_manifest.json` (git SHA, input hashes, seeds),
not by renaming.

## Blocked configurations

**`t2-grouped`** — Tier 2 rows (+Oxford shedding), Darton placebo grouped

> RETIRED 2026-07-31 [Mike, tier2_plan.md]: Tier 2 is individualized-Darton only going
> forward -- no new grouped configuration will be built. Kept in the registry (not
> deleted) so the reason stays visible, same principle as every other blocked/retired
> rung. The grouped-vs-individualized contrast this entry existed for is already made by
> t1-grouped vs t1-indiv; t2-grouped would only duplicate that contrast one increment
> later. See t2-indiv / t2-indiv-vax for the live Tier 2 configurations (eta + psi
> implemented, tier2_plan.md).
> 

`fit_tier.R` refuses a blocked configuration unless `--allow-blocked` is passed.
`Rscript fit_tier.R --list` prints this table and these reasons from the registry.

