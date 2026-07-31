<!-- GENERATED FILE — do not hand-edit. Regenerate:
     Rscript -e 'source("tier_specs.R"); write_tier_ladder_md()' -->
# Tier ladder — generated from tier_specs.R

Generated 2026-07-31 from `tier_specs.R` + `dose_response_data.csv` +
`../analysis_data/darton_individual_endpoints.csv` at git `d6160d7`.

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

## Naming rule

A configuration is named by its two DATA switches and nothing else:
`<row set>-<Darton representation>`. `indiv` means the Darton placebo arm enters as
per-subject n=1 rows (the issue-#15 cascade: `ox_inf_indiv` + `ox_fevginf_indiv`)
rather than as a grouped binomial. (`cascade` is deliberately not used as a
configuration name -- in this repo it already denotes the P_inf x P_fev|inf
factorization.)

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

> SCIENTIFIC DECISION, not a code defect. (a) eta Option A (parametric eta_lo/kappa, what
> the .stan implements) vs Option C (fixed eta_fixed_optC, a CSV column no code reads) is
> undecided; eta_detection() is monotone DECREASING in dose while eta_fixed_optC is
> non-monotone (1.00@1e3, 0.62@1e4, 0.94@1.82e4, 0.92@2e4). Because eta multiplies P_inf
> and shares N50_inf in its exponent, a misfit moves the BIOLOGICAL parameters instead of
> failing visibly. (b) ../joint_inference_plan.md Sec 2.6 EXCLUDES Oxford shedding on
> treatment-truncation grounds while the Tier 2 design restores it with eta -- unresolved
> tension. (c) psi (Sec 2.8, adopted 34aac76) is unimplemented, and psi_stool is
> confounded with eta at the single Darton dose. Unblock deliberately with allow_blocked =
> TRUE.
> 

**`t2-indiv`** — Tier 2 rows (+Oxford shedding), Darton placebo individualized

> Everything blocking t2-grouped, plus: Darton contributes no group-2 row here (the
> grouped D-I-plac is dropped as a double count of the 30 ox_inf_indiv rows for the same
> volunteers), so eta is identified by 5 rows -- W-I-3/4 and the three Jin arms. N_obs 85
> = 86 - 1 for that drop.
> 

`fit_tier.R` refuses a blocked configuration unless `--allow-blocked` is passed.
`Rscript fit_tier.R --list` prints this table and these reasons from the registry.

