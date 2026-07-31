# Tier lock — what runs, what it is called, what is blocked

Scope: the dose-response calibration in this directory. **Binding where this file and
the narrative "Step 1 / 1.5 / 2 / 3" ladder disagree.**

> ## PROPOSED LOCK 2026-07-31 — awaiting [Mike]
> *(This blockquote is the lock. The lock ends where the blockquote ends. Stated
> explicitly because `CALIBRATION_WORKFLOW.md`'s Step-2 lock has no end marker, which
> makes it undecidable whether the lines after it are inside it.)*
>
> **1. A fit configuration is named by its two DATA switches, and nothing else.**
> `<row set>-<Darton representation>`:
>
> | key | `tier_col` | `individualize_darton` |
> |---|---|---|
> | `t1-grouped` | `tier1_active` | `FALSE` |
> | `t1-indiv`   | `tier1_active` | `TRUE`  |
> | `t2-grouped` | `tier2_active` | `FALSE` |
> | `t2-indiv`   | `tier2_active` | `TRUE`  |
>
> `indiv` = the Darton placebo arm enters as per-subject n=1 rows (the issue-#15
> cascade: groups `ox_inf_indiv` + `ox_fevginf_indiv`) instead of a grouped binomial.
> **`cascade` is deliberately NOT a spec name:** in this repo it already denotes the
> P_inf × P_fev|inf factorization, and reusing it would mint a new ambiguity while
> fixing an old one. `t1`/`t2` are the CSV column names, so the label is derived from
> a code identifier rather than from prose.
>
> **2. Run directories are `<data config>__<model stage>`,** e.g.
> `results/t1-indiv__phi/`, with `-prior` appended for the prior predictive. The stage
> token names the parameter increments that are ACTIVE, is **derived** from the
> `.stan`'s `parameters{}`, and is **asserted** against the registry:
> `validate_tier_spec()` hard-errors when the two disagree. A directory therefore
> cannot misdescribe its own model. Current tokens:
>
> | token | parameter set |
> |---|---|
> | `phi` | `phi0_a`, `phi0_b`, `beta_phi` pinned — **the current `.stan`** |
> | `phi-rho` | + `grand_overdispersion_rho` (Step 2: locked, unimplemented) |
> | `phi-eta` | + `eta_lo`, `kappa` reached, which needs a group-2 row |
> | `…-psi` | + `psi_stool`, `psi_late` (Sec 2.8: adopted `34aac76`, unimplemented) |
>
> An increment counts only when the `.stan` declares its parameters AND the tier's
> data actually reach them: `eta_lo`/`kappa` have been declared since before Tier 1
> but appear in exactly one likelihood branch (group 2), so **declared-but-unreached
> is not a model stage.** That is why `t1-*` read `phi` and `t2-*` read `phi-eta`.
>
> **3. A stage is never reused.** A new stage is a new directory, so the pre-increment
> fit is not clobbered. `results/` accretes one directory per (data config × stage).
> This replaces the only naming rule the repo had — "don't overwrite `results/tier1`"
> (`tier1.5_plan.md`) — which the old driver violated on every run. Within a stage,
> runs are distinguished by the manifest's git SHA and input hashes, not by renaming.
>
> **4. Observation counts are GENERATED, never typed.** `tier_specs.R` is the
> registry; `TIER_LADDER.md` is generated from it by calling `build_stan_data()` per
> spec and counting rows, and the generator hard-errors when a computed count
> disagrees with the spec's declared `expect` (both `N_obs` and the per-group table).
> **No other document in this repo states a tier observation count** — they cite
> `TIER_LADDER.md`. It is tracked in git precisely because `results/` is not, so a
> change to the data or to a default appears as a diff.
>
> **5. Every run directory carries `run_manifest.json`** — tier spec, ordered
> `obs_id` (which IS the `p_pred` column order), input md5s, resolved `pr_*` prior
> scalars, seeds and sampler settings, git SHA + dirty flag, R/CmdStan/package
> versions. It is written even when a caller supplies none, so a run cannot be
> anonymous. `audit_run_dirs()` reports which directories are regenerable and from
> what code state; **absence of a manifest is the unambiguous legacy marker.**
>
> **6. The fit that existed before this lock is `t1-indiv`.** [observed] 80
> observations = 24 grouped rows + 30 `ox_inf_indiv` + 26 `ox_fevginf_indiv`, over 6
> likelihood groups. It was labelled "Tier 1", reported as 25 observations, and
> written to `results/tier1/` because `individualize_darton = TRUE` is the
> `build_stan_data()` default and no driver overrode it.
>
> **7. Historical run directories keep their on-disk names.** They are referenced from
> dated lab-notebook entries; renaming them would falsify an append-only log.
>
> | on disk | what it is |
> |---|---|
> | `results/tier1` | the `t1-indiv` fit as of 2026-07-31; superseded by `results/t1-indiv__phi` |
> | `results/tier1_prior` | its prior-predictive companion |
> | `results/tier1_minimal_phi` | C0+C1+C2, with the retired scalar `phi_md` (2026-06-23) |
> | `results/tier1_pre_eumL` | pre-EU/mL dimensionless-CoP fit (2026-06-23) |
> | `results/recovery/{point,tier1}` | known-truth recovery runs, not data fits |
> | `results/scenarios/*` | row-filter / prior-override sensitivities on `t1-indiv` |
>
> None of these carry a manifest, so none can be tied to a code state. Refit or
> delete; do not cite their numbers as current.
>
> **8. Runnable vs blocked.**
> - `t1-indiv` — **RUNNABLE.** This is the fit.
> - `t1-grouped` — **RUNNABLE** as of this lock. It is the diagnostic contrast showing
>   what individualizing Darton buys. It could not previously be figured at all
>   (item 9.2).
> - `t1-indiv` + `grand_overdispersion_rho` — designed and LOCKED
>   (`CALIBRATION_WORKFLOW.md` Step 2, `../joint_inference_plan.md` §5.5), **not
>   implemented**: the `.stan` has no such parameter.
> - `t2-grouped`, `t2-indiv` — **DECLARED AND BLOCKED.** Not "not started". The
>   registry carries the reason and `fit_tier.R` refuses them without
>   `--allow-blocked`.
>
> **9. What blocked the tiers, and what was fixed.** [observed]
> 1. **Darton double count — FIXED.** `build_stan_data()` dropped only the literal
>    `"D-F-plac"` when individualizing, so at `tier2_active` the grouped `D-I-plac`
>    row (n=30, y=19, stool shedding, η-corrected) sat in the likelihood beside the 30
>    `ox_inf_indiv` rows (`bact_or_stool`, 26/30) for the **same 30 volunteers**. Both
>    grouped rows are now replaced (`DARTON_PLACEBO_GROUPED_OBS`), and
>    individualization runs before `drop_obs`/`keep_obs` — which also closed a trap
>    where `drop_obs = "D-F-plac"` silently removed all 56 per-subject rows.
>    `t2-indiv` is therefore 85, not 86. `t1-indiv` is unaffected (`D-I-plac` is
>    `tier1_active == 0`) and is **bit-identical** to the pre-lock fit.
> 2. **Figure suite hard-errored on every grouped Darton row — FIXED.**
>    `validate_curve_specs()` requires each `obs_id` to match exactly one grouping,
>    and the Darton spec's regex was `"^D-(I|FgI)-plac-"` — trailing hyphen, no `F` —
>    so `D-F-plac` and `D-I-plac` matched **zero** groupings. That blocked
>    `t1-grouped`, `t2-grouped` and `t2-indiv`: three of four rungs, not just Tier 2.
>    Added a `$`-anchored grouped column and `drop_empty_columns`.
> 3. **η missing from the figure math — OPEN, mechanical.** Group 2's likelihood is
>    `eta(D) × P_inf` but the curve layer has no η branch, so a `t2-*` grid would draw
>    an uncorrected P(infection|D) against η-corrected data. `mm_obs_prob()` already
>    has η, so the parity gate passes and only the drawn curve is wrong.
> 4. **η Option A vs C — OPEN, SCIENTIFIC.** Parametric `eta_lo`/`kappa` (what the
>    `.stan` implements) vs fixed `eta_fixed_optC` (a CSV column **no code reads**).
>    `eta_detection()` is monotone *decreasing* in dose; `eta_fixed_optC` is
>    non-monotone (1.00@1e3, **0.62@1e4**, 0.94@1.82e4, 0.92@2e4). Because η
>    multiplies `P_inf` and shares `N50_inf` in its exponent, a misfit **moves the
>    biological parameters instead of failing visibly.**
> 5. **Prior question — OPEN.** `../joint_inference_plan.md` §2.6 *excludes* Oxford
>    shedding on treatment-truncation grounds while the Tier 2 design *restores* it
>    with an η correction. Those are in tension and this lock does not resolve it.
>
> **10. Unlock condition for `t2-*`:** 9.3 fixed with a gate that fails if it
> regresses (a `validate_curve_specs()` pass over the `t2-*` obs tables, and a parity
> check that includes group 2 — note `obs_prob_R()` currently *refuses* group 2, so η
> has only a two-implementation check on exactly the branch Tier 2 turns on); and 9.4
> plus 9.5 decided and written into this file.
>
> **11. INERT parameters, stated once.** [observed]
> - **`sigma_study` is inert at EVERY spec**, not just Tier 1: declared in the `.stan`
>   with a prior and used in **zero** likelihood terms. There is no study index
>   anywhere in the Stan data. The cohort random effect it was for is LOCKED-WITHDRAWN
>   (`cohort_random_effects_design.md`), and Step 2's lock already says to delete it.
> - `eta_lo`, `kappa` are inert unless a group-2 (`ox_inf`) row is active — i.e. in
>   every `t1-*` spec.
> - `beta_phi` is pinned to 1 in code.
> - **`phi_md` does not exist.** Retired at C3 (`d1880a8`), replaced by
>   `phi0_a`/`phi0_b`. Any document naming it is describing history.
> - **7 likelihood groups, not 5** (`.GROUP_CODE`, `data_prep.R`).
>
> **12. How to falsify anything above.** Every count is [observed] — computed
> 2026-07-31 from `dose_response_data.csv` and
> `../analysis_data/darton_individual_endpoints.csv` through `build_stan_data()`. Run
> `Rscript -e 'source("tier_specs.R"); validate_all_tiers()'`; it re-derives every
> count and fails loudly on any disagreement. Item 6's count is stated here on purpose
> so the lock is self-contained and falsifiable; it is the one count outside the
> generated table.
>
> *(end of lock)*

## Why this lock exists

`individualize_darton` (`data_prep.R`, default `TRUE`) appeared in **zero** `.md`
files in this repo. [observed] Every document stating 25 or 31 observations was
describing `FALSE`; every document describing the Darton cascade was describing
`TRUE`. The fitted model was the 80-observation configuration while being labelled
and filed as the 25-observation one. One undocumented default was the whole ambiguity.

Counts as of 2026-07-31 [observed]: (`tier1_active`, TRUE) = 80;
(`tier1_active`, FALSE) = 25; (`tier2_active`, TRUE) = 85 after the double-count fix;
(`tier2_active`, FALSE) = 31.

## What this lock does NOT decide

1. `CALIBRATION_WORKFLOW.md` Step 2 calls `grand_overdispersion_rho` a "single shared
   **concentration**", while three lines later — and in `../joint_inference_plan.md`
   §5.5 — ρ is defined as the **ICC** with `k = (1-ρ)/ρ` derived, and the same
   section argues explicitly *against* attaching "overdispersion" to a concentration.
   A self-contradiction inside a lock. [Mike] to reword.
2. That lock has no end marker, so its scope is undecidable. Process suggestion:
   delimit locks with `>` blockquotes, as `cohort_random_effects_design.md` and this
   file do.
3. `../joint_inference_plan.md` disagrees with itself in the section that claims its
   counts are verified: one place says 24 grouped observations, another 25 active,
   another "19 to 25". Which was right in March, and against which row set, is
   [Mike]'s to settle.
4. Whether `sigma_study` is deleted now or at Step 2 (two existing locks say Step 2).
5. Whether the `t2-*` rungs should be fit at all (items 9.4, 9.5).
6. Whether Levine's paired endpoints are factorized into marginal + conditional.
   `Lev-F-k` and `Lev-I-k` are the same volunteers (same `cohort_id`, same n,
   fever ⊆ infection in all four trials) entered as **independent binomials** —
   the one place the repo's own no-double-counting rule is not applied. Hornick and
   Darton already comply. ψ is the precondition, and the row *count* would not change
   (4 rows convert), but the group table and those rows' `n`/`y` would, and the
   fever-among-infected numerators are not in the CSV.
