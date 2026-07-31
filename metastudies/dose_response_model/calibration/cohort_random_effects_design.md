# Cohort random effects: evaluation of joint_inference_plan.md §5.5, and a revision

Written 2026-07-31. Triggered by the Gilman residual surfaced in
`results/tier1/dose_response_grid_maryland.png`.

> ## RESOLVED — LOCKED 2026-07-31 [Mike]
>
> **The study/cohort random effect is WITHDRAWN, not deferred.** Tier 2 extra-binomial
> variation is handled by **one overdispersion parameter: the beta-binomial ICC `grand_overdispersion_rho`**
> (concentration `grand_concentration_k = (1-rho)/rho` derived; rho = 0 is exactly the binomial).
>
> Applied to: `joint_inference_plan.md` §5.5 and the note at §4.2;
> `CALIBRATION_WORKFLOW.md` Step 2 and the tier ladder. `sigma_study` should be
> deleted from the `.stan` when Step 2 is implemented.
>
> `cohort_id` was still added to `dose_response_data.csv` (2026-07-31). It is not passed
> to Stan and no likelihood term reads it; post-fit it keys the LOO units so model
> comparison does not count the same volunteers twice (80 rows -> 46 units).
>
> The analysis below is retained as the record of why.

**Verdict: §5.5 as written should not be implemented.** Its index is wrong, its
primary formulation has almost no leverage, and its motivating premise does not
survive a test. Its *fallback* recommendation (beta-binomial) is right, but
for a reason the plan does not state. Details and a revision below.

---

## 1. The motivating premise — PARTLY REFUTED, see correction at end of section

§5.5 and line 518 both rest on: *"The trial-to-trial variability (25% to 55%) at
the same dose suggests either batch effects, temporal shifts in cohort immunity,
or random variation. This informs the study-level random effect."*

Tested directly on the four Levine trials (all 10^5 CFU, all T>=38.3):

| obs_id | year | n | y | rate |
|---|---|---|---|---|
| Lev-F-1 | 1970 | 26 | 13 | 0.500 |
| Lev-F-2 | 1971 | 33 | 10 | 0.303 |
| Lev-F-3 | 1972 | 22 | 12 | 0.545 |
| Lev-F-4 | 1973 | 16 |  4 | 0.250 |

- chi-square homogeneity: **X2 = 5.80, df = 3, p = 0.122**
- Fisher exact 2x4: **p = 0.125**
- adding `Gil-F-rest` (5 cohorts, same dose/threshold): **p = 0.207**

Logit-scale variance decomposition over those 5 cohorts: total 0.305, mean
binomial sampling variance 0.236, **excess 0.069 -> sigma_cohort ~ 0.26**.

So on FEVER the 25%-55% range is what n = 16-33 binomial sampling looks like.

> **CORRECTION 2026-07-31.** The above tests fever only. The INFECTION endpoint, same
> four Levine trials, same men, **rejects homogeneity: X2 = 10.71, df = 3, p = 0.013**
> (rates 0.731 / 0.455 / 0.773 / 0.375 against a flat model prediction of 0.676).
> Cohort heterogeneity IS detectable; this section's original conclusion was wrong
> because it tested the one endpoint where it was not.
>
> The fever/infection residuals are perfectly rank-concordant across the four trials,
> but that is largely MECHANICAL (fever is nested inside infection in the cascade, and
> the counts respect it: 13<=19, 10<=15, 12<=17, 4<=6), so it is not independent
> evidence. The p = 0.013 stands on its own.
>
> **The lock is unaffected**, because it never rested on this argument -- Sections 4
> and 5 (identifiability: 16 cohorts / 24 observations / 10 singletons; Hornick's
> cohorts ARE the dose ladder) carry it. But nobody should cite "p = 0.122, no
> heterogeneity" as settled.

## 2. The index is wrong (this is the substantive error)

§5.5 specifies `epsilon_s ~ Normal(0, sigma_study)` at the **study** level, while
motivating it with variation **within** Levine. A study-level effect gives Levine
one offset and cannot absorb its own between-trial spread. The plan's motivation
and its index are inconsistent.

`study x year` does not fix it either:

| study | levels by `study` | levels by `study x year` | actual challenge groups |
|---|---|---|---|
| Levine | 1 | 4 | 4 (one per trial; fever + infection rows are the SAME men) |
| Hornick | 1 | **1** (all rows are year 1965) | **5** (one per dose: 1e3, 1e5, 1e7, 1e8, 1e9) |

Hornick's rows are one paper, one year, but five different groups of volunteers
distinguished only by dose. **The cohort is not derivable from any existing column.**
It needs an explicit `cohort_id`.

Two further structural facts the plan does not encode:
- `Lev-F-k` and `Lev-I-k` are the same volunteers (e.g. 13/26 fever and 19/26
  infection in trial 1). Currently independent binomials.
- `H-FgI-7` (28) is nested inside `H-I-7` (30). The cascade factorization already
  handles this one correctly.
- `H-F-5` (n = 116) is itself a pool of many challenges across years, so it is not
  a single cohort even in principle. Any `cohort_id` assigned to it is a fiction.

## 3. The primary formulation has no leverage

§5.5's first equation perturbs immunity: `CoP * exp(epsilon_s)`. At the fitted
`gamma_inf ~ 0.18` the CoP channel is heavily compressed — measured on
`results/tier1`, at the Maryland 10^5 dose, moving CoP from 1 to 53 buys only a
2.40x change in P(fever). A CoP-scale offset therefore needs implausibly large
`epsilon` to move a rate at all.

**Reviewer 2's logit-scale revision is correct and should become the primary form.**
The CoP-scale version should be deleted, not retained as the headline equation.

## 4. Identifiability: the part the plan never assesses

Tier 1 has **24 group-level observations** (plus 56 Darton individual n=1 rows).
Enumerating actual challenge groups (verified programmatically, 2026-07-31):

| cohort | rows |
|---|---|
| Gilman (3 H-strata + infection) | 4 |
| Hornick 1e7 (infection + fever\|inf, nested) | 2 |
| Levine 1 / 2 / 3 / 4 (fever + infection each) | 2 each |
| Hornick 1e3 / 1e5 / 1e8 / 1e9 | 1 each |
| Waddington 1e3 / 1e4 | 1 each |
| Jin ctrl / ViPS / ViTT | 1 each |
| Gibani | 1 |
| **total** | **16 cohorts, 24 rows** |

**16 cohorts for 24 group-level observations, 10 of them singletons.** A free
per-cohort offset is close to one parameter per datum, and for those 10 singletons
the offset absorbs the residual exactly (up to shrinkage) — it is not estimating
heterogeneity there, it is interpolating.

Only 6 cohorts have >1 row, and 5 of those 6 get their second row from the
fever/infection pairing rather than from replication at a shared dose.

The specific danger: **Hornick's five cohorts ARE the dose ladder** (1e3, 1e5, 1e7,
1e8, 1e9). A free offset per Hornick cohort competes directly with `N50_inf` and
`alpha_inf` for the same variance — the random effect can flatten the dose-response
and improve the fit while destroying the quantity the study exists to estimate.

Meanwhile sigma is identified in exactly **one** cell: the Maryland 10^5 / T=38.3
replicate group (Levine x4 + Gilman). So sigma would be learned in one place and
applied where it eats signal.

## 5. What it will and will not fix

It will **not** fix the Gilman residual that prompted this. `Gil-F-Hlo` (22/36 =
0.611 observed vs 0.429 fitted) and `Gil-F-Hhi` (4/17 = 0.235 vs 0.258 fitted,
which fits well: P(Bin(17, 0.258) <= 4) = 0.54) are **strata of one cohort**. A
cohort offset moves all three Gilman rows together — it raises the stratum that is
already fitting while chasing the one that is not. A within-cohort contrast cannot
be fixed by a between-cohort term.

Scale check: the `Gil-F-Hlo` gap is 0.451 - (-0.286) = **0.737 logits**, against
sigma ~ 0.26. That is a 2.8-sigma cohort excursion. It will not be absorbed.

Separately: the Gilman stratum contrast is **not** rejected by the data anyway.
Testing H0 that the true susceptible:immune ratio equals the model's 1.67x gives
**p = 0.33**; the observed 2.60x has 95% CI 1.06-6.36 and the analysis was post-hoc.
There is no stratification anomaly to fix.

---

## 6. Proposed revision

### 6a. Add `cohort_id` to `dose_response_data.csv` regardless

Independent of any random effect, the dataset should record which rows share
volunteers. It is provenance, it is currently unrecoverable from the CSV, and it
is a precondition for any future hierarchical term. Assign explicitly, and mark
`H-F-5` as a known pooled row that is not a true cohort.

### 6b. Do overdispersion with ONE parameter first, not per-cohort offsets

Given 16 cohorts / 24 observations with 10 singletons, a per-cohort random effect
has too much freedom. A single shared overdispersion parameter inflates variance
without giving any individual dose point its own free offset, so it cannot flatten
the dose-response. This is the plan's fallback recommendation and it is the right
starting point — but the plan justifies it as "simpler", when the actual reason is
**identifiability**, and that reason should be recorded.

Restrict it to the Maryland group-level rows, where the replication that motivates
it lives. Do not apply it to the Darton n=1 individual rows, where overdispersion
on a Bernoulli is not identified.

### 6c. The per-cohort random effect is withdrawn [LOCKED 2026-07-31]

Not "deferred pending more data" — withdrawn. Two of the routes that might have
revived it are closed:

- **Woodward's veteran stratification cannot be added.** It is a REVIEW of the same
  16-year program: its ~305 controls at 10^5 overlap the 277 already in the fit
  (Hornick 116, Levine 97, Gilman 64), so adding it double-counts. It also publishes
  percentages only, and never defines its clinical endpoint, so T (and therefore
  phi(T,D)) is undetermined. Its value is concordance on the mixture weight
  (200/305 = 0.66 vs Gilman's 36/53 = 0.68) — a DIFFERENT stratification reaching a
  similar number, documented in priors.yaml, not an interchangeable anchor.
- **Cohort-varying pi cannot span the observed spread.** At the fitted parameters the
  mixture's entire reachable range at Maryland 10^5 is 0.260-0.430 (pi = 0 to pi = 1).
  Observed points run 0.235-0.611: four above the ceiling, two below the floor, one
  inside. The band's width is set by CoP_imm^gamma = 11.3^0.18 = 1.52, i.e. by gamma,
  not by pi. (Caveat: only `Gil-F-Hlo` departs individually significantly, p = 0.022.)
  And gamma is not free to change — the model deliberately shares it across eras on
  the assumption that human immunology is the same in both [Mike].

For the record, had it been adopted the form would have been:

    logit(p_i) = logit(p_model_i) + epsilon_{c(i)}
    epsilon_c  = sigma_cohort * z_c,   z_c ~ Normal(0, 1)     # NON-CENTERED

- **logit scale**, not CoP scale (§3).
- **indexed by `cohort_id`**, shared across a cohort's fever / infection /
  conditional rows (§2).
- **non-centered.** This model has prior form here: the ~99% divergence pathology
  diagnosed in `tier1_pathology_diagnosis.md` came from a density cliff. A centered
  hierarchical term with 16 weakly-identified groups will reintroduce divergences.
- **prior `HalfNormal(0, 0.5)`** — the implied sigma ~ 0.26 sits at 0.52 sd, well
  supported without being at the boundary. Do not loosen it. `HalfNormal(0, 0.3)`
  is defensible given how weak the signal is.

### 6d. Gates (this must not be adopted on improved fit alone)

1. **`log10_N50_inf` and `alpha_inf` must not move materially.** This is the
   primary gate. Hornick's cohorts are the dose ladder; if the dose-response
   parameters shift when the effect is added, it is absorbing signal, not noise.
2. **priorsense on `sigma_cohort`.** Given p = 0.12 on the premise, expect it to be
   prior-dominated. If it is, say so: it is then an uncertainty-honesty device, not
   a finding.
3. **LOO** must improve, using the existing `compute_loo_units()` grouping (which
   already merges the Hornick marginal + conditional into one unit).
4. **Gilman residual is NOT a gate.** Pre-register that it is expected to persist
   (§5). If it disappears, something else changed and that needs explaining.

---

## 7. Summary table

| §5.5 element | verdict |
|---|---|
| motivating premise (Levine 25-55%) | **fails** — p = 0.122, ~77% is sampling noise |
| study-level index | **wrong** — motivated by within-study variation |
| `study x year` as a fix | **also wrong** — collapses Hornick's 5 dose cohorts to 1 |
| CoP-scale `exp(epsilon)` form | **no leverage** at gamma ~ 0.18 |
| Reviewer 2 logit-scale revision | **correct** — promote to primary |
| beta-binomial fallback | **ADOPTED** — right call, wrong reason: identifiability, not simplicity |
| identifiability assessment | **absent** — 16 cohorts / 24 obs, 10 singletons |
| will it fix the Gilman residual | **no** — that is a within-cohort contrast |
