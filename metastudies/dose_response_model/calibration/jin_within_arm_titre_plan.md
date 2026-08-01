# Within-arm titre quadrature: replace the arm-GMT plug-in for Oxford grouped rows

**Status:** DRAFT 2026-08-01, awaiting lock [Mike]. Builds on `t2-indiv-vax`, stage
`phi-rho-eta-psi-vax`, N_obs=182. Scope is `obs_prob()` groups 1 (`ox_fev`) and 2
(`ox_inf`) only. No new rows, no new free parameters, no change to any tier's row
membership.

**Decision recorded 2026-08-01 [Mike]:** adopt the quadrature. Do **not** run the
`rho = 0` refit or the grouped-Darton refit; `grand_overdispersion_rho` stays as
specified.

## Rationale

`joint_inference_plan.md` §5.3 has always specified

$$p_j = \int P_k(D, \text{CoP})\, f_j(\text{CoP})\, d\text{CoP}$$

and has always been implemented as the plug-in $P_k(D, \text{GMT}_j)$. This is a
documented approximation in the current code. This increment removes it.

**This is a small change, not a model fix.** At the
`t2-indiv-vax` posterior medians, integrating over the within-arm spread moves the
fitted Jin probabilities by 0.6–2.2% relative (table in `joint_inference_plan.md`
§5.3). `gamma_inf` and `gamma_fevginf` will not move meaningfully. Two reasons to
make the change:

1. The approximation's bias direction depends on the regime. The
   dose-response is convex in $\log\text{CoP}$ where $\alpha L / \text{CoP}^\gamma < 1$
   and concave above it, and Jin's arms straddle that boundary), so the plug-in's
   error cannot be assigned a sign in advance and can change sign if $\gamma$ or
   the doses change.
2. The same quadrature code is needed for a future within-arm likelihood term. The
   information lost by grouping is the *joint* of titre and outcome; see
   `joint_inference_plan.md` §5.3, "What quadrature does not recover." That term is
   flagged but not adopted.

## 1. Data layer

Add one column to `dose_response_data.csv`:

| column | meaning | default |
|---|---|---|
| `CoP_sd_log10` | within-arm sd of $\log_{10}$ CoP for this row | `0` |

`0` means "point mass at `CoP`". The quadrature then collapses to the current
plug-in **exactly**, so every row and every tier that does not set it is
bit-identical to the current fit. This mirrors the `psi_active` gating pattern
(`tier2_plan.md` decision B) and the `n=1` beta-binomial no-op.

Values, backed out of Jin Suppl. Table S2 reported GMT + 95% CI + arm $n$, assuming
a normal-theory CI on the log geometric mean:

$$\sigma_j = \frac{\log_{10}(\text{hi}_j/\text{lo}_j)}{2 \times 1.96}\sqrt{n_j}$$

| obs_id | GMT (EU/mL) | 95% CI | n | `CoP_sd_log10` |
|---|---|---|---|---|
| `J-F-ctrl`, `J-I-ctrl` | 8.0 | 5.2–12.2 | 31 | 0.526 |
| `J-F-ViPS`, `J-I-ViPS` | 140.5 | 91.0–216.9 | 35 | 0.569 |
| `J-F-ViTT`, `J-I-ViTT` | 562.9 | 396.9–798.8 | 37 | 0.471 |
| all other rows | — | — | — | 0 |

Waddington, Gibani naive, and every Maryland row stay at `CoP = 1`,
`CoP_sd_log10 = 0`. Darton enters per-subject (groups 6/7) and is untouched.

**Provenance note for the extract.** These σ are `[derived]`, not `[reported]`. Jin
tabulates GMT and CI, not sd. The derivation assumes a normal-theory CI on
$\log_{10}$ GM with the per-protocol $n$. Record this in `extracts/Jin_2017.md`
with the source table.

## 2. Left-censoring at the naive floor

`CoP = 1` **is** the assay LLD (3.7 EU/mL imputed), not an interior point. The
within-arm distribution is left-censored with an atom at 1, exactly as the Darton
individual data shows (64 of 90 subjects at exactly `CoP = 1.0`). Quadrature nodes
must therefore be floored:

$$\text{CoP}_m = \max\left(10^{\log_{10}\text{CoP}_j + \sigma_j z_m},\; 1\right)$$

This matters only for the Jin control arm ($\text{CoP} = 2.16$, $\sigma = 0.53$ —
a large fraction of nodes fall below the floor). Vi-PS and Vi-TT sit far enough
above LLD that flooring is inert.

**Known bias, not resolved here.** Jin's reported control GMT is itself computed on
LLD-imputed values, so σ = 0.526 is the sd of the *imputed* distribution and
understates the latent spread. Flooring the nodes is the right correction for the
observation model but does not undo that. An alternative for the control arm —
Jin reports 38% of controls had detectable anti-Vi — is an explicit two-component
mixture: 62% at `CoP = 1`, 38% lognormal above. `[option, not adopted]`

## 3. Stan changes

New data:

```stan
vector<lower=0>[N_obs] CoP_sd_log10;   // 0 => point mass at CoP (current behaviour)
int<lower=1> N_quad;                    // Gauss-Hermite nodes (9 is ample)
vector[N_quad] quad_z;                  // standard-normal nodes
simplex[N_quad] quad_w;                 // normalized weights
```

`obs_prob()` gains `CoP_sd` and the quadrature arrays, and groups 1 and 2 become
averages over nodes. Everything else in `obs_prob()` is unchanged.

**Implementation constraints:**

1. **Group 1 must integrate the product, not multiply the integrals.**
   $E[P_{\text{inf}} \cdot P_{\text{fev|inf}}] \neq E[P_{\text{inf}}] \cdot E[P_{\text{fev|inf}}]$.
   Both factors depend on the same subject's CoP.
2. **η factors out of group 2.** `eta_detection()` depends on dose, not CoP, so
   `eta * E[P_inf]` is correct and cheaper than integrating the product.
3. **Do not integrate φ(T, D).** `phi_TD()` evaluates the naive (CoP = 1) fever
   curve by construction. Immunity already acts upstream, so integrating this
   term would count it twice (`typhoid_dose_response.stan:139`). The function has
   no CoP argument. Groups 3/5 are Maryland and are unaffected.

Groups 3, 4, 5, 6, 7 take `CoP_sd = 0` and must be provably unchanged.

## 4. The three-implementation invariant

`obs_prob()` exists three times by deliberate design
(`test_obs_prob_parity.R` header). All three change together:

1. `calibration/typhoid_dose_response.stan` — `obs_prob()`
2. `calibration/model_math.R` — `mm_obs_prob()` (what the figures run)
3. `calibration/test_obs_prob_parity.R` — `obs_prob_R()` (the independent
   hand transcription; **do not** refactor it to call `model_math.R`)

## 5. Verification plan

Mirrors `tier2_plan.md` §5.

1. **Bit-identity gate (blocking).** With `CoP_sd_log10 = 0` on every row, refit
   `t2-indiv-vax` and require the posterior to match the committed run to sampler
   tolerance. This checks that quadrature is a strict generalization and that no
   existing tier moved.
2. **Parity gate.** Extend `test_obs_prob_parity.R` with σ > 0 covariate rows for
   groups 1 and 2, at several fixed parameter vectors. Stan `p_pred`,
   `mm_obs_prob()`, and `obs_prob_R()` must agree to `TOL_P = 1e-7`.
3. **Quadrature convergence.** `N_quad` ∈ {5, 9, 15, 31} must agree to < 1e-6 on
   every affected row at the posterior median. If 9 nodes is not enough, inspect
   the integrand or increase `N_quad`.
4. **Analytic floor check.** At σ = 0 the floored node set collapses to
   `{max(CoP, 1)}`; assert this in the test rather than trusting it.
5. **Refit `t2-indiv-vax` with the derived σ.** Report the delta on `gamma_inf`,
   `gamma_fevginf`, and the six Jin fitted probabilities against the committed run.

## 6. Falsifiable prediction

**Stated before the fit:** `gamma_inf` and `gamma_fevginf` move by less than
0.005 (well inside MCSE), and the six Jin fitted probabilities move by less than
0.02 absolute, with the control arm's fever row moving most (predicted ratio
0.978 at the posterior median, §5.3 table).

**If either γ moves by more than 0.02, treat this as an implementation failure**
until the group-1 product integration (§3 item 1) and row assignments have been
checked and gate 1 has been rerun. Do not report the movement as a finding before
that check.

## 7. Open questions before implementation

1. Is the normal-theory σ derivation acceptable, or do you want the control arm's
   two-component alternative (§2) in from the start?
2. `N_quad = 9` fixed, or exposed as a stage token?
3. Does this get its own tier key (e.g. `t2-indiv-vax-q`) for A/B against the
   committed run, or does it replace `t2-indiv-vax` in place after gate 1 passes?
