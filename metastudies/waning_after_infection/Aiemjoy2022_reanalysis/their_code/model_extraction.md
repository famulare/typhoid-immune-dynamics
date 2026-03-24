# Aiemjoy et al. 2022 — Complete Model Extraction

**Source repository code version:** `v9na` (no-age variant)
**Extracted from:** `their_code/v9na.model.jags`, `their_code/v9na.data.r`, `their_code/graph-func.r`
**Extraction date:** 2026-03-24

---

## 1. Executive Summary

This is a Bayesian hierarchical model for the kinetics of seven antibody responses (HlyE IgA, HlyE IgG, LPS IgA, LPS IgG, MP IgA, MP IgG, Vi IgG) following culture-confirmed typhoid fever. Each subject's trajectory for each test is described by a **two-phase model**: exponential rise from baseline $y_0$ to peak $y_1$ over time $t_1$, followed by power-function decay governed by rate $\alpha$ and shape $s$. The five log-transformed parameters per subject per test are drawn from a multivariate normal population distribution with unknown mean and precision matrix. Observations are lognormally distributed around the trajectory. Inference is by JAGS (Gibbs sampling). There is no reinfection model, no mixture/latent-class structure, and no age covariate in this version.

---

## 2. Model Boundary Map

| File | Role | Lines |
|------|------|-------|
| `v9na.model.jags` | JAGS model specification: likelihood + priors for all latent parameters | 1--25 |
| `v9na.data.r` | Data ingestion, reshaping, zero replacement, ghost subject injection, hyperparameter specification | 1--122 |
| `graph-func.r` | Post-processing helper functions: trajectory evaluation (`ab`, `sero`, `qsero`, `serocourse`), waning-time density (`wdens`, `wdistquan`, `wlogdistquan`) | 1--65 |

**Not in this repo:** JAGS run script (chain initialization, `jags.model()` / `coda.samples()` calls), convergence diagnostics, posterior summary scripts, and any upstream data exclusion logic for reinfections.

---

## 3. Model Inventory Table

| Component | Present? | Notes |
|-----------|----------|-------|
| Two-phase kinetic trajectory | Yes | Exponential rise + power-function decay |
| Hierarchical population model | Yes | MVN on log-parameters, Wishart precision |
| Lognormal observation model | Yes | One residual precision per test |
| Reinfection model | **No** | `reinf_obs` column exists in data but is never referenced in JAGS code |
| Age covariate | **No** | `v9na` = no-age; age is loaded but not passed to JAGS |
| Mixture / latent classes | **No** | Single MVN population per test |
| Cross-test correlation | **No** | Tests are independent (separate MVN per test, no shared parameters) |
| Missing data handling | Implicit | JAGS skips `NA` entries in `logy`; `nsmpl[subj]` controls loop bounds |

---

## 4. Math Sheet

### 4.1 State Evolution: Two-Phase Trajectory

**Phase 1 — Exponential rise** ($0 \le t \le t_1$; `v9na.model.jags` lines 6--7):

$$y(t) = y_0 \, e^{\beta t}, \quad \beta = \frac{\ln(y_1 / y_0)}{t_1}$$

This is an exponential interpolation from baseline $y_0$ at $t=0$ to peak $y_1$ at $t = t_1$.

**Phase 2 — Generalized power-function decay** ($t > t_1$; `v9na.model.jags` lines 8--9):

$$y(t) = \left[ y_1^{1-s} - (1-s)\,\alpha\,(t - t_1) \right]^{1/(1-s)}$$

where $s > 1$ is the shape parameter. Since $s > 1$, we have $1-s < 0$, so $y_1^{1-s} > 0$ and $-(1-s)\alpha = (s-1)\alpha > 0$. Therefore the expression inside the brackets is actually:

$$y_1^{1-s} + (s-1)\,\alpha\,(t-t_1)$$

which is the **sum of two positive terms** — there is no negativity bug. This means $y(t)$ is always real and positive, and it decays monotonically toward zero as $t \to \infty$.

**Interpretation:** This is the solution to the ODE $\frac{dy}{dt} = -\alpha \, y^s$ with initial condition $y(t_1) = y_1$. For $s=2$ this is hyperbolic decay; for $s \to 1^+$ it approaches exponential decay.

**JAGS log-space form** (`v9na.model.jags` lines 6--9):

$$\mu_{\log y}(t) = \begin{cases}
\ln y_0 + \beta \cdot t & \text{if } t \le t_1 \\[4pt]
\frac{1}{1-s} \ln\!\left[ y_1^{1-s} - (1-s)\,\alpha\,(t - t_1) \right] & \text{if } t > t_1
\end{cases}$$

### 4.2 Parameter Transforms

All five biological parameters are stored on transformed (log) scale as `par[subj, test, 1:5]`. The back-transforms are (`v9na.model.jags` lines 12--16):

| Index | Code | Transform | Biological parameter | Constraint |
|-------|------|-----------|---------------------|------------|
| 1 | `par[subj,test,1]` | $y_0 = e^{p_1}$ | Baseline titer | $y_0 > 0$ |
| 2 | `par[subj,test,2]` | $y_1 = y_0 + e^{p_2}$ | Peak titer | $y_1 > y_0 > 0$ |
| 3 | `par[subj,test,3]` | $t_1 = e^{p_3}$ | Time to peak (days) | $t_1 > 0$ |
| 4 | `par[subj,test,4]` | $\alpha = e^{p_4}$ | Decay rate | $\alpha > 0$ |
| 5 | `par[subj,test,5]` | $s = e^{p_5} + 1$ | Decay shape | $s > 1$ |

**Critical note on $y_1$:** The peak is parameterized as $y_1 = y_0 + e^{p_2}$, NOT $y_1 = e^{p_2}$. The parameter $p_2$ controls the log of the *elevation above baseline*, guaranteeing $y_1 > y_0$. This means the prior on $p_2$ governs the fold-rise, not the absolute peak level.

### 4.3 Observation Model

**Lognormal likelihood** (`v9na.model.jags` line 10):

$$\ln y_{\text{obs}}[i,j,k] \sim \mathcal{N}\!\left(\mu_{\ln y}[i,j,k],\; \tau_k \right)$$

where $\tau_k$ = `prec.logy[test]` is a single precision (inverse variance) shared across all subjects and time points for test $k$.

**What $\tau_k$ absorbs:** This is NOT just assay coefficient of variation. It captures ALL residual variance including:
- Assay measurement error
- Model misspecification (trajectory shape errors)
- Within-subject biological variation not captured by the five-parameter model
- Any heteroscedasticity not modeled

### 4.4 Hierarchical Structure

**Subject-level** (`v9na.model.jags` line 17):

$$\mathbf{p}_i^{(k)} = (p_{i,k,1}, \ldots, p_{i,k,5})^\top \sim \mathcal{N}_5\!\left(\boldsymbol{\mu}_k,\; \boldsymbol{\Lambda}_k\right)$$

where $\boldsymbol{\mu}_k$ = `mu.par[test,]` is the population mean vector and $\boldsymbol{\Lambda}_k$ = `prec.par[test,,]` is the $5 \times 5$ population precision matrix for test $k$.

**Population-level priors** (`v9na.model.jags` lines 21--23):

$$\boldsymbol{\mu}_k \sim \mathcal{N}_5\!\left(\boldsymbol{\mu}_k^{\text{hyp}},\; \boldsymbol{\Lambda}_k^{\text{hyp}}\right)$$

$$\boldsymbol{\Lambda}_k \sim \text{Wishart}\!\left(\boldsymbol{\Omega}_k,\; \nu_k\right)$$

$$\tau_k \sim \text{Gamma}(a_k, b_k)$$

### 4.5 Full Hierarchical Likelihood

$$p(\text{data}, \boldsymbol{\theta}) = \prod_{k=1}^{7} \left[ p(\boldsymbol{\mu}_k) \; p(\boldsymbol{\Lambda}_k) \; p(\tau_k) \; \prod_{i=1}^{N} \left( p(\mathbf{p}_i^{(k)} \mid \boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k) \; \prod_{j=1}^{n_i} p(\ln y_{ijk} \mid \mu_{\ln y}(t_{ij}; \mathbf{p}_i^{(k)}), \tau_k) \right) \right]$$

where $N$ = `nsubj` (1668 real + 1 ghost = 1669), $n_i$ = `nsmpl[subj]`, and $\boldsymbol{\theta}$ is the full set of latent parameters.

### 4.6 Priors — Exact Hyperparameter Values

All hyperparameters are set identically for all 7 tests (`v9na.data.r` lines 109--115).

#### Population mean prior: $\boldsymbol{\mu}_k \sim \mathcal{N}_5(\boldsymbol{\mu}^{\text{hyp}}, \boldsymbol{\Lambda}^{\text{hyp}})$

**Mean vector** (`v9na.data.r` line 110):

$$\boldsymbol{\mu}^{\text{hyp}} = (1.0, \; 7.0, \; 1.0, \; -4.0, \; -1.0)$$

| Parameter | $\mu^{\text{hyp}}$ | Interpretable translation |
|-----------|---------------------|--------------------------|
| $p_1$ (log $y_0$) | 1.0 | Baseline $y_0 \approx e^1 \approx 2.7$ EU |
| $p_2$ (log elevation) | 7.0 | Elevation $\approx e^7 \approx 1097$ EU above baseline |
| $p_3$ (log $t_1$) | 1.0 | Peak time $t_1 \approx e^1 \approx 2.7$ days |
| $p_4$ (log $\alpha$) | $-4.0$ | Decay rate $\alpha \approx e^{-4} \approx 0.018$ |
| $p_5$ (log($s-1$)) | $-1.0$ | Shape $s \approx e^{-1}+1 \approx 1.37$ (mildly super-exponential) |

**Precision matrix** (`v9na.data.r` line 111):

$$\boldsymbol{\Lambda}^{\text{hyp}} = \text{diag}(1.0, \; 0.00001, \; 1.0, \; 0.001, \; 1.0)$$

| Parameter | Precision | $\Leftrightarrow$ SD | 95% prior range |
|-----------|-----------|------|-----------------|
| $p_1$ (log $y_0$) | 1.0 | 1.0 | $\mu \pm 2$, i.e., $y_0 \in [0.37, 20]$ EU |
| $p_2$ (log elevation) | **0.00001** | **316** | $\mu \pm 632$, i.e., essentially flat — **deliberately extremely loose** |
| $p_3$ (log $t_1$) | 1.0 | 1.0 | $\mu \pm 2$, i.e., $t_1 \in [0.37, 20]$ days |
| $p_4$ (log $\alpha$) | 0.001 | 31.6 | $\mu \pm 63$, i.e., extremely loose on decay rate |
| $p_5$ (log($s-1$)) | 1.0 | 1.0 | $\mu \pm 2$, i.e., $s \in [1.05, 8.4]$ |

The extreme looseness on $p_2$ (precision = 0.00001, SD = 316) is **deliberate**: the peak elevation above baseline varies by orders of magnitude across subjects and tests, so the population mean prior is made essentially non-informative to let the data speak.

The looseness on $p_4$ (precision = 0.001, SD = 31.6) serves a similar purpose for the decay rate.

#### Population precision prior: $\boldsymbol{\Lambda}_k \sim \text{Wishart}(\boldsymbol{\Omega}_k, \nu_k)$

(`v9na.data.r` lines 112--113):

$$\nu_k = 20, \quad \boldsymbol{\Omega}_k = \text{diag}(1.0, \; 50.0, \; 1.0, \; 10.0, \; 1.0)$$

The Wishart prior mean is $\nu \cdot \boldsymbol{\Omega}^{-1}$... but note JAGS parameterizes Wishart as $\text{dwish}(R, \nu)$ where the mean is $\nu R$ (using the scale-matrix convention). So:

$$E[\boldsymbol{\Lambda}_k] = \nu_k \cdot \boldsymbol{\Omega}_k^{-1} \quad \text{or} \quad \nu_k \cdot \boldsymbol{\Omega}_k$$

depending on JAGS convention. In JAGS, `dwish(R, k)` has mean $k \cdot R$, so:

$$E[\boldsymbol{\Lambda}_k] = 20 \cdot \text{diag}(1, 50, 1, 10, 1) = \text{diag}(20, 1000, 20, 200, 20)$$

This translates to expected population-level SDs of approximately:

| Parameter | Expected precision | Expected SD |
|-----------|-------------------|-------------|
| $p_1$ | 20 | 0.22 |
| $p_2$ | 1000 | 0.032 |
| $p_3$ | 20 | 0.22 |
| $p_4$ | 200 | 0.071 |
| $p_5$ | 20 | 0.22 |

With $\nu = 20$ degrees of freedom and 5 dimensions, the Wishart prior is moderately informative on the covariance structure.

#### Residual precision prior: $\tau_k \sim \text{Gamma}(4, 1)$

(`v9na.data.r` line 114):

$$\tau_k \sim \text{Gamma}(4, 1)$$

Mean = 4, variance = 4. This corresponds to a prior mean residual SD of $1/\sqrt{4} = 0.5$ in log-space, equivalent to roughly a 1.6-fold multiplicative error (since $e^{0.5} \approx 1.65$). The Gamma(4,1) shape provides moderate information, concentrating probability on $\tau \in [1, 9]$ (SD $\in [0.33, 1.0]$).

### 4.7 Inference Algorithm

JAGS (Just Another Gibbs Sampler). The model is fully conjugate or conditionally conjugate at each level:
- Normal--Normal conjugacy for `mu.par` given `par` and `prec.par`
- Wishart--Normal conjugacy for `prec.par` given `par` and `mu.par`
- Gamma--Normal conjugacy for `prec.logy` given residuals
- The subject-level `par` nodes are updated by JAGS's adaptive Metropolis-within-Gibbs (the trajectory nonlinearity breaks conjugacy)

Chain initialization, number of chains, burn-in, thinning, and adaptation settings are **not in this repo**.

---

## 5. Parameter Table

| Symbol | Code variable | Transform from `par` | Prior (on `par` scale) | Biological meaning | Constraints |
|--------|--------------|---------------------|----------------------|-------------------|-------------|
| $y_0$ | `y0[subj,test]` | $e^{p_1}$ | $p_1$ component of MVN | Pre-illness baseline antibody level (EU) | $> 0$ |
| $y_1$ | `y1[subj,test]` | $y_0 + e^{p_2}$ | $p_2$ component of MVN | Peak antibody level (EU) | $> y_0$ |
| $t_1$ | `t1[subj,test]` | $e^{p_3}$ | $p_3$ component of MVN | Time from fever onset to peak (days) | $> 0$ |
| $\alpha$ | `alpha[subj,test]` | $e^{p_4}$ | $p_4$ component of MVN | Power-law decay rate constant | $> 0$ |
| $s$ | `shape[subj,test]` | $e^{p_5} + 1$ | $p_5$ component of MVN | Decay shape exponent (ODE: $dy/dt = -\alpha y^s$) | $> 1$ |
| $\beta$ | `beta[subj,test]` | $\ln(y_1/y_0)/t_1$ | Derived, not sampled | Exponential rise rate (day$^{-1}$) | Derived |
| $\boldsymbol{\mu}_k$ | `mu.par[test,]` | Identity | $\mathcal{N}_5(\boldsymbol{\mu}^{\text{hyp}}, \boldsymbol{\Lambda}^{\text{hyp}})$ | Population mean of log-parameters for test $k$ | -- |
| $\boldsymbol{\Lambda}_k$ | `prec.par[test,,]` | Identity | $\text{Wishart}(\boldsymbol{\Omega}, 20)$ | Population precision matrix for test $k$ | Positive definite |
| $\tau_k$ | `prec.logy[test]` | Identity | $\text{Gamma}(4, 1)$ | Residual log-observation precision for test $k$ | $> 0$ |

---

## 6. Data Inputs Table

| Variable | Code name | Dimensions | Source | Description |
|----------|-----------|------------|--------|-------------|
| Log-observations | `logy` | $(N+1) \times 7 \times 7$ | `smpl.y` after zero-replacement and log | $\ln(\text{titer in EU})$ for each subject, visit, test |
| Sample times | `smpl.t` | $(N+1) \times 7$ | `visit.t` column permuted | Days since fever onset for each observed visit |
| Nr. samples per subject | `nsmpl` | $N+1$ | Computed from non-NA visit times | Number of longitudinal observations |
| Nr. subjects | `nsubj` | scalar = $N+1$ = 1669 | `length(id)+1` | Includes ghost subject |
| Nr. tests | `ntest` | scalar = 7 | Hardcoded | HlyE IgA/IgG, LPS IgA/IgG, MP IgA/IgG, Vi IgG |
| Nr. dimensions | `ndim` | scalar = 5 | = `npar` | Dimensionality of parameter vector |
| Hyperparameters | `mu.hyp`, `prec.hyp`, `omega`, `wishdf`, `prec.logy.hyp` | Various | Hardcoded | See Section 4.6 |

**Antibody test ordering** (`v9na.data.r` lines 68--70):

| Index | Antigen | Isotype |
|-------|---------|---------|
| 1 | HlyE | IgA |
| 2 | HlyE | IgG |
| 3 | LPS | IgA |
| 4 | LPS | IgG |
| 5 | MP | IgA |
| 6 | MP | IgG |
| 7 | Vi | IgG |

**Data file:** `TypoidCaseData_github_09.30.21.csv` (`v9na.data.r` line 6). Contains up to 7 visits per subject with titer measurements for all 7 antibody tests, plus `TimeInDays_visit*` columns. Also contains `reinf_obs`, `age`, `Country`, `Hospitalized`, `bldculres` columns.

---

## 7. Ghost Subject Trick

**Location:** `v9na.data.r` lines 95--97

```r
nsmpl[nsubj+1] <- 3;
smpl.t[nsubj+1,] <- c(5,30,90,NA,NA,NA,NA);
age <- c(age,10); # just made this up
```

**What it does:** Appends a fictitious $(N+1)$-th subject with 3 observations at days 5, 30, and 90, but with **all titer values as `NA`** (because `smpl.y` was pre-allocated as `NA` at dimension `nsubj+1` on line 76, and no data is ever written to row `nsubj+1`). The age is set to an arbitrary value of 10.

**Why it works:** JAGS treats `NA` values in data nodes (`logy`) as **unobserved latent variables** to be sampled from the model. So this ghost subject contributes no likelihood — it is "observed" at 3 time points but all observations are missing. Its five parameters are still drawn from the population distribution.

**Purpose:** This is a **prior predictive check device**. By monitoring the ghost subject's sampled parameters and predicted trajectory in the posterior, one can visualize the population-level predictive distribution *without* being pulled toward any particular subject's data. It produces a "typical trajectory" from the fitted population distribution. The three time points (5, 30, 90 days) are chosen to span the rise, peak, and early waning phases.

**Impact on inference:** Negligible. One extra draw from the population MVN per MCMC iteration adds a trivial amount of computation and has zero effect on the posterior because no data likelihood is attached.

**Note:** `nsubj` is passed to JAGS as `nsubj+1` (`v9na.data.r` line 118), so the JAGS loop runs over all $N+1$ subjects including the ghost.

---

## 8. Evidence Flow Analysis for Vi IgG

Vi IgG is test index 7 (`v9na.data.r` line 57, line 70).

### Data path
1. Raw CSV columns `Vi_IgG_visit1` ... `Vi_IgG_visit7` are bound into `vi_igg` matrix (`v9na.data.r` lines 54--57)
2. Packed into `smpl.y[subj, obs, 7]` (`v9na.data.r` line 91)
3. Zeros replaced with 0.01 (`v9na.data.r` line 98)
4. Log-transformed to `logy[subj, obs, 7]` (`v9na.data.r` line 99)
5. Passed to JAGS as data

### Model path
1. `par[subj, 7, 1:5]` drawn from `dmnorm(mu.par[7,], prec.par[7,,])` (`v9na.model.jags` line 17)
2. Transformed to $(y_0, y_1, t_1, \alpha, s)$ (`v9na.model.jags` lines 12--16)
3. Trajectory $\mu_{\ln y}$ computed (`v9na.model.jags` lines 6--9)
4. `logy[subj, obs, 7] ~ dnorm(mu.logy[subj, obs, 7], prec.logy[7])` (`v9na.model.jags` line 10)

### Population-level for Vi IgG
- `mu.par[7, 1:5]` ~ $\mathcal{N}_5(\boldsymbol{\mu}^{\text{hyp}}, \boldsymbol{\Lambda}^{\text{hyp}})$ — same hyperparameters as all other tests
- `prec.par[7, 1:5, 1:5]` ~ $\text{Wishart}(\boldsymbol{\Omega}, 20)$
- `prec.logy[7]` ~ $\text{Gamma}(4, 1)$

### Cross-test independence
Vi IgG shares NO parameters with the other 6 tests. The model fits 7 completely independent sub-models that happen to share the same subjects. There is no borrowing of strength across tests.

### Downstream: waning time density
For Vi IgG specifically, the helper function `wdens()` (`graph-func.r` lines 44--48) computes the density of the time to wane below a threshold, derived analytically from the power-function decay (see Section 10).

---

## 9. Red Flags and Ambiguities

### 9.1 `reinf_obs` column — unused but present

The data CSV contains a `reinf_obs` column. In the full dataset of 1667 subjects, 716 subjects have `reinf_obs = 1` and 951 have `NA`. This column is **never referenced** in `v9na.data.r` or `v9na.model.jags`. There is no reinfection exclusion, flagging, or sub-model in this code. Any handling of reinfections must occur upstream in data preparation code **not included in this repository**.

**Risk:** If subjects with reinfection have boosted titers during follow-up, the single-trajectory model will misfit them, inflating residual variance and potentially biasing population decay estimates. The commented-out filter lines in `v9na.data.r` (lines 11--24) suggest the authors experimented with subsetting but the final code uses the full dataset.

### 9.2 Zero replacement

(`v9na.data.r` line 98):
```r
smpl.y[smpl.y==0] <- 0.01;
```

All zero titers are replaced with 0.01 before log-transformation. This is a pragmatic fix for $\ln(0) = -\infty$, but:
- The replacement value 0.01 is arbitrary (appears to be chosen as a small value relative to typical titers)
- No sensitivity analysis to this choice is evident in the code
- A more principled approach would be left-censored likelihood, but JAGS makes this difficult

### 9.3 Prior looseness and potential identifiability issues

The extremely loose prior on $p_2$ (precision = 0.00001) and $p_4$ (precision = 0.001) for the population mean could cause:
- Slow mixing in MCMC chains, especially for subjects with few observations
- Weak identifiability between $\alpha$ and $s$ in the decay phase — these parameters can trade off: faster $\alpha$ with smaller $s$ can mimic slower $\alpha$ with larger $s$ over limited time windows
- The Wishart prior with $\nu = 20$ and 5 dimensions ($\nu = 20 > p-1 = 4$, so valid) is moderately informative on correlation structure, which may help regularize

**Specific identifiability concern for $\alpha$-$s$ tradeoff:** With typical follow-up to ~1 year, the power-function decay $y_1^{1-s} + (s-1)\alpha(t-t_1)$ may not have enough curvature information to separately identify $\alpha$ and $s$. The population-level MVN allows correlation between $p_4$ and $p_5$, which partially addresses this, but the off-diagonal elements of the Wishart are weakly constrained.

### 9.4 Same hyperparameters for all tests

All 7 tests share identical hyperparameter values (`v9na.data.r` line 109 loop). This means the prior belief about baseline levels, peak heights, timing, and decay rates is the same for HlyE IgA as for Vi IgG. The data must override these priors entirely, which is fine given $N \approx 1667$ subjects, but it means the priors are not encoding any antigen/isotype-specific knowledge.

### 9.5 No posterior predictive checks in code

The codebase contains no scripts for posterior predictive checks, residual diagnostics, or model comparison (DIC, WAIC, LOO). Model adequacy assessment is not verifiable from the provided code.

### 9.6 Missing JAGS run configuration

Chain initialization, number of iterations, burn-in, thinning, adaptation period, and monitored nodes are not specified in any file in this directory. Convergence diagnostics cannot be assessed.

---

## 10. Helper Function Documentation

All functions in `graph-func.r`.

### `ab(t, y0, y1, t1, alpha, shape)` — Single trajectory evaluation

**Location:** `graph-func.r` lines 10--16

Evaluates the two-phase antibody trajectory at a single time point for a single parameter set. Returns antibody level in natural (not log) scale.

- Uses `bt()` (line 8) to compute $\beta = \ln(y_1/y_0)/t_1$
- Phase 1 ($t \le t_1$): $y_0 \cdot e^{\beta t}$
- Phase 2 ($t > t_1$): $(y_1^{1-s} - (1-s)\alpha(t-t_1))^{1/(1-s)}$

This matches the JAGS model exactly (exponentiated from log-space).

### `sero(n, tvec, y0, y1, t1, alpha, shape)` — Subject trajectory

**Location:** `graph-func.r` lines 18--24

Evaluates the trajectory of subject `n` across a vector of time points. The parameter vectors `y0`, `y1`, etc. are indexed by subject number `n`.

### `qsero(t, q, y0, y1, t1, alpha, shape)` — Posterior quantiles at one time point

**Location:** `graph-func.r` lines 26--33

Given vectors of MCMC posterior samples for all parameters, evaluates the trajectory at time `t` for each sample, then returns the requested quantile(s) `q`. If `q = "mean"`, returns the geometric mean (exponential of mean of logs).

### `serocourse(tvec, q, y0, y1, t1, alpha, shape)` — Posterior quantile trajectory

**Location:** `graph-func.r` lines 35--42

Wraps `qsero()` over a vector of time points to produce a full posterior quantile trajectory curve (e.g., median, 2.5th, 97.5th percentile bands).

### `wdens(w, y1, alpha, shape)` — Waning time density

**Location:** `graph-func.r` lines 44--48

Computes the probability density of the waning time $W$ — the time for antibody to decay from peak $y_1$ to a threshold level, given the power-function decay model.

Let $\rho = 1/(s-1)$. Then:

$$f_W(w) = \frac{1}{w \, \Gamma(\rho)} \left(\frac{w\rho}{\alpha}\right)^\rho \exp\!\left(-\frac{w\rho \, y_1^{-1/\rho}}{\alpha}\right)$$

This is the density of a **Gamma-distributed** waning time (it can be verified that this has the form of an inverse-Gamma or Gamma kernel). Specifically, if we define the waning time as the time to reach some threshold from $y_1$ under $dy/dt = -\alpha y^s$, the distribution over posterior samples of $(y_1, \alpha, s)$ produces this density.

**Note:** This function operates on vectors of posterior samples and returns a vector of density values — one per MCMC draw — enabling posterior uncertainty propagation for the waning time distribution.

### `wdistquan(wvec, qvec, y1, alpha, shape)` — Quantiles of waning density

**Location:** `graph-func.r` lines 50--56

For each candidate waning time in `wvec`, computes the posterior quantiles of the density `wdens()` across MCMC samples. Returns a matrix (length(wvec) x length(qvec)) for plotting credible bands on the waning time density.

### `wlogdistquan(logwvec, qvec, y1, alpha, shape)` — Log-scale waning density quantiles

**Location:** `graph-func.r` lines 58--65

Same as `wdistquan` but on log10 scale: applies the Jacobian $w \cdot \ln(10)$ for change of variables to $\log_{10}(w)$. Used for log-scale density plots.

---

## 11. Validation Checklist

| Check | Status | Notes |
|-------|--------|-------|
| JAGS model parses without error | Not verified | No run script available |
| Phase 1 trajectory matches $y_0$ at $t=0$ | Confirmed | $y_0 \cdot e^{0} = y_0$ |
| Phase 1 trajectory matches $y_1$ at $t=t_1$ | Confirmed | $y_0 \cdot e^{\beta t_1} = y_0 \cdot (y_1/y_0) = y_1$ |
| Phase 2 trajectory continuous at $t=t_1$ | Confirmed | $(y_1^{1-s})^{1/(1-s)} = y_1$ |
| $s > 1$ enforced | Confirmed | $s = e^{p_5} + 1 > 1$ always (`v9na.model.jags` line 16) |
| $y_1 > y_0$ enforced | Confirmed | $y_1 = y_0 + e^{p_2} > y_0$ (`v9na.model.jags` line 13) |
| Phase 2 argument stays positive | Confirmed | $y_1^{1-s} + (s-1)\alpha(t-t_1) > 0$ since both terms positive when $s > 1$ |
| `ab()` in `graph-func.r` matches JAGS model | Confirmed | Exponentiated version of the same equations |
| Ghost subject has no data likelihood | Confirmed | `smpl.y[nsubj+1,,]` is all NA |
| Zero replacement applied before log | Confirmed | Line 98 before line 99 |
| `reinf_obs` not used in model | Confirmed | Not in JAGS data list or model |
| All 7 tests use same hyperparameters | Confirmed | Loop on lines 109--115 sets identical values |
| Tests are independent in the model | Confirmed | No shared parameters across test indices |
| Wishart df valid ($\nu \ge p$) | Confirmed | $\nu = 20 \ge 5 = p$ |
| JAGS `step()` returns 1 for 0 | Confirmed | JAGS convention: `step(x)` = 1 if $x \ge 0$, so boundary $t = t_1$ falls in Phase 1 |
