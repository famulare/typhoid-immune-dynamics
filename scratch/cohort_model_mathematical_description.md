---
output:
  html_document: default
  pdf_document: default
---
# Typhoid Cohort Incidence Model: Mathematical Description

## Overview

This model simulates typhoid infection and disease dynamics over individual lifetimes using anti-Vi IgG titers as a correlate of protection (CoP). The core biological insight is that immunity modulates susceptibility through dose-dependent mechanisms—higher bacterial doses can overwhelm prior immunity.

---

## 1. Titer Response Function

Describes CoP dynamics after an immunizing event (infection or vaccination). The function has two phases: an initial rise to peak, then power-law decay.

### 1.1 Age-Dependent Parameters

Both decay dynamics depend on age at infection:

$$T_{d,a} = T_{\text{decay}} \cdot e^{\beta_T \cdot \text{age}}$$

$$\kappa_a = \alpha \cdot e^{\beta_\alpha \cdot \text{age}}$$

where:
- $T_{\text{decay}} = 11$ days (reference decay timescale)
- $\alpha = 1.16$ (reference power-law exponent)
- $\beta_T = 0.057$ (age coefficient for decay time; positive → slower decay with age)
- $\beta_\alpha = -0.060$ (age coefficient for exponent; negative → different decay shape with age)

### 1.2 Rise Phase ($t < t_{\text{peak}}$)

A normalized logistic function interpolates from pre-challenge titer to peak:

$$
\text{titer}(t) = \text{CoP}_{\text{pre}} + (\text{CoP}_{\text{peak}} - \text{CoP}_{\text{pre}}) \cdot \frac{L(t) - L(0)}{L(t_{\text{peak}}) - L(0)}
$$

where the logistic kernel is:

$$
L(t) = \frac{1}{1 + e^{-(t - 5 t_{\text{start}})/T_{\text{rise}}}}
$$

Default parameters:
- $t_{\text{peak}} = 30.4$ days (~1 month to peak)
- $t_{\text{start}} = 3$ days (inflection point scaling)
- $T_{\text{rise}} = 1.5$ days (rise steepness)

The normalization ensures:
- At $t = 0$: $\text{titer} = \text{CoP}_{\text{pre}}$
- At $t = t_{\text{peak}}$: $\text{titer} = \text{CoP}_{\text{peak}}$

**Note:** The rise phase is primarily cosmetic for continuity—the monthly timestep means susceptibility calculations use the titer value at the start of each month, so intra-month dynamics don't affect outcomes.

### 1.3 Decay Phase ($t \geq t_{\text{peak}}$)

Power-law decay from peak back toward baseline:

$$
\text{titer}(t) = \text{CoP}_{\min} + (\text{CoP}_{\text{peak}} - \text{CoP}_{\min}) \left(1 + \frac{t - t_{\text{peak}}}{\kappa_a \cdot T_{d,a}}\right)^{-\kappa_a}
$$

**Biological interpretation:**
- Power-law (not exponential) decay matches observed long-term antibody kinetics
- Older individuals have slower waning (larger $T_{d,a}$), consistent with TCV trial data
- The exponent $\kappa_a$ controls the "heaviness" of the tail

---

## 2. Fold-Rise Model

Determines the peak titer achieved after an immunizing event, given the pre-challenge titer.

$$
\text{fold-rise} = 10^{\mu \cdot \left(1 - \frac{\log_{10}(\text{CoP}_{\text{pre}}) - \log_{10}(\text{CoP}_{\min})}{\log_{10}(\text{CoP}_{\max}) - \log_{10}(\text{CoP}_{\min})}\right)}
$$

For individual-level variation: $\mu \sim \mathcal{N}(\mu_0, \sigma_0)$ with $\mu_0 = 1.25$, $\sigma_0 = 0.5$.

The post-challenge peak is then:
$$
\text{CoP}_{\text{peak}} = \text{CoP}_{\text{pre}} \times \text{fold-rise}
$$

**Key properties:**
- When $\text{CoP}_{\text{pre}} = \text{CoP}_{\min}$: fold-rise $= 10^{\mu_0} \approx 18$
- When $\text{CoP}_{\text{pre}} = \text{CoP}_{\max}$: fold-rise $= 10^0 = 1$ (ceiling reached)
- Default ceiling: $\text{CoP}_{\max} = 10^{3.5} \approx 3162$ EU/ml

**Biological interpretation:** Naive individuals mount large responses; those near the antibody ceiling see minimal boosting. This enforces a physiological maximum titer.

---

## 3. Titer Timeseries Composition Over Multiple Infections

The lifetime titer trajectory is built as a **piecewise function** where each infection event overwrites the future trajectory.

### 3.1 Initialization

At birth ($t = 0$), every individual starts naive:
$$
\text{titer}(t) = \text{CoP}_{\min} = 1 \quad \forall t \in [0, t_{\max}]
$$

### 3.2 Upon Infection at Time $t_k$

When an infection occurs at month $t_k$:

**Step 1: Read current state**
$$
\text{CoP}_{\text{pre}}^{(k)} = \text{titer}(t_k)
$$

**Step 2: Calculate new peak via fold-rise**
$$
\text{CoP}_{\text{peak}}^{(k)} = \text{CoP}_{\text{pre}}^{(k)} \times \text{fold-rise}(\text{CoP}_{\text{pre}}^{(k)})
$$

**Step 3: Overwrite future trajectory**

For all $t \geq t_k$, replace the titer trajectory:
$$
\text{titer}(t) \leftarrow \text{titer\_vs\_time}\left(t - t_k, \; \text{age} = t_k/12, \; \text{CoP}_{\text{pre}}^{(k)}, \; \text{CoP}_{\text{peak}}^{(k)}\right)
$$

The previous trajectory for $t \geq t_k$ is **discarded**—the new response completely supersedes it.

### 3.3 Worked Example: Two Infections

Consider an individual infected at age 2 years ($t_1 = 24$ months) and again at age 7 years ($t_2 = 84$ months):

**Initial state (birth to infection 1):**
$$
\text{titer}(t) = 1 \quad \text{for } t \in [0, 24)
$$

**After infection 1 at $t_1 = 24$:**
- $\text{CoP}_{\text{pre}}^{(1)} = 1$
- $\text{fold-rise}^{(1)} \approx 10^{1.25} \approx 18$ (for a naive person)
- $\text{CoP}_{\text{peak}}^{(1)} = 1 \times 18 = 18$
- For $t \geq 24$: titer follows rise-then-decay curve with peak 18

**After infection 2 at $t_2 = 84$:**
- $\text{CoP}_{\text{pre}}^{(2)} = \text{titer}(84)$ ← read from waned curve (say, ~5 EU/ml after 5 years of decay)
- $\text{fold-rise}^{(2)} \approx 10^{1.25 \cdot (1 - \log_{10}(5)/3.5)} \approx 10$ (smaller than first infection)
- $\text{CoP}_{\text{peak}}^{(2)} = 5 \times 10 = 50$
- For $t \geq 84$: entire future trajectory replaced with new curve from (5 → 50 → waning)

**Resulting piecewise trajectory:**
```
titer(t) =
  ⎧ 1                                           if t < 24
  ⎨ rise_decay(t-24, pre=1, peak=18, age=2)     if 24 ≤ t < 84
  ⎩ rise_decay(t-84, pre=5, peak=50, age=7)     if t ≥ 84
```

### 3.4 Key Implications

1. **Memory through state, not history:** The model doesn't track infection history explicitly—all past immunological experience is encoded in the current titer value.

2. **Boosting from elevated baseline:** Reinfection when titer is still elevated starts the new curve from that elevated baseline, leading to higher peaks (until ceiling is approached).

3. **Age affects decay, not peak:** The fold-rise model is age-independent; only the subsequent waning rate depends on age at infection.

4. **Complete replacement:** There's no superposition of immune responses—each new infection fully determines the future trajectory.

---

## 4. Dose-Response Model

Probability of infection or fever given bacterial dose and immunity. Uses a modified beta-Poisson form where immunity scales the effective number of "hits."

### 4.1 Infection Given Dose

$$
P(\text{infection} | D, \text{CoP}) = 1 - \left(1 + D \cdot \frac{2^{1/\alpha_I} - 1}{N_{50,I}}\right)^{-\alpha_I / \text{CoP}^{\gamma_I}}
$$

Default parameters:
- $N_{50,I} = 695$ (dose for 50% infection probability in naive individuals)
- $\alpha_I = 1.68$ (shape parameter)
- $\gamma_I = 0.2$ (immunity scaling exponent)

### 4.2 Fever Given Dose

$$
P(\text{fever} | D, \text{CoP}) = 1 - \left(1 + D \cdot \frac{2^{1/\alpha_F} - 1}{N_{50,F}}\right)^{-\alpha_F / \text{CoP}^{\gamma_F}}
$$

Default parameters:
- $N_{50,F} = 27800$ (dose for 50% fever probability in naive individuals)
- $\alpha_F = 0.84$ (shape parameter)
- $\gamma_F = 0.4$ (immunity scaling exponent)

### 4.3 Fever Given Infection

Derived as conditional probability:
$$
P(\text{fever} | \text{infection}, D, \text{CoP}) = \frac{P(\text{fever} | D, \text{CoP})}{P(\text{infection} | D, \text{CoP})}
$$

### 4.4 Mechanism of Immunity

The CoP enters the exponent as $\alpha / \text{CoP}^{\gamma}$:
- When $\text{CoP} = 1$ (naive): exponent is $-\alpha$, giving standard beta-Poisson
- When $\text{CoP} > 1$: exponent magnitude decreases → probability decreases
- Higher $\gamma$ means immunity has stronger protective effect

**Biological interpretation:** Immunity reduces the effective "hit rate" of bacterial challenge. At very high doses, even high immunity can be overwhelmed—this is the key distinction from standard logistic CoP models.

---

## 5. Protective Efficacy

Relative risk reduction for an immune individual vs. naive control:

$$
\text{PE}(D, \text{CoP}) = 1 - \frac{P(\text{outcome} | D, \text{CoP})}{P(\text{outcome} | D, \text{CoP} = 1)}
$$

**Critical insight:** PE depends on challenge dose. The same titer provides different apparent protection in low-dose vs. high-dose settings—explaining why observed vaccine efficacy varies across epidemiological contexts.

---

## 6. Cohort Model Structure

Simulates $N$ individuals from birth to `age_max` years at monthly timesteps ($\Delta t = 30.4$ days).

### 6.1 Exposure Process

Exposures occur via Poisson process with age-dependent rate:
$$
n_{\text{exposures}}(t) \sim \text{Poisson}(\lambda \cdot m_{\text{age}}(t) \cdot \Delta t)
$$

where $m_{\text{age}}(t)$ is an age-dependent multiplier (default: reduced in infancy due to diet).

### 6.2 Simulation Algorithm

```
For each individual i = 1 to N:

    Initialize: titer[0:T] = 1  (naive at all future times)

    # Phase 1: Generate all exposure events
    For each month t:
        n_exp[t] ~ Poisson(λ · m_age(t))

    # Phase 2: Process exposures sequentially
    For each month t where n_exp[t] > 0:

        # Infection probability (multiple exposures compound)
        p_once = P(infection | dose, titer[t])
        p_inf = 1 - (1 - p_once)^n_exp[t]

        If rand() < p_inf:  # Infection occurs

            # Fever determination
            p_fever = P(fever | infection, dose, titer[t])
            fever[t] = (rand() < p_fever)

            # Update titer trajectory
            CoP_pre = titer[t]
            CoP_peak = CoP_pre × fold_rise(CoP_pre)

            For all s ≥ t:
                titer[s] = titer_vs_time(s - t, age=t/12, CoP_pre, CoP_peak)
```

### 6.3 Output Variables

| Variable | Description |
|----------|-------------|
| `titer[t]` | CoP value at each monthly timestep |
| `exposure_month` | Times when exposure events occurred |
| `exposure_count` | Number of exposures per event |
| `infected` | Binary: did exposure lead to infection? |
| `fever` | Binary: did infection lead to fever? |

---

## 7. Model Flow Diagram

```
Birth: titer(t) = 1 ∀t
         │
         ▼
┌────────────────────────────────────────────────────────┐
│  For each month t = 0, 1, ..., age_max × 12:           │
│                                                        │
│    n_exp ~ Poisson(λ · m_age(t))                       │
│                                                        │
│    If n_exp > 0:                                       │
│      │                                                 │
│      ├─► p_inf = 1 - (1 - P(inf|dose,titer[t]))^n_exp │
│      │                                                 │
│      ├─► If infected (with prob p_inf):               │
│      │     │                                           │
│      │     ├─► Determine fever (conditional on inf)    │
│      │     │                                           │
│      │     └─► Update future titers:                   │
│      │           CoP_pre = titer[t]                    │
│      │           CoP_peak = CoP_pre × fold_rise        │
│      │           titer[t:end] ← new waning curve       │
│      │                                                 │
│      └─► If not infected: titer unchanged              │
└────────────────────────────────────────────────────────┘
         │
         ▼
    Output: titer traces, infection/fever events by age
```

---

## 8. Key Biological Assumptions

1. **Anti-Vi IgG as sole CoP:** A single serological marker captures protective immunity
2. **Dose-dependent protection:** Efficacy is relative to challenge dose, not absolute
3. **Homogeneous immune kinetics:** All individuals respond identically (within age strata)
4. **Infection required for boosting:** Exposure without infection doesn't affect titers
5. **Complete trajectory replacement:** New infections fully supersede previous immune dynamics
6. **No maternal antibodies:** Individuals start naive at birth
7. **Constant transmission ecology:** Force of infection doesn't change over the simulation
8. **Monthly resolution:** Sub-monthly dynamics (incubation, acute phase) are not resolved
