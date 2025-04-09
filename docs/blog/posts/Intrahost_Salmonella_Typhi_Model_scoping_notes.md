---
draft: false
date:
  created: 2025-04-09
authors:
  - famulare
categories:
  - Immunity models
  - Intrahost models
  - Software architecture
  - Project Management
---

# Intra-host *Salmonella Typhi* Model - scoping notes

Notes prepared from voice dictation with assistance of ChatGPT-4o.

## Scope

These notes provide the skeleton of a model that aims to simulate within-host dynamics of *Salmonella Typhi*, connecting immune responses (especially antibody titers) to protection outcomes. It bridges vaccine efficacy and immunogenicity study field data and controlled human infection models (CHIM), and serves as the foundation for a next-generation Typhoid intrahost model suitable for generalizing across all exposure and immunological histories.

**Note that we may not ever build this model** as it is beyond the scope of current vaccine schedule policy needs. However, it is essential to think through the whole problem to understand how to make incisive tactical simplifications of the science, the project scope, and deliverables.


---

## Core Model Ingredients

### 1. Immune Components
- Serum antibody titers from vaccination (e.g. anti-Vi IgG) or natural exposure (e.g. anti-hemolysin-E IgG).
- Potential mucosal or cellular immunity (to be abstracted, as likely not directly modelable since not measured).
- Waning of immunity over time.

### 2. Host Outcomes
- Probability of infection upon exposure.
- Probability of developing disease given infection.
    - Multiple levels of disease severity?
- Probability of bacteremia given exposure
- Probablity of converse to carrier given exposure (likely mediated by bacteremia (+ host factors)). 

### 3. Vaccination and Exposure
- Vaccine type, dose, and schedule.
    - Primary focus on Typbar-TCV (Vi-TT)
    - also relevant is Vi-polysachharide, other conjugate VI (Vi-rEPA, non-Typbar Vi-TT)
    - Non-Vi-based vaccine (Ty21a) probably well out of scope but may get "for free" with natural infection model
- Time since vaccination or prior infection.
- Natural exposure characteristics (dose, infection frequency).

### 4. Data Sources
- Field efficacy trials.
- Cross-study immunogenicity and serological datasets.
- Controlled Human Infection Model (CHIM) studies.

---

## Key Complexities

### 1. Mixed Immunological Histories
- Real-world individuals often have overlapping natural exposure and vaccine histories.
- Hard to isolate clean correlates of protection.

### 2. Data Harmonization
- Assay variability and inconsistent outcome definitions across studies.
- Need for preprocessing and standardization.

### 3. Multi-scale Integration
- Bridging within-host protection dynamics with population-level outcomes.
- Eventually feeding into stochastic cohort model and transmission model.

### 4. Parameter Identifiability
- Difficult to separate effects of different immune components from limited data.
- Identify data sources for each model component.
    - Current gap to be filled in: boosting.
    - Nice to have: dose response vs immunity via challenge studies

### 5. Model Parsimony vs. Fidelity
- A minimal model should balance simplicity (for calibration and communication) with enough biological realism to be meaningful and extensible.

---

## Key Deliverables

### 1. Regression Model Prototype
- **Goal:** Connect antibody titers to vaccine efficacy.
- **Scope:** Exclude complexities such as mixed cohorts (e.g., vaccination and natural infection).
- **Deadline:** *Thursday April 17, 2025*.
- **Output:** A documented prototype, even if incomplete, embedded in the current modeling framework.

### 2. Data Aggregation and Visualization
- **Goal:** Collect and visualize immunogenicity and efficacy data from multiple studies.
- **Status:** First tangible deliverable to support model design and communication.

### 3. Calibration harness
- **Goal:** calibration framework that links data to be ingested with future model. Defines likelihood functions and model API interface requirements.
- **Status:** to be assigned.

### 4. Stochastic Cohort Simulation 
- **Goal:** Build a simulation that captures variability in infection and protection at the population level.
- **Timeline:** Follows the completion of the regression model and calibration harness.

### 5. Integration into transmission model
- ** Goal:** integrate immunity model into TyphoidSim member of the [Starsim family](https://starsim.org/). 



