# Typhoid immune dynamics

Research repository for typhoid fever immune dynamics modeling. The project develops mechanistic models connecting bacterial dose, pre-existing immunity, and clinical outcomes—calibrated to human challenge study data and validated against field epidemiology.

Docs and blog at <https://famulare.github.io/typhoid-immune-dynamics/>.

## Project Focus

- **Dose-response modeling**: Modified beta-Poisson framework with immunity scaling
- **Correlates of protection**: Anti-Vi IgG as mechanistic predictor of susceptibility
- **Cohort dynamics**: Age-structured models of exposure, infection, and immunity acquisition
- **Vaccine efficacy**: Simulation frameworks for TCV trial design and interpretation

## Folders

- **data/**: Raw data files (seroprevalence, challenge study results)

- **scratch/**: Exploratory R scripts and prototypes

- **metastudies/**: Systematic data extraction and model calibration
    - **metastudies/dose_response_model/**: Challenge study meta-analysis (see folder README)

- **analysis/**: Structured analysis outputs

- **docs/**: Source material for the documentation website
    - **docs/blog/**: Blog posts (rendered from R scripts via knitr::spin)
    - **docs/docs/**: Stand-alone documentation pages