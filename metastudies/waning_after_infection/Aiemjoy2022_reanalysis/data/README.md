# Data Provenance

## Source
Downloaded from: https://github.com/kaiemjoy/TyphoidSeroIncidence
File: `data/TypoidCaseData_github_09.30.21.csv`
Download date: 2026-03-24
Repository has no LICENSE file as of download date.

## Citation
Aiemjoy K, Seidman JC, Saha S, et al. Estimating typhoid incidence from community-based serosurveys: a multicohort study. Lancet Microbe 2022; 3: e578–87. https://doi.org/10.1016/S2666-5247(22)00114-8

Paper data sharing statement: "All analysis code is available at https://github.com/kaiemjoy/TyphoidSeroIncidence.git"

## Contents
Individual-level longitudinal antibody measurements from 1,666 blood-culture-confirmed enteric fever patients enrolled through the SEAP (Bangladesh, Nepal, Pakistan) and SETA (Ghana) studies, 2016–2021.

## Column Dictionary
| Column | Description |
|--------|-------------|
| index_id | Subject identifier |
| age | Age in years |
| Country | Study country |
| bldculres | Blood culture result (serovar: S. Typhi, S. Paratyphi A) |
| Hospitalized | Hospitalization status |
| reinf_obs | Suspected reinfection observed during follow-up |
| {Antigen}_{Isotype}_visit{1-7} | Antibody level (ELISA units) at visit 1-7 |
| TimeInDays_visit{1-7} | Days since fever onset at visit 1-7 |

Antigen-isotype combinations: HlyE IgA, HlyE IgG, LPS IgA, LPS IgG, MP IgA, MP IgG, Vi IgG (7 total).

Vi IgG measured in Nepal only (columns 49-55).

## Data Handling Notes
From their code (`v9na.data.r`): zero antibody values are replaced with 0.01 before log transformation.
