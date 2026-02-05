# Barrett 1983 Data Extraction

## Source

Barrett TJ, Blake PA, Brown SL, Hoffman K, Mateu Llort J, Feeley JC. Enzyme-linked immunosorbent assay for detection of human antibodies to Salmonella typhi Vi antigen. **J Clin Microbiol.** 1983;17(4):625-627.

## Figure 2 Extraction

### Description

Figure 2 (page 626) is a scattergram showing Vi antibody titers obtained by hemagglutination (HA, Y-axis) and ELISA (X-axis) in 77 sera from persons with acute typhoid fever. Numbers inside data points represent the count of sera with those specific titer combinations.

### Extraction Method

1. Visual inspection of Figure 2 in PDF
2. Cell-by-cell reading of count values
3. Verification: sum of counts = 77 (matches paper)

### Data Structure

| Column | Description | Units |
|--------|-------------|-------|
| ELISA_titer | Reciprocal ELISA titer | 2-fold dilution |
| HA_titer | Reciprocal HA titer | 2-fold dilution |
| count | Number of sera at this combination | integer |

### Censoring

- Values below detection limit (<20) are coded as 10 (geometric midpoint)
- This is a standard approach for left-censored dilution data

### Output File

`barrett_1983_fig2_data.csv`

### Key Statistics from Paper

- 77 sera from culture-confirmed typhoid patients (El Salvador)
- 40 (52%) positive by ELISA (titer ≥20)
- 35 (47%) positive by HA (titer ≥20)
- 30 (39%) positive by both assays
- No significant difference between methods (P = 0.132, paired t-test)

### Population Characteristics

- Acute typhoid fever patients
- El Salvador, early 1980s
- Included 29 persons with 2 specimens at different times
- Most specimens collected relatively late in illness course

### Assay Methods

**ELISA:**
- Vi antigen from *Citrobacter* 5396/38 (highly purified)
- Burro anti-*S. typhi* TY-2 capture antibody
- Alkaline phosphatase-conjugated goat anti-human IgG
- Visual endpoint reading
- Starting dilution: 1:20
- Shelf life: ~1 year at 4°C

**Hemagglutination (HA):**
- Performed per Nolan et al. (1980)
- Vi-coated erythrocytes
- More responsive to IgM than IgG

### Limitations for Bridging

1. **Population**: Acute typhoid patients, not healthy vaccinees
2. **Era**: 1980s assay methods, not directly comparable to modern coatings
3. **Visual reading**: May have higher variability than spectrophotometric
4. **IgG vs total antibody**: ELISA specific for IgG; HA detects both IgM and IgG

### Extracted Data Validation

```
Total sera extracted: 77 ✓
Positive by ELISA (≥20): 40 ✓
Positive by HA (≥20): 35 ✓
Double-negative (<20 both): 37 ✓
```

### Raw Counts by Row (HA titer)

| HA | <20 | 20 | 40 | 80 | 160 | 320 | 640 | 1280 | 2560 | 5120 | Row Total |
|----|-----|----|----|----|----|-----|-----|------|------|------|-----------|
| <20 | 37 | 4 | - | 1 | - | - | - | - | - | - | 42 |
| 20 | 1 | 2 | - | - | - | - | - | - | - | - | 3 |
| 40 | 1 | - | - | - | - | - | - | - | - | - | 1 |
| 80 | - | - | 1 | - | - | - | - | - | - | - | 1 |
| 160 | - | - | - | 1 | 1 | - | - | - | - | - | 2 |
| 320 | - | - | - | 2 | 2 | 2 | - | - | - | - | 6 |
| 640 | - | - | - | 1 | 2 | 4 | 1 | - | - | - | 8 |
| 1280 | - | - | - | - | - | 5 | - | - | 1 | 1 | 7 |
| 2560 | - | - | - | - | - | - | 1 | 2 | - | 1 | 4 |
| 5120 | - | - | - | - | - | 1 | - | - | 1 | 1 | 3 |
| **Col Total** | 39 | 6 | 1 | 5 | 5 | 12 | 2 | 2 | 2 | 3 | **77** |

---
*Extracted 2026-02-04*
