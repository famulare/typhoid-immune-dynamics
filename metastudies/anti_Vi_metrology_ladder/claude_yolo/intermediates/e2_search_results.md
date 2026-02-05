# E2 Bridge Data Search Results

**Date**: 2026-02-04
**Objective**: Find any bridging studies between Barrett 1983 ELISA and modern in-house ELISAs

## Search Strategy

### Queries Executed

1. `Vi antibody ELISA comparison validation historical bridging study typhoid`
2. `"Vi polysaccharide" ELISA calibration reference standard comparison Barrett`
3. `WHO NIBSC Vi antibody assay comparison technical report collaborative study`
4. `Rijpkema 2018 Vi polysaccharide international standard collaborative study NIBSC`
5. `"Vi-IgGR1,2011" OR "CBER reference" Vi antibody calibration typhoid`
6. `"hemagglutination" "ELISA" Vi typhoid antibody correlation comparison 1980s 1990s`
7. `Barrett Grubb Gruschkau Vi antibody ELISA typhoid carrier hemagglutination`
8. `"Journal of Clinical Microbiology" 1983 Vi antibody ELISA typhoid Barrett`
9. `Vi antibody assay evolution history standardization 1983 2000 typhoid serology`

### Sources Searched
- PubMed/PMC
- Google Scholar
- NIBSC website
- WHO IRIS database (indirect)

## Key Findings

### Confirming E2 Gap

**No direct bridging study was found** between early 1980s Vi ELISAs (Barrett, Grubb, etc.) and modern in-house ELISAs (2000s-present).

The gap is real and spans approximately 40 years of assay evolution.

### Relevant Papers Found

#### Modern Standardization (Post-2010)

1. **Rijpkema 2018** (PMC6238147)
   - Established NIBSC 16/138 as first International Standard
   - Compared VaccZyme commercial ELISA with in-house ELISAs
   - No mention of 1980s assays
   - Commutability confirmed between modern assays

2. **Szu 2013** (PMC3839630)
   - Established Vi-IgGR1,2011 (CBER reference, 33 μg/ml)
   - Weight-based calibration using precipitin analysis
   - Conversion: 1 ELISA Unit ≈ 1.24 μg/ml
   - No comparison to historical assays

3. **Lee 2020** (PMC7156108)
   - In-house ELISA vs VaccZyme comparison
   - r = 0.991, excellent correlation
   - Provides E3 bridge data

#### Historical Methods (1980s)

4. **Barrett 1983** (JCM, Parts 1 & 2)
   - Original Vi ELISA development
   - Compared to hemagglutination
   - Provides E1 bridge data
   - Used glutaraldehyde-fixed Vi coating

5. **Grubb et al. 1980s**
   - Modified Vi hemagglutination tests
   - CIE vs HA comparisons
   - No ELISA bridging

### Why E2 Doesn't Exist

Historical context suggests reasons for the gap:

1. **Lack of standardization focus**: Before NIBSC 16/138, there was no international standard to anchor comparisons
2. **Method evolution**: Vi coating methods changed (glutaraldehyde → poly-L-lysine → biotinylated)
3. **Declining carrier surveillance**: Reduced need for carrier screening assays in developed countries
4. **Vaccine development focus**: Modern assays developed for TCV trials, not carrier detection

### Potential Mitigation

Both old and new assays can theoretically report in μg/ml when calibrated against:
- Precipitin analysis with purified Vi (historical)
- Vi-IgGR1,2011 (modern)

This provides a **common unit** that could bridge the gap, with caveats:
- Unknown systematic bias between assay generations
- Different coating antigens may have different epitope presentations
- Antibody populations differ (carrier vs vaccinee)

## Conclusion

**E2 is confirmed missing.** No direct comparison study exists.

**Recommended approach**: Treat E2 as a structural assumption with explicit uncertainty (τ parameter) rather than aborting the bridge entirely.

## References

- Barrett TJ et al. (1983). J Clin Microbiol 18:625-632.
- Lee C et al. (2020). PLoS Negl Trop Dis 14(4):e0008171.
- Rijpkema S et al. (2018). Biologicals 56:29-36.
- Szu SC et al. (2013). Vaccine 31(18):2238-2242.
