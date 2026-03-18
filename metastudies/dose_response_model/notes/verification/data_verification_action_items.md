# Data Verification: Consolidated Action Items

**Date**: 2026-03-18
**Sources**: Three independent verification reports:
1. Opus v1 (extracts only) — `data_verification_opus.md`
2. Opus v2 (PDFs) — `data_verification_opus_v2_from_pdfs.md`
3. GPT-5.4 (extracts + PDFs) — `data_verification_gpt-5.4.md`

**Agreement matrix**: Where all three agents agree, confidence is high. Where they disagree, the issue needs human eyes.

---

## CRITICAL — Must fix before implementation

### 1. Waddington 10⁵ arm: does it exist?
- **Opus v2 (PDF)**: Says 10⁵ arm does NOT exist. Claims Waddington tested only 10³ and 10⁴.
- **GPT-5.4**: Says W-F-5 and W-I-5 match the extraction (5/5 and 4/5). No flag.
- **Opus v1 (extracts)**: Matches extraction. No flag.
- **The extraction**: Reports 10⁵ (n=5) with specific breakdown data (bacteremia only=0, fever only=2, both=3, neither=0).
- **Conflict**: Opus v2 says the PDF has no 10⁵ group. But the extraction has detailed per-criterion breakdown that's hard to fabricate. The Waddington paper describes a Phase 1 dose escalation (n=5 per dose at 10³, 10⁴, 10⁵) halted when 10⁵ hit 100%, then Phase 2 (additional enrollment at 10³ and 10⁴). The extraction and the plan's narrative both describe this. Opus v2 may have misread the PDF.
- **Assigned to**: **Mike** — Please check Waddington 2014 PDF Table 2 (p.1234) for the 10⁵ row. The Phase 1 design clearly intended 10⁵; the question is whether the paper reports it.

### 2. Waddington 10⁴ denominators: n=16 or n=20?
- **Opus v2 (PDF)**: Says n=20 at 10⁴ with y=13 TD (65%). Claims Table 1 shows n=20 per-protocol at each dose.
- **GPT-5.4**: Says n=16, y=10 matches the extraction. No flag.
- **Opus v1 (extracts)**: Matches extraction (n=16). No flag.
- **The issue**: Phase 1 had n=5 at 10⁴. Phase 2 added n=11 at 10⁴. Total = 16. But Opus v2 claims the per-protocol total in Table 1 is 20 at each dose. This is the key question: does Table 2 (attack rates) report Phase 1+2 combined (n=16 at 10⁴) or something else (n=20)?
- **Assigned to**: **Mike** — Please check Waddington 2014 Table 2 column headers and verify the denominators at 10⁴. The extraction's narrative (p.1231-1232) describes Phase 1: n=5 per dose, Phase 2: n=15 additional at 10³, n=11 additional at 10⁴. So 10³ total = 20, 10⁴ total = 16. But Opus v2 reads n=20 at both. One of them is wrong.

### 3. Hornick 10³: 0/14 confirmed, extraction was wrong
- **Opus v2 (PDF)**: Both Hornick 1966 Figure 2 AND Hornick 1970 Table 1 show 0/14. The extraction's 9/14 was a hallucination.
- **GPT-5.4**: Notes the conflict between the two extractions but cannot resolve from extractions alone. Cites Hornick 1967 extraction as also showing 0/14.
- **Opus v1 (extracts)**: Notes conflict between extractions; cannot resolve.
- **Resolution**: **0/14 is correct.** The Hornick 1966 extraction is wrong and should be corrected. Prerequisite 1 is resolved. The three-way sensitivity analysis (0/14 vs 9/14 vs excluded) can be simplified to just 0/14.
- **Assigned to**: **Claude** — Fix the Hornick 1966 extraction and remove Prerequisite 1 from the plan.

---

## MAJOR — Should fix before implementation

### 4. Oxford infection definition inconsistency
- **All three agents agree**: Waddington and Jin use shedding-only; Darton uses "bacteremia OR shedding." The plan mixes them.
- **Magnitude**: At ~10⁴ bicarb, Darton placebo shows 87% (bact OR shed) while Jin control shows 71% (shedding only). That's a 16-point gap at comparable doses — enough to bias N50_inf.
- **Options**: (a) Use shedding-only everywhere (need to extract shedding-only from Darton, or drop Darton infection obs). (b) Use bact-OR-shed everywhere (need to compute union for Jin from the separately reported bacteremia and shedding counts). (c) Model the definition difference explicitly.
- **GPT-5.4 additional note**: Jin controls had 24/31 bacteremia and 22/31 shedding. The union (bact OR shed) would be at least 24/31 = 77% (since all bacteremic subjects presumably also shed or are a superset). Need to check for subjects who shed but were not bacteremic.
- **Assigned to**: **Mike** — Decide which infection definition to use. If shedding-only: need to check whether Darton reports shedding-only rates (the extraction's Table 4 / stool shedding section may have it). If bact-OR-shed: need to compute Jin union counts from the per-subject data.

### 5. Maryland fever definition heterogeneity
- **All three agents agree**: Hornick uses ≥103°F/24-36h; Gilman uses fever + culture (treatment at ≥103°F/1d, ≥101°F/3d, ≥100°F/5d); Levine uses ≥101°F + culture.
- **Impact**: Levine's 101°F threshold is 2°F lower than Hornick's 103°F. This likely captures milder cases and could explain why Levine and Gilman controls show higher attack rates (40-53%) than Hornick (28%) at the same 10⁵ dose. The plan currently treats them as the same outcome with a study-level random effect absorbing the heterogeneity.
- **Options**: (a) Note it and let the random effect absorb it. (b) Exclude Levine from the primary fever likelihood and use only for sensitivity. (c) Model a definition offset for Levine/Gilman relative to Hornick. (d) Use only Hornick for the canonical Maryland fever curve and Gilman/Levine for the H-antibody stratification only.
- **Assigned to**: **Mike** — Decide how to handle. Option (d) seems cleanest: Hornick is the canonical multi-dose curve; Gilman/Levine add immunity heterogeneity information but with a different definition that should be acknowledged.

### 6. Hornick Table 1 n=116 vs Table 5 n=104 at 10⁵ — potential hidden overlap
- **Opus v2 and GPT-5.4 both flag this**: The 12-subject gap suggests Table 1 may include controls from later studies (Gilman, Levine era) that the plan also includes separately.
- **Risk**: If Hornick Table 1's 116 includes some Gilman/Levine controls, then H-F-5 double-counts with Gil-F and Lev-F observations.
- **Options**: (a) Use Table 5's n=104 for Hornick unvaccinated at 10⁵ instead of Table 1's n=116. (b) Keep n=116 but flag the potential overlap. (c) Investigate whether Table 1 includes non-Quailes strains (the 1966 paper's 10⁵ data is also n=116, which predates Gilman/Levine).
- **Assigned to**: **Mike** — Check whether the Hornick 1966 Figure 2 also shows n=116 at 10⁵. If so, the 116 predates Gilman/Levine and there's no overlap. If the 116 is specific to the 1970 review, it may have accumulated later data.

### 7. Gil-F-rest (n~37, y~19) is fabricated by subtraction
- **All three agents agree**: This row is derived from 64 - 14 - 13 = 37, not extracted from any source.
- **GPT-5.4**: Notes this is the only row that can't be verified against any source.
- **Opus v2**: Confirms it's internally consistent but approximate (depends on rounding of the stratified counts).
- **Options**: (a) Drop it and use only Gil-F-Hlo + Gil-F-Hhi (n=27 total, losing 37 subjects). (b) Keep it but label as derived and note the approximation. (c) Use the pooled Gil-F-ctrl (n=64, y=31) instead of the stratified observations (loses the immunity information).
- **Assigned to**: **Claude** — Implement option (b) in the plan: keep Gil-F-rest but mark it as derived, add a sensitivity analysis dropping it.

### 8. Gilman H-antibody subsample representativeness
- **Opus v2 and GPT-5.4 agree**: H-antibody was measured in 27 of 64 controls, not all.
- **Prerequisite 2 resolved**: It's a subsample. But is it representative?
- **Assigned to**: **Mike** — Check Gilman 1977 for any description of how the 27 were selected. If it's the first 27 enrolled, or a random sample, it's likely representative. If it's a convenience sample, the pi_susc estimate from this subsample may be biased.

---

## MODERATE — Note and address during implementation

### 9. Gilman shedding definition differs from other studies
- **GPT-5.4 flags**: Gil-I-ctrl uses "late shedding 4-30 days" (26/43), not "any time" shedding.
- **Opus v2 confirms**: Table 4 in Gilman shows 4-30 day shedding specifically.
- **Impact**: "Late shedding" likely underestimates total infection compared to "any time" shedding (Waddington, Levine). Early shedding (0-3 days) is common but may represent pass-through rather than true colonization.
- **Assigned to**: **Claude** — Note this in the plan's data catalog. The 4-30 day definition is arguably a better proxy for true colonization than any-time shedding, but it's different from Waddington's definition.

### 10. Additional Levine shedding data available
- **Opus v2 identifies**: Levine Table 2 has "Stool Anytime" for all 4 trials, not just Trial 1. Trials 2-4 shedding data: 15/33 (46%), 17/22 (77%), 6/16 (38%).
- **GPT-5.4 also notes this** in the missing data inventory.
- **Assigned to**: **Claude** — Add Lev-I-2, Lev-I-3, Lev-I-4 to the data catalog as available but not yet included. Note that these use "any time" shedding, different from Gilman's "4-30 day" definition.

### 11. Section 11 bookkeeping
- **GPT-5.4 flags**: Section 11 shows 19 Oxford + 15 active Maryland = 34 active rows, not 35. The Hornick vaccinated rows (H-V-K5, H-V-L5, H-V-K7, H-V-L7) are in Section 4 but not in Section 11.
- **Assigned to**: **Claude** — Reconcile Section 4 and Section 11. Either add the vaccine rows to Section 11 or explicitly note they are excluded from the primary likelihood.

### 12. Waddington shedding at 10⁴ needs individual-level verification
- **Opus v2 flags**: The per-individual shedding count at 10⁴ (claimed 10/n) cannot be confirmed from the published tables — the PDF doesn't clearly tabulate individual-level shedding by dose group.
- **Assigned to**: **Mike** — If supplementary data (S1 Dataset referenced in Waddington) is accessible, verify the per-dose shedding counts. Otherwise mark W-I-4 shedding numerator as "from extraction, not directly verified against PDF."

---

## MINOR — Document only

### 13. Jin dose uncertainty
- All agents note: Jin reports "1-5 × 10⁴" as target range with no median actual dose. Plan uses "~1e4" which is the low end. The geometric mean of the range [1e4, 5e4] is ~2.2e4.
- **Assigned to**: **Claude** — Note the dose uncertainty. Consider using ~2e4 rather than ~1e4 as the point estimate, or model dose uncertainty explicitly.

### 14. Darton post-vaccination GMT for M01ZH09 and Ty21a not in extraction
- **GPT-5.4 flags**: The plan describes these qualitatively but the extraction doesn't tabulate post-vaccination GMTs for these two vaccines.
- **Assigned to**: **Claude** — Check Darton 2016 PDF for post-vaccination anti-Vi GMT for M01ZH09 and Ty21a arms. If not available (M01ZH09 had minimal Vi response), note that CoP for these arms must come from alternative markers or be treated as weakly constrained.

### 15. Hornick 1966 extraction needs correction
- The extraction reports 9/14 at 10³ which is confirmed wrong by the PDF.
- **Assigned to**: **Claude** — Correct the Hornick_1966.md extraction file.

---

## Summary counts

| Priority | Count | Assigned to Mike | Assigned to Claude |
|----------|-------|-----------------|-------------------|
| CRITICAL | 3 | 2 (items 1-2) | 1 (item 3) |
| MAJOR | 5 | 3 (items 4-6, 8) | 2 (items 7, 8 partial) |
| MODERATE | 4 | 1 (item 12) | 3 (items 9-11) |
| MINOR | 3 | 0 | 3 (items 13-15) |
| **Total** | **15** | **6** | **9** |

Mike's items require human judgment (definition choices, PDF verification) or access to primary data.
Claude's items are mechanical fixes to the plan, extraction, or data catalog.
