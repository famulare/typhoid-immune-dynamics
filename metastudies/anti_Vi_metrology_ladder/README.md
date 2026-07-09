# Anti-Vi Metrology Bridge — a Claude-assisted study, start to finish

This folder is the complete record of an attempt to build a **metrology bridge**: an
empirical conversion cascade from 1960s–80s Vi hemagglutination (HA) titers to modern
anti-Vi IgG VaccZyme ELISA units, so historical typhoid-immunity data can be used in
contemporary dose-response models. It was built as a deliberate experiment in
**autonomous ("YOLO") scientific work** — Claude Opus 4.5 running a pre-agreed contract
without stopping until blocked — and it is shared here for two reasons at once: the
**workflow** (how a Claude+human team ran it) and the **product** (what it found).

**The finding is quantified impossibility.** A continuous bridge is *not* defensibly
achievable: there is a ~40-year data void (edge **E2**) between the 1983 Barrett ELISA
and modern in-house ELISAs, with no published study spanning it. The run bridged that
void with an admitted-fabricated structural assumption, propagated the uncertainty
honestly, and the result is a conversion whose 95% prediction intervals span ~2 orders
of magnitude — **E2 alone contributes ~70% of the total variance**. The useful output
is not the conversion table; it is the precise map of *why we can't do better and where
to invest to fix it*. Full results table and figure: [05_inference/README.md](05_inference/README.md#L104).

```
        E1 (Barrett 1983, n=77)      E2  no data          E3 (Lee 2020, r=0.991)
Vi HA ─────────────────────────▶ Barrett ELISA ⇢⇢⇢⇢▶ Modern in-house ─────────────▶ VaccZyme
 titer      data-supported          titer   [τ assumption]   ELISA (µg/ml)              (EU/ml)
```
Detailed graph, edge tables, and uncertainty budget: [05_inference/bridge_graph.md](05_inference/bridge_graph.md).

---

## How to read this folder

The folders are numbered in the **order the work actually happened** (see the commit
trail below). Read top-to-bottom for the story; open any folder to zoom in.

| # | Folder | What's in it |
|---|--------|--------------|
| — | [README.md](README.md) | this narrative (the only document written for the share) |
| 00 | [00_contract/](00_contract/) | the pre-agreed [design contract](00_contract/anti_Vi_metrology_ladder_contract.md), the [onboarding one-pager](00_contract/onboarding_one_pager.md) (science framing), and the [planned 8-phase checklist](00_contract/progress_checklist.md) |
| 01 | [01_reference_model/](01_reference_model/) | setup-era conceptual scaffolding: the [reference/ideal bridge spec](01_reference_model/metrology_bridge_specification.md), [assay→node mapping](01_reference_model/assay_mapping.md), and two planning placeholders never filled in ([cross-cutting](01_reference_model/cross_cutting_observations.md), [triage](01_reference_model/paper_triage.md)) |
| 02 | [02_input_studies/](02_input_studies/) | the study corpus assembled for the run (6 PDFs + a data table) |
| 03 | [03_gap_search/](03_gap_search/) | the run's **first act**: [confirming the E2 void](03_gap_search/e2_search_results.md) |
| 04 | [04_extraction/](04_extraction/) | structured data extracted from the papers (+ the [extraction schema](04_extraction/_TEMPLATE.md)) |
| 05 | [05_inference/](05_inference/) | the models, the Monte Carlo cascade, the [decision log](05_inference/decision_log.md), outputs, and the [caveats/stop-criteria synthesis](05_inference/caveats_and_stop_criteria.md) |
| 06 | [06_review/](06_review/) | an adversarial self-[review](06_review/reviewer2_critique.md) (Major Revision) and the [author response](06_review/author_response.md) |
| 07 | [07_reflection/](07_reflection/) | post-run [reflection and the human↔Claude dialogue](07_reflection/reflection.md) |

The R scripts in `05_inference/` read paths relative to the **repo root** — run them from
there (or from the RStudio project), e.g. `Rscript metastudies/anti_Vi_metrology_ladder/05_inference/end_to_end_bridge.R`. They need only `tidyverse`.

---

## The commit trail — what happened, when

The chronology is grounded in git, not reconstructed from memory. Six commits:

| Commit | Time (PST) | Phase | Autonomy | Contents |
|--------|-----------|-------|----------|----------|
| [`0a8cba0`](https://github.com/famulare/typhoid-immune-dynamics/commit/0a8cba0) | 2026-02-04 17:01 | setup | collaborative | contract + all scaffolding (→ 00, 01) |
| [`636e030`](https://github.com/famulare/typhoid-immune-dynamics/commit/636e030) | 2026-02-04 17:02 | setup | collaborative | input study PDFs (→ 02) |
| [`7b140d6`](https://github.com/famulare/typhoid-immune-dynamics/commit/7b140d6) | 2026-02-04 17:32 | **the run** | **autonomous** | gap search + extraction + models + decision log + caveats + outputs — all in one commit (→ 03, 04, 05) |
| [`373d87b`](https://github.com/famulare/typhoid-immune-dynamics/commit/373d87b) | 2026-02-04 17:37 | post-run | autonomous subagent | Reviewer 2 critique + author response (→ 06) |
| [`a5210c0`](https://github.com/famulare/typhoid-immune-dynamics/commit/a5210c0) | 2026-02-05 09:07 | next day | collaborative | reflection (→ 07) |
| [`067fde9`](https://github.com/famulare/typhoid-immune-dynamics/commit/067fde9) | 2026-02-05 09:27 | next day | collaborative | post-reflection dialogue (→ 07) |

The autonomous scientific work — gap search, extraction, and all modeling — took **~30
minutes** (17:02 → 17:32) and landed as a single commit.

**Autonomy gradient, not a human/agent boundary.** Every commit is co-authored by
Claude, and the contract was co-designed over several prior sessions. Nothing here was
done "by hand" in a Claude-free sense; the whole study was Claude-assisted throughout.
What varies is the *degree of autonomy* — steered collaboration in setup (00–02), an
unsupervised run (03–05), a spawned self-review, and next-day joint reflection. **The
session logs from this era are lost**, so the exact steering granularity inside the
collaborative phases is unrecoverable. Where a step reads as "done collaboratively" vs.
"honestly can't tell," this record leaves it at that rather than reconstructing it.

---

## The story, in order

### Collaborative setup — [`0a8cba0`](https://github.com/famulare/typhoid-immune-dynamics/commit/0a8cba0), [`636e030`](https://github.com/famulare/typhoid-immune-dynamics/commit/636e030)

Before the autonomous run, Mike and Claude co-designed the [8-phase contract](00_contract/anti_Vi_metrology_ladder_contract.md) — extraction templates, a reference
bridge model, prior specifications, and, critically, **pre-committed stop criteria** for
abandoning the continuous bridge if the data couldn't support it. The contract's
decision-tagging convention (`[USER-LOCKED]` / `[ASSISTANT-PROPOSED]` / `[OPEN]`) and the
[reference/ideal model](01_reference_model/metrology_bridge_specification.md) became the
scaffolding the run then built against. The study corpus (Barrett 1983 ×2, Lee 2020,
Rijpkema 2018, the NIBSC 16/138 IFU) was assembled in [02_input_studies/](02_input_studies/).

### The autonomous run — [`7b140d6`](https://github.com/famulare/typhoid-immune-dynamics/commit/7b140d6)

This single commit bundles the whole run, so git can't order the steps within it. The
internal order below is from the run's own narrative and the artifacts, **not** from git.

1. **Gap search — the first move.** Rather than start from the papers it had, the run
   first went looking for the *missing* one: any published comparison between a 1980s Vi
   ELISA and a modern platform. Nine queries across PubMed, Scholar, and NIBSC reports
   returned nothing. [The E2 void was systematically confirmed, not assumed](03_gap_search/e2_search_results.md) — and it became the defining constraint on everything after.

2. **Extraction.** The papers it *did* have were parsed into [structured extractions](04_extraction/): Barrett 1983 Figure 2 (77 paired HA/ELISA sera, manually counted
   cell-by-cell to [barrett_1983_fig2_data.csv](04_extraction/barrett_1983_fig2_data.csv)),
   the Lee 2020 in-house↔VaccZyme regression, and the NIBSC commutability statement.

3. **Inference.** Three edges modeled ([e1](05_inference/e1_ha_elisa_bridge.R),
   [e3](05_inference/e3_inhouse_vacczyme_bridge.R), and the full Monte Carlo cascade in
   [end_to_end_bridge.R](05_inference/end_to_end_bridge.R)). E1 and E3 are data-driven;
   E2 is not. To cross E2 the run wrote an **admittedly fabricated calibration anchor** —
   [`titer 160 ≈ 10 µg/ml`](05_inference/end_to_end_bridge.R#L148-L160) — and a
   [structural uncertainty parameter τ = 1.0 log unit](05_inference/decision_log.md#L30)
   chosen by "round it up conservatively," not from data. Every such choice is in the
   [decision log](05_inference/decision_log.md). The [caveats/stop-criteria synthesis](05_inference/caveats_and_stop_criteria.md#L88) is blunt about it: *"E2 Is Made Up."*

### Adversarial self-review — [`373d87b`](https://github.com/famulare/typhoid-immune-dynamics/commit/373d87b)

The run then spawned an Opus subagent as "Reviewer 2 (adversarial but fair)" to critique
its own work. The [review](06_review/reviewer2_critique.md) returned **Major Revision**
with 6 major + 8 minor concerns — the fabricated calibration, the arbitrary τ, and a
population mismatch (E1 carriers vs. E3 vaccinees) it argued might be *fatal*, not just a
bias term. The [author response](06_review/author_response.md) accepted the major
concerns and proposed ~20 revisions. The separate reviewing agent, lacking sunk cost,
caught what the implementer glossed over.

### Reflection & dialogue — [`a5210c0`](https://github.com/famulare/typhoid-immune-dynamics/commit/a5210c0), [`067fde9`](https://github.com/famulare/typhoid-immune-dynamics/commit/067fde9)

The next morning, Claude wrote a [reflection](07_reflection/reflection.md) — including
that [*"the tau parameter is a confession, not a measurement"*](07_reflection/reflection.md#L21) — followed by a captured human↔Claude dialogue in which Mike reframed the
fabricated anchor as a legitimate way to "work out the full shape of the thing" so long
as its provenance is honest. That exchange is the meta-capstone of the whole exercise.

---

## Status — what is finished, and what is not

**The science is preserved exactly as the autonomous run left it. Nothing has been
re-run, corrected, or completed for this share.** Specifically:

- The review cycle ends at **Major Revision + a revision plan**. The [~20 agreed revisions](06_review/author_response.md) — presenting τ as a scenario rather than a default,
  adding a separate calibration-uncertainty term (κ), re-verifying the Barrett Figure 2
  extraction, renaming the misleadingly-named `parameter_estimates.csv`, and others —
  **were never executed.** This is a snapshot of where the run stopped, not a finished
  analysis.
- `05_inference/outputs/parameter_estimates.csv` keeps its original (Reviewer-2-flagged)
  name; it actually holds bridge *predictions*, not parameter estimates.
- The [reference spec](01_reference_model/metrology_bridge_specification.md) still says
  "Sections 5–7 to be filled"; the two [01_reference_model/](01_reference_model/)
  placeholders are empty planning scaffolds. They are kept as authentic artifacts of the
  plan-vs-execution gap, not tidied away.

**A provenance "honestly can't tell":** the run summary in
[05_inference/README.md](05_inference/README.md) lists a `reviews/opus_review_phase1.md`,
but no file by that name appears on any reachable git ref (checked with `git log --all
--follow` / `-S`); only `reviewer2_critique.md` and `author_response.md` were ever
committed. Whether an early opus-review file existed transiently and uncommitted during
the (log-lost) session is unrecoverable. It is flagged here rather than asserted either
way.

---

## What changed in this reorganization

Everything above `05_inference/README.md` and the numbered folders is the original YOLO
material, **moved verbatim** from a former `claude_yolo/` subdirectory (`git mv`, history
preserved). The **only** file authored for this share is this `README.md`. The **only**
edits to moved files are mechanical: the three R scripts' path strings (updated for the
new locations; the `here` dependency was dropped in favor of repo-root-relative paths)
and, outside this folder, the deep-links in the project's blog post.
