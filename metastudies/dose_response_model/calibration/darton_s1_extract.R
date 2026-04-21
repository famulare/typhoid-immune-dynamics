#' ---
#' title: "Darton 2016 S1 Individual-Level Data Extraction"
#' output:
#'   md_document
#' ---
#'
#' # Darton 2016 S1 Data Extraction
#'
#' Extract individual-level data from the Darton 2016 supplementary dataset
#' to resolve the Oxford infection definition inconsistency (shedding-only vs
#' bacteremia-OR-shedding) and enable individual-level CoP analysis.
#'
#' Source: Darton et al. (2016) PLoS NTD, S1 Dataset.

library(tidyverse)
library(readxl)

#' ## Setup

xlsx_path <- file.path(
  "metastudies/dose_response_model/input_papers",
  "Darton et al._2016_Using a Human Challenge Model of Infection to Measure Vaccine Efficacy A Randomised, Controlled Tri - S1 data.xlsx"
)
output_dir <- "metastudies/dose_response_model/analysis_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(xlsx_path))

#' ## Step 1: Read and clean Endpoints sheet

endpoints_raw <- read_excel(xlsx_path, sheet = "Endpoints")

#' The sheet has metadata in the rightmost columns (column 17-18) and some
#' formula strings in TTT38. Clean it up.

endpoints <- endpoints_raw |>
  select(1:16) |>
  rename(
    group = Group,
    day_td = DayTD,
    itt_dx = `ITT Dx`,
    ppp_dx = `PPP Dx`,
    td_criteria = TDCriteria,
    tt_dx_hrs = TTDxHrs,
    tt_t38 = TTT38,
    tt_pos_bc = TTPosBC,
    t38_then_pos_bc = T38thenPosBC,
    t37 = T37,
    t37_5 = T37.5,
    t38 = T38,
    t38_5 = T38.5,
    t39 = T39,
    pos_bc_or_stool = PosBCorStool,
    comments = Comments
  ) |>
  mutate(
    row_id = row_number(),
    #' Derive bacteremia: TTPosBC > 0 means had a positive blood culture
    #' TTPosBC = 0 means never had positive BC
    bacteremia = as.integer(!is.na(tt_pos_bc) & tt_pos_bc > 0),
    #' Per-protocol: exclude the withdrawn subject
    is_ppp = !is.na(ppp_dx)
  )

cat("Endpoints: ", nrow(endpoints), "total rows\n")
cat("Per-protocol:", sum(endpoints$is_ppp), "subjects\n")
cat("Withdrawn:", sum(!endpoints$is_ppp), "subject(s)\n")
cat("Groups (all):", table(endpoints$group), "\n")

#' Filter to per-protocol population
endpoints_ppp <- endpoints |> filter(is_ppp)

cat("\nPer-protocol by group:\n")
print(table(endpoints_ppp$group))

#' ## Step 2: Compute shedding from Microbiology sheet

micro_raw <- read_excel(xlsx_path, sheet = "Microbiology results")

micro <- micro_raw |>
  select(1:9) |>
  rename(
    subject_id = SubjectID,
    specimen = Specimen,
    culture_result = CultureResult,
    s_typhi = STyphi,
    contaminant = Contaminant,
    group = Group,
    td_day_no = TDDayNo,
    sample_day_no = SampleDayNo,
    td_sample_day_no = TDSampleDayNo
  )

#' Get all unique subjects
all_subjects <- micro |>
  distinct(subject_id, group)
cat("\nMicrobiology: ", nrow(micro), "samples from",
    nrow(all_subjects), "unique subjects\n")

#' Compute per-subject stool shedding (any positive faecal culture)
stool_shedding <- micro |>
  filter(specimen == "FAECES") |>
  group_by(subject_id, group) |>
  summarise(
    n_stool_samples = n(),
    n_stool_positive = sum(s_typhi == 1, na.rm = TRUE),
    stool_positive = as.integer(any(s_typhi == 1, na.rm = TRUE)),
    .groups = "drop"
  )

#' Also compute per-subject blood culture positivity from microbiology
#' (as a cross-check against the Endpoints bacteremia field)
blood_culture <- micro |>
  filter(specimen == "BLOOD") |>
  group_by(subject_id, group) |>
  summarise(
    n_blood_samples = n(),
    n_blood_positive = sum(s_typhi == 1, na.rm = TRUE),
    blood_positive = as.integer(any(s_typhi == 1, na.rm = TRUE)),
    .groups = "drop"
  )

#' Merge stool and blood into per-subject microbiology summary
micro_per_subject <- all_subjects |>
  left_join(stool_shedding |> select(subject_id, stool_positive),
            by = "subject_id") |>
  left_join(blood_culture |> select(subject_id, blood_positive),
            by = "subject_id") |>
  mutate(
    stool_positive = replace_na(stool_positive, 0L),
    blood_positive = replace_na(blood_positive, 0L)
  )

cat("\nShedding by group (from Microbiology):\n")
micro_per_subject |>
  group_by(group) |>
  summarise(
    n = n(),
    n_stool_positive = sum(stool_positive),
    n_blood_positive = sum(blood_positive),
    .groups = "drop"
  ) |>
  print()

#' ## Step 3: Map Subject IDs to Endpoints rows
#'
#' The Endpoints sheet lacks SubjectID. Check if we can map by group order.

cat("\n--- Subject ID Mapping ---\n")
cat("Endpoints rows per group (all, incl withdrawn):\n")
print(table(endpoints$group))
cat("Microbiology unique subjects per group:\n")
print(table(micro_per_subject$group))

#' The Blood tests sheet also has SubjectID + Group — use it to get
#' the definitive subject list with ordering.
blood_tests <- read_excel(xlsx_path, sheet = "Blood tests") |>
  select(1:2) |>
  rename(subject_id = SubjectID, group = Group)

subject_order <- blood_tests |>
  distinct(subject_id, group) |>
  arrange(match(group, c("M01ZH09", "Placebo", "Ty21a")), subject_id) |>
  group_by(group) |>
  mutate(within_group_order = row_number()) |>
  ungroup()

cat("\nSubjects per group from Blood tests:\n")
print(table(subject_order$group))

#' Now check: do we have the same N per group as Endpoints?
#' Endpoints has M01ZH09=33, Placebo=30, Ty21a=30 (including 1 withdrawn M01ZH09)
#' Blood tests should have the same subjects.

#' Strategy: assign subject_ids to Endpoints rows by matching within-group order.
#' Both sheets are ordered by group (M01ZH09 first, then Placebo, then Ty21a).

endpoints_with_id <- endpoints |>
  group_by(group) |>
  mutate(within_group_order = row_number()) |>
  ungroup() |>
  left_join(
    subject_order |> select(subject_id, group, within_group_order),
    by = c("group", "within_group_order")
  )

n_matched <- sum(!is.na(endpoints_with_id$subject_id))
cat("\nSubject ID mapping: ", n_matched, "of", nrow(endpoints), "matched\n")

if (n_matched < nrow(endpoints)) {
  cat("WARNING: Not all endpoints rows matched to subject IDs.\n")
  cat("Unmatched rows:\n")
  endpoints_with_id |> filter(is.na(subject_id)) |> select(row_id, group, td_criteria) |> print()
}

#' ## Step 4: Read ELISA data

elisa_raw <- read_excel(xlsx_path, sheet = "ELISA data")

elisa <- elisa_raw |>
  rename(
    group = Group,
    vi_igg_day_neg28 = `Vi IgG Day-28`,
    vi_igg_day0 = `Vi IgG Day 0`
  ) |>
  select(group, vi_igg_day_neg28, vi_igg_day0) |>
  mutate(
    elisa_row = row_number(),
    #' Flag below-LOD values (3.7 = LOD/2, the LOD is 7.4 EU/mL)
    vi_baseline_detectable = as.integer(vi_igg_day_neg28 > 7.4),
    vi_prechallenge_detectable = as.integer(vi_igg_day0 > 7.4)
  )

cat("\nELISA: ", nrow(elisa), "rows\n")

#' The ELISA sheet is NOT in the same row order as Endpoints — it interleaves
#' groups by enrollment order rather than group-blocked. We must merge by
#' (group, within_group_order) instead of raw row order.
elisa <- elisa |>
  group_by(group) |>
  mutate(within_group_order = row_number()) |>
  ungroup()

#' Verify per-group counts match
elisa_counts <- table(elisa$group)
endpt_counts <- table(endpoints$group)
cat("\nELISA per-group counts:", elisa_counts, "\n")
cat("Endpoints per-group counts:", endpt_counts, "\n")
if (all(elisa_counts == endpt_counts)) {
  cat("Per-group counts MATCH — safe to merge by (group, within_group_order).\n")
} else {
  cat("WARNING: Per-group counts DO NOT MATCH.\n")
}

#' ## Step 5: Merge into individual-level dataset

individual <- endpoints_with_id |>
  #' Add ELISA data by (group, within_group_order) — NOT by raw row order
  left_join(
    elisa |> select(group, within_group_order,
                     vi_igg_day_neg28, vi_igg_day0,
                     vi_baseline_detectable, vi_prechallenge_detectable),
    by = c("group", "within_group_order")
  ) |>
  rename(
    vi_igg_baseline = vi_igg_day_neg28,
    vi_igg_prechallenge = vi_igg_day0
  ) |>
  #' Add microbiology shedding by subject_id
  left_join(
    micro_per_subject |> select(subject_id, stool_positive, blood_positive),
    by = "subject_id"
  )

#' Filter to per-protocol and select output columns
individual_out <- individual |>
  filter(is_ppp) |>
  select(
    subject_id, group, ppp_dx, td_criteria,
    t37, t37_5, t38, t38_5, t39,
    bacteremia, stool_positive, blood_positive, pos_bc_or_stool,
    vi_igg_baseline, vi_igg_prechallenge,
    vi_baseline_detectable, vi_prechallenge_detectable
  ) |>
  rename(
    fever_td = ppp_dx,
    fever_37 = t37,
    fever_37_5 = t37_5,
    fever_38 = t38,
    fever_38_5 = t38_5,
    fever_39 = t39,
    bact_or_stool = pos_bc_or_stool
  )

write_csv(individual_out, file.path(output_dir, "darton_individual_endpoints.csv"))
cat("\nWrote darton_individual_endpoints.csv:", nrow(individual_out), "rows\n")

#' ## Step 6: Group summaries

group_summary <- individual_out |>
  group_by(group) |>
  summarise(
    n_ppp = n(),
    n_td = sum(fever_td, na.rm = TRUE),
    n_fever_37 = sum(fever_37, na.rm = TRUE),
    n_fever_37_5 = sum(fever_37_5, na.rm = TRUE),
    n_fever_38 = sum(fever_38, na.rm = TRUE),
    n_fever_38_5 = sum(fever_38_5, na.rm = TRUE),
    n_fever_39 = sum(fever_39, na.rm = TRUE),
    n_bacteremia = sum(bacteremia, na.rm = TRUE),
    n_stool_positive = sum(stool_positive, na.rm = TRUE),
    n_bact_or_stool = sum(bact_or_stool, na.rm = TRUE),
    gmt_vi_igg_day0 = exp(mean(log(vi_igg_prechallenge), na.rm = TRUE)),
    frac_vi_detectable_baseline = mean(vi_baseline_detectable, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(group_summary, file.path(output_dir, "darton_group_summaries.csv"))
cat("\nWrote darton_group_summaries.csv\n")
print(group_summary)

#' ## Step 7: Cross-tabulations

cross_tab <- individual_out |>
  group_by(group) |>
  summarise(
    stool_pos_fever_pos = sum(stool_positive == 1 & fever_td == 1, na.rm = TRUE),
    stool_pos_fever_neg = sum(stool_positive == 1 & fever_td == 0, na.rm = TRUE),
    stool_neg_fever_pos = sum(stool_positive == 0 & fever_td == 1, na.rm = TRUE),
    stool_neg_fever_neg = sum(stool_positive == 0 & fever_td == 0, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(cross_tab, file.path(output_dir, "darton_cross_tabulation.csv"))
cat("\nWrote darton_cross_tabulation.csv\n")
cat("\nCross-tabulation (stool shedding x TD diagnosis):\n")
print(cross_tab)

#' Also compute cross-tab for stool x fever_38 and stool x fever_39
cross_tab_38 <- individual_out |>
  group_by(group) |>
  summarise(
    stool_pos_f38_pos = sum(stool_positive == 1 & fever_38 == 1, na.rm = TRUE),
    stool_pos_f38_neg = sum(stool_positive == 1 & fever_38 == 0, na.rm = TRUE),
    stool_neg_f38_pos = sum(stool_positive == 0 & fever_38 == 1, na.rm = TRUE),
    stool_neg_f38_neg = sum(stool_positive == 0 & fever_38 == 0, na.rm = TRUE),
    .groups = "drop"
  )
cat("\nCross-tabulation (stool x fever>=38C):\n")
print(cross_tab_38)

cross_tab_39 <- individual_out |>
  group_by(group) |>
  summarise(
    stool_pos_f39_pos = sum(stool_positive == 1 & fever_39 == 1, na.rm = TRUE),
    stool_pos_f39_neg = sum(stool_positive == 1 & fever_39 == 0, na.rm = TRUE),
    stool_neg_f39_pos = sum(stool_positive == 0 & fever_39 == 1, na.rm = TRUE),
    stool_neg_f39_neg = sum(stool_positive == 0 & fever_39 == 0, na.rm = TRUE),
    .groups = "drop"
  )
cat("\nCross-tabulation (stool x fever>=39C):\n")
print(cross_tab_39)

#' ## Step 8: Verification checks
#'
#' Compare computed values to published Darton 2016 Table 2.

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("VERIFICATION CHECKS vs Darton 2016 Table 2\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

expected <- tribble(
  ~group,     ~metric,          ~expected, ~source,
  "Placebo",  "n_ppp",          30,        "Table 2",
  "M01ZH09",  "n_ppp",          31,        "Table 2",
  "Ty21a",    "n_ppp",          30,        "Table 2",
  "Placebo",  "n_td",           20,        "Table 2",
  "M01ZH09",  "n_td",           18,        "Table 2",
  "Ty21a",    "n_td",           13,        "Table 2",
  "Placebo",  "n_bact_or_stool", 26,       "Table 2 sensitivity",
  "M01ZH09",  "n_bact_or_stool", 21,       "Table 2 sensitivity",
  "Ty21a",    "n_bact_or_stool", 16,       "Table 2 sensitivity",
  "Placebo",  "n_fever_38",     18,        "Table 2 sensitivity",
  "M01ZH09",  "n_fever_38",     16,        "Table 2 sensitivity",
  "Ty21a",    "n_fever_38",     9,         "Table 2 sensitivity",
  "Placebo",  "n_fever_39",     9,         "Table 2 sensitivity",
  "M01ZH09",  "n_fever_39",     9,         "Table 2 sensitivity",
  "Ty21a",    "n_fever_39",     3,         "Table 2 sensitivity"
)

gs_long <- group_summary |>
  pivot_longer(-group, names_to = "metric", values_to = "computed")

checks <- expected |>
  left_join(gs_long, by = c("group", "metric")) |>
  mutate(
    match = computed == expected,
    status = if_else(match, "PASS", "FAIL")
  )

print(checks |> select(group, metric, expected, computed, status))

n_pass <- sum(checks$match)
n_fail <- sum(!checks$match)
cat("\n", n_pass, "PASS,", n_fail, "FAIL out of", nrow(checks), "checks\n")

if (n_fail > 0) {
  cat("\nFAILED CHECKS:\n")
  checks |> filter(!match) |> print()
}

#' ### Key new value: Darton placebo SHEDDING-ONLY rate
cat("\n--- KEY NEW RESULT ---\n")
cat("Darton placebo shedding-only:",
    group_summary |> filter(group == "Placebo") |> pull(n_stool_positive),
    "/",
    group_summary |> filter(group == "Placebo") |> pull(n_ppp),
    "=",
    round(group_summary |> filter(group == "Placebo") |> pull(n_stool_positive) /
            group_summary |> filter(group == "Placebo") |> pull(n_ppp) * 100, 1),
    "%\n")
cat("(vs published bact-OR-stool: 26/30 = 86.7%)\n")
cat("(vs Jin control shedding-only: 22/31 = 71.0%)\n")

#' ### Baseline anti-Vi check
cat("\n--- BASELINE ANTI-Vi CHECK ---\n")
for (g in c("M01ZH09", "Placebo", "Ty21a")) {
  frac <- group_summary |> filter(group == g) |> pull(frac_vi_detectable_baseline)
  cat(g, ": ", round(frac * 100, 0), "% with baseline Vi > 7.4 EU/mL\n")
}
cat("Expected from Darton p.13: M01ZH09 19%, Placebo 40%, Ty21a 28%\n")

#' ## Step 9: Summary plots

#' ### Plot 1: Attack rates by group and endpoint
plot_data <- group_summary |>
  mutate(across(starts_with("n_"), ~ . / n_ppp, .names = "rate_{.col}")) |>
  select(group, starts_with("rate_")) |>
  pivot_longer(-group, names_to = "endpoint", values_to = "rate") |>
  mutate(
    endpoint = str_remove(endpoint, "rate_n_") |>
      fct_recode(
        "TD composite" = "td",
        "Fever >=37" = "fever_37",
        "Fever >=37.5" = "fever_37_5",
        "Fever >=38" = "fever_38",
        "Fever >=38.5" = "fever_38_5",
        "Fever >=39" = "fever_39",
        "Bacteremia" = "bacteremia",
        "Stool shedding" = "stool_positive",
        "Bact OR Stool" = "bact_or_stool"
      ) |>
      fct_relevel("Stool shedding", "Bacteremia", "Bact OR Stool",
                  "TD composite", "Fever >=37", "Fever >=37.5",
                  "Fever >=38", "Fever >=38.5", "Fever >=39"),
    group = fct_relevel(group, "Placebo", "M01ZH09", "Ty21a")
  )

p1 <- ggplot(plot_data, aes(x = endpoint, y = rate, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Darton 2016: Attack rates by endpoint and vaccine group",
    subtitle = "Individual-level data from S1 Dataset (per-protocol population)",
    x = NULL, y = "Attack rate", fill = "Group"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#' ### Plot 2: Individual anti-Vi IgG vs outcome
p2 <- individual_out |>
  mutate(
    group = fct_relevel(group, "Placebo", "M01ZH09", "Ty21a"),
    td_outcome = if_else(fever_td == 1, "Diagnosed", "Not diagnosed")
  ) |>
  ggplot(aes(x = group, y = vi_igg_prechallenge, color = td_outcome)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.7, size = 2.5) +
  scale_y_log10(labels = scales::comma) +
  geom_hline(yintercept = 7.4, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.5, y = 7.4, label = "LOD (7.4)", hjust = 0, vjust = -0.5,
           color = "grey50", size = 3) +
  labs(
    title = "Darton 2016: Pre-challenge anti-Vi IgG by outcome",
    subtitle = "Individual-level data; dashed line = limit of detection",
    x = NULL, y = "Anti-Vi IgG at challenge (EU/mL, log scale)",
    color = "TD Outcome"
  ) +
  theme_minimal(base_size = 11)

#' ### Plot 3: Cross-tabulation
cross_tab_plot <- cross_tab |>
  pivot_longer(-group, names_to = "cell", values_to = "count") |>
  mutate(
    stool = if_else(str_detect(cell, "stool_pos"), "Stool +", "Stool -"),
    fever = if_else(str_detect(cell, "fever_pos"), "TD +", "TD -"),
    group = fct_relevel(group, "Placebo", "M01ZH09", "Ty21a")
  )

p3 <- ggplot(cross_tab_plot, aes(x = stool, y = fever, fill = count)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = count), size = 5, fontface = "bold") +
  facet_wrap(~group) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Darton 2016: Cross-tabulation of shedding x TD diagnosis",
    subtitle = "Per-protocol population, individual-level",
    x = "Stool shedding", y = "TD diagnosis", fill = "Count"
  ) +
  theme_minimal(base_size = 11)

#' Save plots
pdf(file.path(output_dir, "darton_s1_plots.pdf"), width = 10, height = 7)
print(p1)
print(p2)
print(p3)
dev.off()

#' Also save as individual PNGs for easy viewing
ggsave(file.path(output_dir, "darton_s1_attack_rates.png"), p1,
       width = 10, height = 6, dpi = 150)
ggsave(file.path(output_dir, "darton_s1_anti_vi_vs_outcome.png"), p2,
       width = 8, height = 6, dpi = 150)
ggsave(file.path(output_dir, "darton_s1_cross_tabulation.png"), p3,
       width = 10, height = 5, dpi = 150)

cat("\nPlots saved to", output_dir, "\n")
cat("\nDone!\n")
