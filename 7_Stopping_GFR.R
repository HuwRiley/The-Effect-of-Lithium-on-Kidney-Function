# =============================================================================
# Script 7: eGFR at Time of Lithium Discontinuation
# =============================================================================
# Examines the distribution of kidney function at the point when patients
# stopped lithium, to explore whether clinical decision-making around
# discontinuation is influenced by CKD stage thresholds.
#
# For each patient classified as "stopping lithium" (i.e. their last observed
# lithium status before a gap of >12 months), the most recent creatinine
# result at or before discontinuation is extracted as a proxy for kidney
# function at the time of stopping.
#
# Outputs produced (reported in results):
#   - Summary statistics for eGFR at discontinuation
#   - Distribution of CKD stage at stopping (flextable; Figure 5 in dissertation)
#
# Requires: Script 4 (gfr_slopes with lithium_status column)
#
# Outputs:
#   stopping_egfr_all      : all creatinine rows labelled "stopping lithium"
#   stopping_egfr_patient  : one row per patient (last eGFR before stopping)
#   egfr_at_stopping_stats : summary statistics for eGFR at discontinuation
# =============================================================================


# =============================================================================
# Step 1: Extract all creatinine measurements classified as "stopping lithium"
# =============================================================================
# These are measurements where the patient had a positive lithium test in the
# prior 12 months but no positive test in the following 12 months, within an
# observable window (see Script 4 classification rules).

stopping_egfr_all <- gfr_slopes %>%
  filter(lithium_status == "stopping lithium", !is.na(egfr)) %>%
  arrange(study_id, collection_date)


# =============================================================================
# Step 2: Take the last eGFR measurement per patient before stopping
# =============================================================================
# slice_tail(n = 1) gives the final creatinine result before discontinuation,
# used as the best available measure of kidney function at the time of stopping.

stopping_egfr_patient <- stopping_egfr_all %>%
  group_by(study_id) %>%
  slice_tail(n = 1) %>%
  ungroup()


# =============================================================================
# Step 3: Summary statistics for eGFR at discontinuation
# =============================================================================

egfr_at_stopping_stats <- stopping_egfr_patient %>%
  summarise(
    n          = n(),
    mean_egfr  = mean(egfr,   na.rm = TRUE),
    median_egfr = median(egfr, na.rm = TRUE),
    sd_egfr    = sd(egfr,     na.rm = TRUE),
    min_egfr   = min(egfr,    na.rm = TRUE),
    max_egfr   = max(egfr,    na.rm = TRUE)
  )

egfr_at_stopping_stats


# =============================================================================
# Step 4: Classify each patient into a CKD stage at the time of stopping
# =============================================================================
# CKD stages are based on KDIGO eGFR thresholds:
#   Normal (no CKD) : eGFR >= 90
#   Mild CKD        : eGFR 60–89  (stage G2)
#   Moderate CKD    : eGFR 30–59  (stages G3a/G3b)
#   Severe CKD      : eGFR < 30   (stages G4/G5)

stopping_egfr_patient <- stopping_egfr_patient %>%
  mutate(
    ckd_stage = case_when(
      egfr >= 90 ~ "Normal(>90)",
      egfr >= 60 ~ "Mild(60-90)",
      egfr >= 30 ~ "Moderate(30-60)",
      egfr <  30 ~ "Severe(<30)"
    )
  )


# =============================================================================
# Step 5: CKD stage distribution table (Figure 5)
# =============================================================================

ckd_stage_counts <- stopping_egfr_patient %>%
  count(ckd_stage) %>%
  mutate(proportion = n / sum(n))

ckd_stage_counts %>%
  mutate(
    ckd_stage   = factor(ckd_stage,
                         levels = c("Normal(>90)", "Mild(60-90)",
                                    "Moderate(30-60)", "Severe(<30)")),
    proportion  = paste0(round(proportion * 100, 1), "%")
  ) %>%
  arrange(ckd_stage) %>%
  rename(
    "CKD Stage"  = ckd_stage,
    "N"          = n,
    "Proportion" = proportion
  ) %>%
  flextable() %>%
  set_caption("CKD Stage Distribution at Lithium Discontinuation") %>%
  theme_vanilla() %>%
  autofit()


# =============================================================================
# End of Script 7
# Outputs: stopping_egfr_all, stopping_egfr_patient, egfr_at_stopping_stats
# =============================================================================
