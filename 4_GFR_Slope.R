# =============================================================================
# Script 4: Lithium Status Classification and eGFR Slope Analysis
# =============================================================================
# Assigns a lithium exposure status to every creatinine (eGFR) measurement,
# calculates rolling 2-year eGFR slopes, and fits a linear mixed effects
# model to assess the association between lithium exposure and eGFR slope.
#
# Lithium status classification:
#   Each eGFR measurement is assigned one of four categories based on whether
#   a positive lithium test (>0.3 mmol/L) was recorded within 12 months
#   before and/or after the eGFR date:
#
#     "on lithium"       : positive test within 12m before AND after
#     "starting lithium" : positive test within 12m after but not before
#                          (only assigned if the 12m before is observable)
#     "stopping lithium" : positive test within 12m before but not after
#                          (only assigned if the 12m after is observable)
#     "not on lithium"   : no positive test within 12m either side
#
#   Edge cases: measurements near the start/end of follow-up, or near death,
#   are conservatively labelled "on lithium" to avoid misclassification.
#
# eGFR slope:
#   A rolling 2-year linear regression slope is computed at each eGFR
#   measurement using calculate_local_gfr_slope() from Script 1.
#   Slopes beyond ±10 ml/min/1.73m2/year are excluded as likely reflecting
#   acute illness rather than chronic decline.
#   One median slope per patient per calendar year is used in models to
#   prevent patients with more frequent tests having disproportionate weight.
#
# Requires: Scripts 1, 2, 3
#
# Outputs:
#   all_tests_tagged          : all_tests_with_lithium with lithium_status added
#   gfr_slopes                : creatinine rows with rolling slope calculated
#   binary_lithium_labels     : gfr_slopes with binary on/off lithium flag
#   yearly_slope_summary      : one row per patient per year (used in model)
#   binary_lme_model_summary  : summary of the binary exposure LME model
#   binary_lme_model_ci       : 95% confidence intervals for model coefficients
# =============================================================================


# =============================================================================
# Step 1: Separate creatinine and lithium tests for date comparison
# =============================================================================

# All creatinine (eGFR) measurements from the full lithium-exposed cohort
creatinine_tests <- all_tests_with_lithium %>%
  filter(test == "Creatinine")

# All positive lithium tests (> 0.3 mmol/L threshold)
positive_lithium_tests <- all_tests_with_lithium %>%
  filter(test == "Lithium", test_data > 0.3) %>%
  select(study_id, collection_date)

# Overall date range of the dataset — used to determine whether the ±12m
# window around each eGFR measurement is observable at the start/end of
# follow-up
follow_up_end_date   <- max(all_tests$collection_date, na.rm = TRUE)
follow_up_start_date <- min(all_tests$collection_date, na.rm = TRUE)


# =============================================================================
# Step 2: Cross-reference every eGFR date against every lithium test date
# =============================================================================
# For each eGFR measurement, flag whether any positive lithium test occurred:
#   - within 12 months BEFORE the eGFR date (lithium_before)
#   - on the SAME day as the eGFR date (lithium_same)
#   - within 12 months AFTER the eGFR date (lithium_after)
# Also flag whether a 12-month window backwards/forwards is fully observable
# within the overall follow-up period.

egfr_lithium_crossref <- creatinine_tests %>%
  left_join(positive_lithium_tests, by = "study_id",
            suffix = c("_egfr", "_li")) %>%
  mutate(
    diff_days         = as.numeric(collection_date_li - collection_date_egfr),
    lithium_before    = diff_days >= -365 & diff_days < 0,
    lithium_same      = diff_days == 0,
    lithium_after     = diff_days > 0 & diff_days <= 365,
    # Can we observe 12m before/after this eGFR date within the dataset?
    observable_after  = collection_date_egfr + 365 <= follow_up_end_date,
    observable_before = collection_date_egfr - 365 >= follow_up_start_date
  )

# Last positive lithium test date per patient (used in edge-case classification)
last_lithium_date_per_patient <- positive_lithium_tests %>%
  group_by(study_id) %>%
  summarise(last_lithium_date = max(collection_date), .groups = "drop")


# =============================================================================
# Step 3: Summarise flags per eGFR measurement
# =============================================================================

egfr_flags_summary <- egfr_lithium_crossref %>%
  group_by(study_id, collection_date_egfr) %>%
  summarise(
    lithium_before    = any(lithium_before,    na.rm = TRUE),
    lithium_same      = any(lithium_same,      na.rm = TRUE),
    lithium_after     = any(lithium_after,     na.rm = TRUE),
    observable_after  = first(observable_after),
    observable_before = first(observable_before),
    dod               = first(dod),   # date of death, if recorded
    .groups = "drop"
  ) %>%
  left_join(last_lithium_date_per_patient, by = "study_id")


# =============================================================================
# Step 4: Assign lithium status to each eGFR measurement
# =============================================================================

egfr_lithium_status <- egfr_flags_summary %>%
  mutate(
    # Days since the patient's last ever positive lithium test
    days_since_last_lithium = as.numeric(collection_date_egfr - last_lithium_date),
    # Flag if patient died within 12m after this eGFR measurement
    # (treated as still on lithium to avoid misclassifying as stopping)
    died_within_12m = !is.na(dod) & dod <= collection_date_egfr + 365,

    lithium_status = case_when(

      # --- Standard classifications ---

      # Clearly on lithium: a positive test within 12m on both sides
      lithium_before & lithium_after
        ~ "on lithium",

      # Starting lithium: a test coming up but not one in the past 12m,
      # and we have enough history to confirm no prior test
      (lithium_same | lithium_after) & !lithium_before & observable_before
        ~ "starting lithium",

      # Stopping lithium: a test in the past 12m but none coming up,
      # and we can observe the following 12m (patient did not die soon after)
      (lithium_before | lithium_same) & !lithium_after & observable_after & !died_within_12m
        ~ "stopping lithium",

      # --- Edge cases: insufficient observation window ---

      # Cannot observe 12m backwards — conservatively treat as ongoing
      (lithium_same | lithium_after) & !observable_before
        ~ "on lithium",

      # Cannot observe 12m forward AND last test was recent — still on lithium
      (lithium_before | lithium_same) & !observable_after & days_since_last_lithium <= 365
        ~ "on lithium",

      # Cannot observe 12m forward AND last test was >1 year ago — likely stopped
      (lithium_before | lithium_same) & !observable_after & days_since_last_lithium > 365
        ~ "not on lithium",

      # Died within 12m of eGFR — classify as on lithium rather than stopping
      (lithium_before | lithium_same) & !lithium_after & died_within_12m
        ~ "on lithium",

      # Default: no positive lithium test within 12m on either side
      TRUE ~ "not on lithium"
    )
  ) %>%
  select(study_id, collection_date_egfr, lithium_status)

# Join lithium status back onto the full test dataset
all_tests_tagged <- all_tests_with_lithium %>%
  left_join(
    egfr_lithium_status,
    by = c("study_id", "collection_date" = "collection_date_egfr")
  )


# =============================================================================
# Step 5: Calculate rolling 2-year eGFR slope at each creatinine measurement
# =============================================================================
# Uses calculate_local_gfr_slope() from Script 1. Slopes are only calculated
# where at least 3 creatinine measurements exist in the prior 2 years.
# Slopes exceeding ±10 ml/min/1.73m2/year are set to NA as they likely
# reflect acute kidney injury or data instability rather than chronic decline.

gfr_slopes <- all_tests_tagged %>%
  filter(test == "Creatinine") %>%
  arrange(study_id, collection_date) %>%
  group_by(study_id) %>%
  mutate(
    egfr_slope_2yr = map_dbl(
      collection_date,
      ~ calculate_local_gfr_slope(
        dates       = collection_date,
        gfr         = egfr,
        centre_date = .x
      )
    )
  ) %>%
  ungroup() %>%
  # Remove physiologically implausible slopes
  mutate(
    egfr_slope_2yr = if_else(abs(egfr_slope_2yr) > 10, NA_real_, egfr_slope_2yr)
  )


# =============================================================================
# Step 6: Prepare data for the binary lithium exposure model
# =============================================================================
# For the LME model, lithium status is simplified to a binary variable:
#   TRUE  = on lithium (includes "on lithium" and "stopping lithium"
#           as patients in the latter group were likely on lithium for most
#           of the 2-year slope window)
#   FALSE = off lithium (includes "not on lithium" and "starting lithium"
#           as patients in the latter group were likely off lithium for most
#           of the 2-year slope window)

binary_lithium_labels <- gfr_slopes %>%
  mutate(
    calendar_year = year(collection_date),
    lithium_on    = lithium_status %in% c("on lithium", "stopping lithium")
  )

# Collapse to one median slope per patient per calendar year.
# This prevents patients with more frequent monitoring from having
# disproportionate influence on the model.
yearly_slope_summary <- binary_lithium_labels %>%
  group_by(study_id, calendar_year) %>%
  summarise(
    lithium_on_year  = any(lithium_on,        na.rm = TRUE),
    egfr_slope_year  = median(egfr_slope_2yr, na.rm = TRUE),
    age_year         = median(age_at_test,     na.rm = TRUE),
    sex              = first(sex),
    n_measurements   = n(),
    .groups = "drop"
  )


# =============================================================================
# Step 7: Fit binary lithium exposure LME model
# =============================================================================
# Linear mixed effects model assessing the association between binary lithium
# exposure (on vs off) and eGFR slope, adjusting for age and sex.
# A random intercept per patient accounts for repeated measurements.
#
# Fixed effects:
#   lithium_on_year : TRUE = on lithium (reference: FALSE = not on lithium)
#   age_year        : age in years (covariate)
#   sex             : male/female (covariate)
#
# Random effect: (1 | study_id) = patient-level intercept

binary_lme_model <- lmer(
  egfr_slope_year ~ lithium_on_year + age_year + sex + (1 | study_id),
  data = yearly_slope_summary
)

binary_lme_model_summary <- summary(binary_lme_model)
binary_lme_model_summary

# Profile 95% confidence intervals for all fixed effect coefficients
binary_lme_model_ci <- confint(binary_lme_model, method = "profile")
binary_lme_model_ci


# =============================================================================
# End of Script 4
# Outputs: all_tests_tagged, gfr_slopes, binary_lithium_labels,
#          yearly_slope_summary, binary_lme_model_summary, binary_lme_model_ci
# =============================================================================
