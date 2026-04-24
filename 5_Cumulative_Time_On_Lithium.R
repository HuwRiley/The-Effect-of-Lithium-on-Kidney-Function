# =============================================================================
# Script 5: Cumulative Lithium Exposure and eGFR Slope Analysis
# =============================================================================
# Calculates each patient's cumulative years of lithium exposure up to every
# eGFR measurement date, and tests whether greater cumulative exposure is
# associated with a steeper rate of eGFR decline.
#
# Cumulative exposure is computed using calculate_cumulative_lithium_time()
# from Script 1, based on treatment intervals derived from positive lithium
# tests. Each interval is extended by 365 days beyond the final test to
# account for ongoing treatment after the last recorded level.
#
# A linear mixed effects model then tests the association between cumulative
# lithium years and eGFR slope, adjusting for age and sex.
#
# Requires: Scripts 1, 2, 3, 4 (uses binary_lithium_labels from Script 4)
#
# Outputs:
#   lithium_intervals              : treatment intervals per patient (start/end dates)
#   yearly_slope_cumulative_summary : one row per patient per year, with
#                                     cumulative_lithium_years added
#   cumulative_lme_model_summary    : summary of the cumulative exposure model
# =============================================================================


# =============================================================================
# Step 1: Define lithium treatment intervals per patient
# =============================================================================
# Positive lithium tests are grouped into episodes using the same 365-day
# gap rule as Script 3. Each episode is then extended by 365 days past the
# final positive test to represent the likely ongoing treatment period.

# Extract all positive lithium tests from the full cohort
all_positive_lithium_tests <- all_tests_with_lithium %>%
  filter(test == "Lithium", as.numeric(test_data) > 0.3) %>%
  arrange(study_id, collection_date)

# Label each test with its episode number (same logic as Script 3)
lithium_test_episodes <- all_positive_lithium_tests %>%
  group_by(study_id) %>%
  arrange(collection_date) %>%
  mutate(
    gap_days    = collection_date - lag(collection_date),
    new_episode = if_else(is.na(gap_days) | gap_days > 365, 1, 0),
    episode_id  = cumsum(new_episode)
  )

# Summarise each episode to a start date and end date (+365 days)
lithium_intervals <- lithium_test_episodes %>%
  group_by(study_id, episode_id) %>%
  summarise(
    start_date = min(collection_date),
    end_date   = max(collection_date) + days(365),  # extend beyond final test
    .groups = "drop"
  )


# =============================================================================
# Step 2: Calculate cumulative lithium exposure at each eGFR measurement date
# =============================================================================
# For each row in binary_lithium_labels (one row per creatinine test),
# calculate_cumulative_lithium_time() sums the days of lithium exposure
# that occurred before that eGFR measurement date, across all intervals.

egfr_with_cumulative_exposure <- binary_lithium_labels

egfr_with_cumulative_exposure$cumulative_lithium_years <- map2_dbl(
  egfr_with_cumulative_exposure$study_id,
  egfr_with_cumulative_exposure$collection_date,
  ~ calculate_cumulative_lithium_time(
    patient_id = .x,
    gfr_date   = .y,
    intervals  = lithium_intervals
  )
)


# =============================================================================
# Step 3: Collapse to one observation per patient per calendar year
# =============================================================================
# Takes the maximum cumulative exposure within each year (i.e. end-of-year
# total) and the median eGFR slope, to avoid patients with more frequent
# testing having disproportionate influence on the model.

yearly_slope_cumulative_summary <- egfr_with_cumulative_exposure %>%
  group_by(study_id, calendar_year) %>%
  summarise(
    lithium_on_year          = any(lithium_on,             na.rm = TRUE),
    egfr_slope_year          = median(egfr_slope_2yr,      na.rm = TRUE),
    cumulative_lithium_years = max(cumulative_lithium_years, na.rm = TRUE),
    age_year                 = median(age_at_test,          na.rm = TRUE),
    sex                      = first(sex),
    n_measurements           = n(),
    .groups = "drop"
  )


# =============================================================================
# Step 4: Fit cumulative exposure LME model
# =============================================================================
# Tests whether each additional year of cumulative lithium exposure is
# associated with a greater rate of eGFR decline (more negative slope).
#
# Fixed effects:
#   cumulative_lithium_years : total years on lithium up to each observation
#   age_year                 : age in years (covariate)
#   sex                      : male/female (covariate)
#
# Random effect: (1 | study_id) = patient-level intercept

cumulative_lme_model <- lmer(
  egfr_slope_year ~ cumulative_lithium_years + age_year + sex + (1 | study_id),
  data = yearly_slope_cumulative_summary
)

cumulative_lme_model_summary <- summary(cumulative_lme_model)
cumulative_lme_model_summary

# Profile 95% confidence intervals
confint(cumulative_lme_model, method = "profile")


# =============================================================================
# End of Script 5
# Outputs: lithium_intervals, yearly_slope_cumulative_summary,
#          cumulative_lme_model_summary
# =============================================================================
