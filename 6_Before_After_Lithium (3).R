# =============================================================================
# Script 6: Before / During / After Lithium Phase Analysis
# =============================================================================
# Compares eGFR slope across three phases of lithium treatment:
#   "before"  : creatinine tests before the patient's first positive lithium test
#   "during"  : creatinine tests during the lithium treatment episode
#   "after"   : creatinine tests after the patient's final positive lithium test
#
# This analysis is restricted to patients with a single uninterrupted lithium
# episode (from Script 3) and a minimum episode length of 5 years, to ensure
# all three phases are meaningfully represented.
#
# Rolling 2-year eGFR slopes are computed as in Script 4. Slopes are then
# collapsed to one median per patient per year within each phase, and a linear
# mixed effects model estimates the mean slope in each phase.
#
# Requires: Scripts 1, 2, 3 (uses all_tests_exclusions with episode dates)
#
# Outputs:
#   phase_data                : creatinine tests with lithium phase labels
#   episode_length_stats      : summary statistics for episode lengths (sense check)
#   yearly_phase_summary      : one row per patient per year per phase
#   phase_counts_per_patient  : number of patient-years per phase (reported in results)
#   phase_lme_model_summary   : summary of the before/during/after LME model
#   phase_lme_model_ci        : 95% confidence intervals for model coefficients
# =============================================================================


# =============================================================================
# Step 1: Restrict to creatinine tests and attach lithium phase labels
# =============================================================================
# Phase is determined by comparing each test date against the episode start
# and end dates joined onto all_tests_exclusions in Script 3.

creatinine_single_episode <- all_tests_exclusions %>%
  filter(test == "Creatinine")

phase_data <- creatinine_single_episode %>%
  mutate(
    lithium_phase = case_when(
      collection_date < episode_start                                   ~ "before",
      collection_date >= episode_start & collection_date <= episode_end ~ "during",
      collection_date > episode_end                                     ~ "after",
      TRUE ~ NA_character_
    ),
    # Factor with before as reference level for the LME model
    lithium_phase = factor(lithium_phase, levels = c("before", "during", "after"))
  )


# =============================================================================
# Step 2: Episode length summary (sense check — not reported in results)
# =============================================================================

episode_length_stats <- phase_data %>%
  group_by(study_id) %>%
  summarise(episode_length = first(episode_length_years), .groups = "drop") %>%
  summarise(
    mean_episode_length = mean(episode_length, na.rm = TRUE),
    sd_episode_length   = sd(episode_length,   na.rm = TRUE),
    iqr_episode_length  = IQR(episode_length,  na.rm = TRUE)
  )

episode_length_stats


# =============================================================================
# Step 3: Restrict to patients with a lithium episode of >= 5 years
# =============================================================================
# Patients with shorter episodes are unlikely to have sufficient data in the
# before and after phases to contribute meaningfully to those model estimates.

phase_data_five_yr <- phase_data %>%
  filter(episode_length_years >= 5)


# =============================================================================
# Step 4: Calculate rolling 2-year eGFR slope at each creatinine test
# =============================================================================
# As in Script 4, slopes are computed using calculate_local_gfr_slope()
# and values beyond ±10 ml/min/1.73m2/year are excluded.

phase_slopes <- phase_data_five_yr %>%
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
  mutate(
    egfr_slope_2yr = if_else(abs(egfr_slope_2yr) > 10, NA_real_, egfr_slope_2yr)
  )


# =============================================================================
# Step 5: Collapse to one median slope per patient per year per phase
# =============================================================================
# year_bucket groups each test into an integer year within its phase:
#   - before/during: years since episode start
#   - after: years since episode end
# This ensures the grouping is clinically meaningful within each phase.

yearly_phase_summary <- phase_slopes %>%
  mutate(
    year_bucket = case_when(
      lithium_phase == "before" ~ floor(as.numeric(collection_date - episode_start) / 365.25),
      lithium_phase == "during" ~ floor(as.numeric(collection_date - episode_start) / 365.25),
      lithium_phase == "after"  ~ floor(as.numeric(collection_date - episode_end)   / 365.25)
    )
  ) %>%
  group_by(study_id, lithium_phase, year_bucket) %>%
  summarise(
    slope_year     = median(egfr_slope_2yr, na.rm = TRUE),
    age_year       = median(age_at_test,    na.rm = TRUE),
    sex            = first(sex),
    n_measurements = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(slope_year))


# =============================================================================
# Step 6: Patient count and mean follow-up per phase (reported in results)
# =============================================================================

phase_counts_per_patient <- yearly_phase_summary %>%
  group_by(lithium_phase) %>%
  summarise(
    n_patient_years = n(),
    n_patients      = n_distinct(study_id),
    mean_years      = n() / n_distinct(study_id),
    .groups = "drop"
  )

phase_counts_per_patient


# =============================================================================
# Step 7: Fit before/during/after LME model
# =============================================================================
# "before" is the reference phase (set by factor level ordering in Step 1).
# Model coefficients therefore represent the change in mean eGFR slope
# during and after lithium, relative to before lithium.
#
# Fixed effects:
#   lithium_phase : before (reference), during, after
#   age_year      : age in years (covariate)
#   sex           : male/female (covariate)
#
# Random effect: (1 | study_id) = patient-level intercept
# Optimizer: bobyqa used to improve convergence for this larger model

phase_lme_model <- lmer(
  slope_year ~ lithium_phase + age_year + sex + (1 | study_id),
  data    = yearly_phase_summary,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

phase_lme_model_summary <- summary(phase_lme_model)
phase_lme_model_summary

# Profile 95% confidence intervals
phase_lme_model_ci <- confint(phase_lme_model, method = "profile")
phase_lme_model_ci


# =============================================================================
# End of Script 6
# Outputs: phase_data, episode_length_stats, yearly_phase_summary,
#          phase_counts_per_patient, phase_lme_model_summary, phase_lme_model_ci
# =============================================================================
