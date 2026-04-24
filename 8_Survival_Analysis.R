# =============================================================================
# Script 8: Survival Analysis — Time to Lithium Discontinuation by eGFR Group
# =============================================================================
# Investigates whether kidney function thresholds are associated with the
# likelihood of lithium discontinuation, using two complementary approaches:
#
#   1. Kaplan-Meier (KM) curves: visualise the probability of remaining on
#      lithium over time for three groups stratified by eGFR threshold crossing.
#
#   2. Time-varying Cox proportional hazards model: quantifies the hazard ratio
#      for lithium discontinuation associated with eGFR category, allowing
#      eGFR group to change over the course of follow-up.
#
# Three KM groups:
#   "eGFR >= 60"           : patients whose eGFR remained >= 60 throughout;
#                            follow-up starts at their first eGFR observation
#                            and is censored when (if) eGFR drops below 60.
#   "eGFR <60 (never <30)" : patients who crossed below 60 but never below 30;
#                            follow-up starts at the date of first crossing.
#   "eGFR <30"             : patients who crossed below 30;
#                            follow-up starts at the date of first crossing.
#
# Only patients on lithium at the time of their group entry date are included.
# The event of interest is the first date classified as "stopping lithium".
# Patients are censored at date of death or last recorded observation.
#
# Requires: Script 4 (gfr_slopes with egfr, lithium_status, dod columns)
#
# Outputs:
#   thresholds      : per-patient eGFR threshold crossing dates
#   km_df           : combined KM dataset (all three groups)
#   km_fit          : survfit object for KM plot (used in summary figures)
#   cox_tv          : time-varying Cox model object
# =============================================================================


# =============================================================================
# Step 1: Build ancillary lookup tables
# =============================================================================

# Last observed creatinine date per patient — used as censoring date if no
# discontinuation or death occurs during follow-up
last_observation_date <- gfr_slopes %>%
  filter(!is.na(collection_date)) %>%
  group_by(study_id) %>%
  summarise(last_obs_date = max(collection_date, na.rm = TRUE), .groups = "drop")

# Date of death per patient (if recorded)
death_date_lookup <- gfr_slopes %>%
  select(study_id, dod) %>%
  distinct()

# Date lithium was first classified as "stopping" per patient — used as the
# event date if the patient discontinued lithium during follow-up
lithium_stop_date_lookup <- gfr_slopes %>%
  filter(lithium_status == "stopping lithium") %>%
  arrange(study_id, collection_date) %>%
  group_by(study_id) %>%
  summarise(stop_date = first(collection_date), .groups = "drop")

# Per-patient dates of first eGFR crossing below 60 and below 30,
# plus first observation date (used as group entry for the >= 60 reference group)
egfr_threshold_dates <- gfr_slopes %>%
  filter(!is.na(egfr)) %>%
  arrange(study_id, collection_date) %>%
  group_by(study_id) %>%
  summarise(
    first_below60_date = collection_date[which(egfr < 60)[1]],
    first_below30_date = collection_date[which(egfr < 30)[1]],
    first_obs_date     = first(collection_date),
    .groups = "drop"
  ) %>%
  mutate(
    # Flag patients whose eGFR never fell below 30 (used to define the
    # <60 group as "never <30")
    never_below_30 = is.na(first_below30_date)
  )


# =============================================================================
# Step 2: Helper function to build a KM cohort from a threshold crossing date
# =============================================================================
# For the <60 and <30 groups, follow-up begins at the threshold crossing date.
# Only patients who were still on lithium at that date are included.
# Events and censoring are defined relative to that threshold date.

build_km_cohort <- function(threshold_col, group_label,
                             restrict_never_below30 = FALSE) {

  df <- egfr_threshold_dates %>%
    filter(!is.na(.data[[threshold_col]]))

  # For the <60 group, exclude patients who also reached <30
  # (they form the separate <30 group)
  if (restrict_never_below30) {
    df <- df %>% filter(never_below_30)
  }

  df %>%
    left_join(lithium_stop_date_lookup, by = "study_id") %>%
    left_join(death_date_lookup,        by = "study_id") %>%
    left_join(last_observation_date,    by = "study_id") %>%
    mutate(threshold_date = .data[[threshold_col]]) %>%

    # Exclude patients who had already stopped lithium before their threshold date
    filter(is.na(stop_date) | stop_date > threshold_date) %>%

    mutate(
      # Stopping event: first "stopping lithium" date after the threshold
      stop_after = if_else(
        !is.na(stop_date) & stop_date > threshold_date,
        stop_date, as.Date(NA)
      ),
      event = as.integer(!is.na(stop_after)),

      # Event date: stopping date if event occurred; death or last obs otherwise
      event_date = case_when(
        !is.na(stop_after)                  ~ stop_after,
        !is.na(dod) & dod > threshold_date  ~ dod,
        TRUE                                ~ last_obs_date
      ),
      time_years = as.numeric(
        difftime(event_date, threshold_date, units = "days")
      ) / 365.25,
      group = group_label
    ) %>%
    filter(time_years > 0)
}


# =============================================================================
# Step 3: Build each of the three KM cohorts
# =============================================================================

# Reference group: patients whose eGFR remained >= 60 throughout.
# Follow-up starts at first eGFR observation; censored when eGFR drops below 60
# (or at last observation/death if eGFR never drops below 60).
km_above60 <- egfr_threshold_dates %>%
  left_join(lithium_stop_date_lookup, by = "study_id") %>%
  left_join(death_date_lookup,        by = "study_id") %>%
  left_join(last_observation_date,    by = "study_id") %>%
  filter(!is.na(first_obs_date)) %>%
  # Restrict to patients who had at least one eGFR >= 60
  semi_join(gfr_slopes %>% filter(egfr >= 60), by = "study_id") %>%
  # Only include patients on lithium at first observation
  filter(is.na(stop_date) | stop_date > first_obs_date) %>%
  mutate(
    threshold_date = first_obs_date,
    # Censor at the date eGFR first dropped below 60, or last observation
    censor_date = case_when(
      !is.na(first_below60_date) ~ first_below60_date,
      TRUE                       ~ last_obs_date
    ),
    stop_after = if_else(
      !is.na(stop_date) & stop_date > threshold_date & stop_date <= censor_date,
      stop_date, as.Date(NA)
    ),
    event = as.integer(!is.na(stop_after)),
    event_date = case_when(
      !is.na(stop_after)                                          ~ stop_after,
      !is.na(dod) & dod > threshold_date & dod <= censor_date    ~ dod,
      TRUE                                                        ~ censor_date
    ),
    time_years = as.numeric(
      difftime(event_date, threshold_date, units = "days")
    ) / 365.25,
    group = "eGFR \u226560"
  ) %>%
  filter(time_years > 0)

# eGFR <60 group: patients who crossed below 60 but never below 30
km_below60 <- build_km_cohort(
  threshold_col         = "first_below60_date",
  group_label           = "eGFR <60 (never <30)",
  restrict_never_below30 = TRUE
)

# eGFR <30 group: patients who crossed below 30
km_below30 <- build_km_cohort(
  threshold_col = "first_below30_date",
  group_label   = "eGFR <30"
)

# Combine all three groups into a single dataset for survival analysis
km_df <- bind_rows(km_above60, km_below60, km_below30) %>%
  mutate(
    group = factor(group,
                   levels = c("eGFR \u226560",
                              "eGFR <60 (never <30)",
                              "eGFR <30"))
  )


# =============================================================================
# Step 4: Kaplan-Meier survival fit and plot
# =============================================================================
# km_fit is also used in summary_of_r_figures.R to produce the final figure 7.

km_fit <- survfit(Surv(time_years, event) ~ group, data = km_df)

# Log-rank test for difference between groups (reported in results)
survdiff(Surv(time_years, event) ~ group, data = km_df)

# KM plot (interim version; final formatted version is in summary_of_r_figures.R)
ggsurvplot(
  km_fit,
  data              = km_df,
  risk.table        = TRUE,
  pval              = TRUE,
  xlim              = c(0, 10),
  break.time.by     = 2,
  xlab              = "Years from eGFR threshold crossing",
  ylab              = "Probability of remaining on lithium",
  legend.title      = "eGFR group",
  legend.labs       = c("eGFR \u226560", "eGFR <60 (never <30)", "eGFR <30"),
  ggtheme           = theme_classic(),
  size              = 1.2,
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  palette           = c("#00A087", "#4DBBD5", "#E64B35"),
  censor            = FALSE
)


# =============================================================================
# Step 5: Time-varying Cox proportional hazards model
# =============================================================================
# Splits each patient's follow-up at the dates their eGFR first crossed
# below 60 and below 30, creating a long-format dataset where eGFR_group
# acts as a time-varying covariate. This allows the model to account for
# the fact that a patient's eGFR category can change over time.

long_format_df <- egfr_threshold_dates %>%
  left_join(lithium_stop_date_lookup, by = "study_id") %>%
  left_join(death_date_lookup,        by = "study_id") %>%
  left_join(last_observation_date,    by = "study_id") %>%
  mutate(
    # Final date of follow-up for each patient
    final_date = case_when(
      !is.na(stop_date) ~ stop_date,
      !is.na(dod)       ~ dod,
      TRUE              ~ last_obs_date
    ),
    event = as.integer(!is.na(stop_date))
  )

# Split each patient's follow-up into intervals at each threshold crossing date.
# Within each interval, the patient is assigned to the eGFR group that applied
# at the start of that interval.
long_format_df <- long_format_df %>%
  split(.$study_id) %>%
  map_dfr(function(row) {
    row <- row[1, ]

    # The boundaries of each follow-up interval for this patient
    boundaries <- sort(unique(na.omit(as.Date(c(
      row$first_obs_date,
      row$first_below60_date,
      row$first_below30_date,
      row$final_date
    )))))

    if (length(boundaries) < 2) return(NULL)

    tibble(
      study_id = row$study_id,
      start    = head(boundaries, -1),
      stop     = tail(boundaries, -1)
    ) %>%
      mutate(
        # Assign eGFR group based on which thresholds have been crossed
        # by the start of each interval
        egfr_group = case_when(
          !is.na(row$first_below30_date) & start >= row$first_below30_date ~ "eGFR <30",
          !is.na(row$first_below60_date) & start >= row$first_below60_date ~ "eGFR <60",
          TRUE ~ "eGFR \u226560"
        ),
        event      = if_else(stop == row$final_date & row$event == 1L, 1L, 0L),
        baseline   = row$first_obs_date,
        start_time = as.numeric(difftime(start, baseline, units = "days")) / 365.25,
        stop_time  = as.numeric(difftime(stop,  baseline, units = "days")) / 365.25
      )
  }) %>%
  filter(stop_time > start_time) %>%
  mutate(
    egfr_group = factor(egfr_group,
                        levels = c("eGFR \u226560", "eGFR <60", "eGFR <30"))
  )

# Fit the time-varying Cox model
# Reference category: eGFR >= 60
cox_time_varying <- coxph(
  Surv(start_time, stop_time, event) ~ egfr_group,
  data = long_format_df
)

summary(cox_time_varying)

# =============================================================================
# Step 6: Proportional hazards assumption check
# =============================================================================
# Schoenfeld residuals test: a non-significant result supports the assumption
# that the hazard ratio between groups is constant over time.

cox.zph(cox_time_varying)
ggcoxzph(cox.zph(cox_time_varying))


# =============================================================================
# End of Script 8
# Outputs: egfr_threshold_dates, km_df, km_fit, cox_time_varying
# =============================================================================
