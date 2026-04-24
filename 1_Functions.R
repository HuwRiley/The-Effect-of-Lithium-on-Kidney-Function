# =============================================================================
# Script 1: Helper Functions
# =============================================================================
# Defines three functions used throughout the analysis pipeline:
#   - calculate_egfr()         : computes eGFR from creatinine, age, and sex
#   - calculate_local_gfr_slope() : fits a rolling 2-year eGFR slope
#   - calculate_cumulative_lithium_time() : sums lithium exposure up to a date
#
# No inputs required. Run this script first before any other script.
# =============================================================================


# -----------------------------------------------------------------------------
# Function 1: calculate_egfr
# -----------------------------------------------------------------------------
# Calculates estimated Glomerular Filtration Rate (eGFR) using the CKD-EPI
# equation (Levey et al., 2009). This equation was chosen over MDRD as it
# is more accurate at higher (healthier) levels of kidney function.
#
# Ethnicity is not included in line with current recommendations to remove
# race-based adjustments from eGFR estimation.
#
# Arguments:
#   sex : character, "male" or "female"
#   age : numeric, age in years at time of test
#   scr : numeric, serum creatinine in mg/dL
#         (note: NHS Lothian reports creatinine in umol/L; divide by 88.4
#          before passing to this function)
#
# Returns: numeric eGFR in ml/min/1.73m2
# -----------------------------------------------------------------------------
calculate_egfr <- function(sex, age, scr) {

  # Sex-specific constants from the CKD-EPI equation
  kappa <- case_when(
    sex == "female" ~ 0.7,
    sex == "male"   ~ 0.9
  )

  alpha <- case_when(
    sex == "female" ~ -0.329,
    sex == "male"   ~ -0.411
  )

  sex_constant <- case_when(
    sex == "female" ~ 1.018,
    sex == "male"   ~ 1.000
  )

  egfr <- 141 *
    (pmin(scr / kappa, 1) ^ alpha) *
    (pmax(scr / kappa, 1) ^ -1.209) *
    (0.993 ^ age) *
    sex_constant

  return(egfr)
}


# -----------------------------------------------------------------------------
# Function 2: calculate_local_gfr_slope
# -----------------------------------------------------------------------------
# Fits a linear regression to all eGFR measurements within a rolling 2-year
# window ending at a given centre date, and returns the slope (rate of change
# in eGFR per year).
#
# Requires at least 3 measurements in the window to fit a reliable slope;
# returns NA if fewer are available.
#
# Arguments:
#   dates       : vector of Date values for all measurements for a patient
#   gfr         : numeric vector of corresponding eGFR values
#   centre_date : Date, the measurement date being evaluated (right edge of window)
#   window_days : numeric, size of the lookback window in days (default = 730.5,
#                 i.e. 2 years)
#
# Returns: numeric slope in ml/min/1.73m2/year, or NA if insufficient data
# -----------------------------------------------------------------------------
calculate_local_gfr_slope <- function(dates, gfr, centre_date, window_days = 730.5) {

  # Identify measurements that fall within the 2-year window
  in_window <- dates > (centre_date - window_days) & dates <= centre_date
  has_egfr  <- in_window & !is.na(gfr)

  # Require at least 3 measurements for a meaningful slope
  if (sum(has_egfr) < 3) return(NA_real_)

  # Convert dates to time in years relative to the centre date
  time_years <- as.numeric(dates[has_egfr] - centre_date) / 365.25

  # Fit simple linear regression: eGFR ~ time
  model <- lm(gfr[has_egfr] ~ time_years)
  slope <- coef(model)[2]

  return(slope)
}


# -----------------------------------------------------------------------------
# Function 3: calculate_cumulative_lithium_time
# -----------------------------------------------------------------------------
# Calculates the total number of years a patient has been exposed to lithium
# up to a given date, based on their lithium treatment intervals.
#
# Treatment intervals are defined in script 5 as date ranges during which
# a patient had active lithium levels (>0.3 mmol/L), extended by 365 days
# past the final positive test to account for ongoing treatment.
#
# Arguments:
#   patient_id : numeric study_id for the patient
#   gfr_date   : Date, the eGFR measurement date to calculate exposure up to
#   intervals  : data frame with columns study_id, start_date, end_date
#                (one row per treatment interval per patient)
#
# Returns: numeric cumulative years of lithium exposure up to gfr_date
# -----------------------------------------------------------------------------
calculate_cumulative_lithium_time <- function(patient_id, gfr_date, intervals) {

  # Extract this patient's treatment intervals
  patient_intervals <- intervals %>% filter(study_id == patient_id)

  # Return 0 if no intervals exist for this patient
  if (nrow(patient_intervals) == 0) return(0)

  total_days <- 0

  for (i in seq_len(nrow(patient_intervals))) {

    interval_start <- patient_intervals$start_date[i]
    interval_end   <- patient_intervals$end_date[i]

    # Only count exposure that occurred before the eGFR measurement date
    if (gfr_date > interval_start) {

      # The overlap ends at whichever comes first: interval end or eGFR date
      overlap_end <- min(interval_end, gfr_date)

      if (overlap_end > interval_start) {
        total_days <- total_days + as.numeric(overlap_end - interval_start)
      }
    }
  }

  # Convert days to years
  return(total_days / 365.25)
}
