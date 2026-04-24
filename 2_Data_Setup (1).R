# =============================================================================
# Script 2: Data Setup and Preprocessing
# =============================================================================
# Loads the raw NHS Lothian laboratory dataset and performs all initial
# cleaning and variable derivation required before analysis.
#
# Steps performed:
#   1. Load raw data and remove invalid test results
#   2. Assign anonymous study IDs (replacing CHI numbers)
#   3. Derive sex and date of birth from the CHI number
#   4. Calculate age at time of test and age at the 2012 census date
#   5. Calculate eGFR from each creatinine measurement
#
# Requires: Script 1 (1_Functions.R) to be run first for calculate_egfr()
#
# Input:  lithium_data.csv  (raw laboratory extract; not included in repo
#         due to NHS data governance — contains patient CHI numbers)
# Output: all_tests  (one row per laboratory test per patient)
# =============================================================================

library(tidyverse)   # data wrangling and ggplot2
library(janitor)     # clean_names() for consistent column naming
library(here)        # here() for relative file paths
library(purrr)       # map functions used in scripts 4 and 5
library(lme4)        # linear mixed effects models (scripts 4, 5, 6)
library(lmerTest)    # p-values for lme4 models
library(survival)    # survival analysis objects (script 8)
library(survminer)   # KM plots and ggcoxzph (script 8)
library(flextable)   # formatted tables (scripts 7 and summary figures)
library(ggsci)       # Lancet colour palette for all figures
library(scales)      # comma() for formatted axis labels in figures
library(officer)     # fp_border() for flextable borders in summary figures


# -----------------------------------------------------------------------------
# Step 1: Load raw data and remove invalid test results
# -----------------------------------------------------------------------------
# The raw file contains both lithium and creatinine test results for all
# patients who had a non-zero lithium test recorded in NHS Lothian in 2012.
# We remove any rows where test_data is non-numeric (e.g. "Insufficient
# sample", quoted strings, or asterisks) as these cannot be used in analysis.

all_tests <- read_csv(here("lithium_data.csv")) %>%
  clean_names() %>%
  filter(
    !grepl("[A-Za-z]", test_data),  # remove rows containing letters
    !test_data == '""',              # remove empty quoted strings
    !test_data == "*"                # remove placeholder asterisks
  )


# -----------------------------------------------------------------------------
# Step 2: Assign anonymous study IDs
# -----------------------------------------------------------------------------
# CHI (Community Health Index) numbers are NHS Scotland patient identifiers
# and cannot be shared publicly. Each unique CHI is mapped to an integer
# study_id. All downstream analysis uses study_id only.

study_id_lookup <- all_tests %>%
  distinct(chi) %>%
  arrange(chi) %>%
  mutate(study_id = row_number())

all_tests <- all_tests %>%
  left_join(study_id_lookup, by = "chi")


# -----------------------------------------------------------------------------
# Step 3: Derive sex from the CHI number
# -----------------------------------------------------------------------------
# The 9th digit of the CHI number encodes sex: odd = male, even = female.

all_tests <- all_tests %>%
  mutate(
    sex_digit = as.numeric(substr(chi, 9, 9)),
    sex = case_when(
      sex_digit %% 2 == 0 ~ "female",
      sex_digit %% 2 == 1 ~ "male"
    )
  ) %>%
  select(-sex_digit)


# -----------------------------------------------------------------------------
# Step 4: Derive date of birth from the CHI number
# -----------------------------------------------------------------------------
# The first 6 digits of the CHI encode date of birth as DDMMYY.
# Two-digit years <= 00 are assumed to be 21st century (20xx);
# all others are assumed to be 20th century (19xx).

all_tests <- all_tests %>%
  mutate(
    dob_day_month = substr(chi, 1, 4),
    dob_year_2dig = substr(chi, 5, 6),
    dob_full_year = case_when(
      as.numeric(dob_year_2dig) <= 0 ~ paste0("20", dob_year_2dig),
      as.numeric(dob_year_2dig) >  0 ~ paste0("19", dob_year_2dig)
    ),
    date_of_birth = as.Date(
      paste0(dob_day_month, dob_full_year),
      format = "%d%m%Y"
    )
  ) %>%
  select(-dob_day_month, -dob_year_2dig, -dob_full_year)


# -----------------------------------------------------------------------------
# Step 5: Calculate age at test and age at the 2012 census date
# -----------------------------------------------------------------------------
# age_at_test      : used in eGFR calculation and as a covariate in models
# age_2012_01_01   : used for cohort description (Figure 3)
# The cohort was defined using 2012 as the census year to allow sufficient
# follow-up for both lithium monitoring and creatinine measurements.

all_tests <- all_tests %>%
  mutate(
    collection_date = as.Date(collection_date),
    age_at_test = floor(
      time_length(interval(date_of_birth, collection_date), "years")
    ),
    age_at_2012 = floor(
      time_length(interval(date_of_birth, as.Date("2012-01-01")), "years")
    )
  )


# -----------------------------------------------------------------------------
# Step 6: Calculate eGFR from creatinine measurements
# -----------------------------------------------------------------------------
# eGFR is calculated using the CKD-EPI equation (defined in Script 1).
# NHS Lothian reports creatinine in umol/L; CKD-EPI requires mg/dL,
# so values are divided by 88.4 before passing to calculate_egfr().
# eGFR is set to NA for non-creatinine tests (i.e. lithium tests).

all_tests <- all_tests %>%
  mutate(
    egfr = if_else(
      test == "Creatinine",
      calculate_egfr(
        sex = sex,
        age = age_at_test,
        scr = as.numeric(test_data) / 88.4   # convert umol/L to mg/dL
      ),
      NA_real_
    )
  )


# =============================================================================
# End of Script 2
# Output: all_tests
# =============================================================================
