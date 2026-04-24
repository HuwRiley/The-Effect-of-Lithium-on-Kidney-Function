# =============================================================================
# Script 3: Cohort Exclusions and Episode Classification
# =============================================================================
# Applies two sequential exclusion steps to define the analysis cohorts:
#
#   Exclusion 1: Remove patients with no positive lithium test (>0.3 mmol/L)
#                during follow-up. This threshold is slightly below the
#                therapeutic range to maximise sensitivity for detecting
#                ongoing lithium use while excluding subtherapeutic levels.
#
#   Exclusion 2: Identify lithium treatment episodes for each patient.
#                A new episode is defined when a gap of >365 days occurs
#                between consecutive positive lithium tests. Episode start and
#                end dates and durations are calculated for each episode.
#
# Requires: Script 1 (functions) and Script 2 (all_tests)
#
# Outputs:
#   all_tests_with_lithium  : all tests for patients with >= 1 positive result
#                             (n = 1,045; used in scripts 4, 5, and 8)
#   all_tests_exclusions    : all tests for single-episode patients only,
#                             with episode dates joined on
#                             (n = 513 after 5-year filter in script 6)
#   episode_numbers         : lithium test-level episode labels (used in
#                             Random_Sample.R)
# =============================================================================


# =============================================================================
# Exclusion 1: Retain only patients with at least one positive lithium test
# =============================================================================
# A positive test is defined as a serum lithium level > 0.3 mmol/L.
# 1,069 patients had any non-zero lithium test in 2012; 1,045 had at least
# one result above the 0.3 mmol/L threshold.

all_tests_with_lithium <- all_tests %>%
  group_by(study_id) %>%
  filter(any(test == "Lithium" & test_data > 0.3)) %>%
  ungroup()



# =============================================================================
# Exclusion 2: Identify and label lithium treatment episodes
# =============================================================================
# For each patient, consecutive positive lithium tests are grouped into
# episodes. A gap of more than 365 days between tests is treated as the
# start of a new episode. This allows us to identify patients who stopped
# and restarted lithium during follow-up.

# Label each positive lithium test with its episode number
episode_numbers <- all_tests_with_lithium %>%
  filter(test == "Lithium", test_data > 0.3) %>%
  group_by(study_id) %>%
  arrange(collection_date) %>%
  mutate(
    # Days since the previous lithium test for this patient
    gap_days    = collection_date - lag(collection_date),
    # Flag as new episode if gap > 365 days or first test for this patient
    new_episode = if_else(is.na(gap_days) | gap_days > 365, 1, 0),
    # Cumulative sum of new_episode flags gives the episode number
    episode_id  = cumsum(new_episode)
  )


# -----------------------------------------------------------------------------
# Count episodes per patient and identify single-episode patients
# -----------------------------------------------------------------------------
# For the before/during/after analysis (Script 6), only patients with a single
# uninterrupted lithium episode can be classified into clear phases.

episode_counts <- episode_numbers %>%
  group_by(study_id) %>%
  summarise(
    n_episodes = n_distinct(episode_id),
    .groups = "drop"
  )

# study_ids of patients who had exactly one treatment episode
single_episode_patient_ids <- episode_counts %>%
  filter(n_episodes == 1)

# Retain only single-episode patients in the exclusions dataset
all_tests_exclusions <- all_tests_with_lithium %>%
  filter(study_id %in% single_episode_patient_ids$study_id)


# -----------------------------------------------------------------------------
# Calculate episode start date, end date, and duration
# -----------------------------------------------------------------------------
# Episode length is calculated as the time from first to last positive test,
# plus 365 days to account for the ongoing treatment period after the final
# recorded test (reflecting the ±1 year exposure window used in Script 4).

episode_lengths <- episode_numbers %>%
  group_by(study_id, episode_id) %>%
  summarise(
    episode_start        = min(collection_date),
    episode_end          = max(collection_date),
    episode_length_days  = as.numeric((episode_end - episode_start) + 365.25),
    episode_length_years = episode_length_days / 365.25,
    .groups = "drop"
  ) %>%
  select(study_id, episode_start, episode_end,
         episode_length_days, episode_length_years)

# Join episode dates onto the single-episode patient dataset
all_tests_exclusions <- all_tests_exclusions %>%
  left_join(episode_lengths, by = "study_id")


# =============================================================================
# End of Script 3
# Outputs: all_tests_with_lithium, all_tests_exclusions, episode_numbers
# =============================================================================
