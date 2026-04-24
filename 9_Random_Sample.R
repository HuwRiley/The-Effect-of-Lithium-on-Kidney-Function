# =============================================================================
# Script 9: Random Sample for Manual Record Review
# =============================================================================
# Generates the sample of 100 patients used in the exploratory manual review
# of NHS TrakCare electronic patient records to identify documented reasons
# for lithium discontinuation.
#
# Sampling approach:
#   All patients classified as "stopping lithium" at least once during
#   follow-up are randomly ordered using a fixed seed (set.seed(18)) for
#   reproducibility. The first 100 patients in this randomly ordered list
#   are selected. This matches the method described in the dissertation.
#
# The output CSV includes CHI numbers to allow TrakCare lookup and is
# therefore NOT included in the public GitHub repository.
#
# After manual review, reasons for stopping are entered into the
# reason_for_stopping column of the CSV and results are visualised
# using the hardcoded counts at the bottom of this script.
#
# Requires: Scripts 3 (episode_numbers) and 7 (stopping_egfr_patient)
#
# Output:
#   random_stopped_lithium.csv  : DO NOT COMMIT — contains CHI numbers
# =============================================================================

set.seed(18)  # fixed seed ensures the same 100 patients are selected each time


# =============================================================================
# Step 1: Retrieve CHI number for each study_id
# =============================================================================
# CHI numbers are needed to look up patients in NHS TrakCare.
# They are pulled from all_tests (which retains CHI before anonymisation)
# and will be removed from any publicly shared version of the output.

chi_lookup <- all_tests %>%
  group_by(study_id) %>%
  summarise(chi = first(chi), .groups = "drop")


# =============================================================================
# Step 2: Build wide-format episode date table
# =============================================================================
# For each patient, episode start and end dates are pivoted to wide format
# so that each patient appears on one row with columns for each episode.
# This gives reviewers the context needed to assess TrakCare records.

episode_dates_wide <- episode_numbers %>%
  group_by(study_id, episode_id) %>%
  summarise(
    episode_start        = min(collection_date),
    episode_end          = max(collection_date),
    episode_length_years = as.numeric(
      max(collection_date) - min(collection_date) + 365.25
    ) / 365.25,
    .groups = "drop"
  ) %>%
  group_by(study_id) %>%
  mutate(episode_number = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    id_cols     = study_id,
    names_from  = episode_number,
    values_from = c(episode_start, episode_end, episode_length_years),
    names_glue  = "episode{episode_number}_{.value}"
  )


# =============================================================================
# Step 3: Build the sampling frame and randomly order it
# =============================================================================
# The sampling frame is all patients whose last observed lithium status was
# "stopping lithium" (i.e. stopping_egfr_patient from Script 7).

random_stopped_sample <- stopping_egfr_patient %>%
  select(study_id, collection_date, egfr, egfr_slope_2yr, ckd_stage) %>%
  left_join(episode_dates_wide, by = "study_id") %>%
  left_join(chi_lookup,         by = "study_id") %>%
  rename(
    lithium_stop_date     = collection_date,
    egfr_at_stopping      = egfr,
    egfr_slope_at_stopping = egfr_slope_2yr
  ) %>%
  # Add empty column for reviewers to fill in during TrakCare review
  mutate(reason_for_stopping = NA_character_) %>%
  # Arrange columns: identifiers first, then episode dates, then clinical data
  select(
    chi,
    reason_for_stopping,
    starts_with("episode1_episode_start"), starts_with("episode1_episode_end"),
    starts_with("episode2_episode_start"), starts_with("episode2_episode_end"),
    starts_with("episode3_episode_start"), starts_with("episode3_episode_end"),
    starts_with("episode4_episode_start"), starts_with("episode4_episode_end"),
    starts_with("episode5_episode_start"), starts_with("episode5_episode_end"),
    everything()
  ) %>%
  # Randomly shuffle all rows (set.seed above ensures reproducibility)
  slice_sample(prop = 1)

# Select the first 100 rows from the randomly ordered dataset
# This is equivalent to a simple random sample of 100 patients
sample_for_review <- head(random_stopped_sample, 100)

# Export to CSV for TrakCare review
# WARNING: this file contains CHI numbers — do not commit to GitHub
write.csv(sample_for_review, "random_stopped_lithium.csv", row.names = FALSE)


# =============================================================================
# Step 4: Visualise reasons for discontinuation
# =============================================================================
# Counts are hardcoded from the completed manual review of 100 records.
# Categories were defined prior to review (see dissertation methods).

review_results <- data.frame(
  reason = c(
    "Did not stop lithium",
    "Unknown",
    "Non-renal side effects",
    "Renal side effects",
    "Unsuitable medication",
    "Lithium toxicity\n(unspecified cause)"
  ),
  n = c(36, 25, 12, 15, 8, 4)
) %>%
  mutate(
    percent = (n / sum(n)) * 100,
    # Order by ascending percentage for horizontal bar chart
    reason  = factor(reason, levels = reason[order(percent)])
  )

# Horizontal bar chart — interim version; final formatted version is
# figure_8 in summary_of_r_figures.R
ggplot(review_results, aes(x = reason, y = percent)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = paste0("n = ", n)),
    hjust = -0.1,
    size  = 4
  ) +
  coord_flip() +
  labs(
    title = "Reasons for Lithium Discontinuation (n = 100)",
    x     = NULL,
    y     = "Percentage (%)"
  ) +
  ylim(0, max(review_results$percent) + 5) +
  theme_minimal() +
  theme(
    plot.title  = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 11)
  )


# =============================================================================
# End of Script 9
# Output: random_stopped_lithium.csv (DO NOT COMMIT — contains CHI numbers)
# =============================================================================
