# =============================================================================
# Script 10: Dissertation Figures (Figures 3–8)
# =============================================================================
# Produces all figures included in the dissertation results section using a
# consistent Lancet colour palette (ggsci::pal_lancet). All packages are
# loaded in Script 2 (2_Data_Setup.R).
#
# Figures produced:
#   Figure 3 : Age distribution of cohort at January 2012 (histogram)
#   Figure 4 : eGFR slope by lithium status category (boxplot)
#   Figure 5 : CKD stage distribution at lithium discontinuation (flextable)
#   Figure 6 : eGFR at lithium discontinuation (histogram with CKD thresholds)
#   Figure 7 : Kaplan-Meier plot — time to discontinuation by eGFR group
#   Figure 8 : Reasons for lithium discontinuation (horizontal bar chart)
#
# Source objects required (run scripts in order 1 → 9 before this script):
#   all_tests_with_lithium  (Script 3) — Figure 3
#   gfr_slopes              (Script 4) — Figure 4
#   stopping_egfr_patient   (Script 7) — Figures 5 and 6
#   km_fit, km_df           (Script 8) — Figure 7
#   hardcoded review counts            — Figure 8
# =============================================================================

# Lancet colour palette: 8 colours used consistently across all figures
# #00468B #ED0000 #42B540 #0099B4 #925E9F #FDAF91 #AD002A #ADB6B6
lancet_cols <- pal_lancet("lanonc")(8)


# =============================================================================
# Figure 3: Age Distribution of Cohort at January 2012
# =============================================================================
# One row per patient (slice(1) within each study_id) to avoid double-counting.
# age_at_2012 was derived from CHI in Script 2.

cohort_one_row_per_patient <- all_tests_with_lithium %>%
  group_by(study_id) %>%
  slice(1) %>%
  ungroup()

age_distribution_stats <- cohort_one_row_per_patient %>%
  summarise(
    n        = n(),
    mean_age = mean(age_at_2012, na.rm = TRUE),
    sd_age   = sd(age_at_2012,   na.rm = TRUE),
    min_age  = min(age_at_2012,  na.rm = TRUE),
    max_age  = max(age_at_2012,  na.rm = TRUE)
  )

figure_3 <- ggplot(cohort_one_row_per_patient, aes(x = age_at_2012)) +
  geom_histogram(
    binwidth  = 5,
    fill      = lancet_cols[1],
    colour    = "white",
    linewidth = 0.3
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_x_continuous(breaks = seq(20, 90, by = 10)) +
  labs(
    title    = "Age Distribution of Cohort at January 2012",
    subtitle = sprintf(
      "n = %d  |  Mean = %.1f years (SD = %.1f)  |  Range: %d\u2013%d years",
      age_distribution_stats$n,
      age_distribution_stats$mean_age,
      age_distribution_stats$sd_age,
      age_distribution_stats$min_age,
      age_distribution_stats$max_age
    ),
    x = "Age at 1 January 2012 (years)",
    y = "Number of patients"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 15, colour = "black"),
    plot.subtitle = element_text(size = 12, colour = "grey30", margin = margin(b = 8)),
    axis.title    = element_text(face = "bold", colour = "black"),
    axis.text     = element_text(colour = "black"),
    plot.margin   = margin(16, 16, 16, 16)
  )

figure_3


# =============================================================================
# Figure 4: eGFR Slope by Lithium Status (Boxplot)
# =============================================================================
# Shows the distribution of 2-year rolling eGFR slopes across the four
# lithium status categories. Outliers beyond ±10 ml/min/1.73m2/year were
# already excluded in Script 4, so are not shown here.
# n labels are placed at the bottom of each box for clarity.

slope_boxplot_data <- gfr_slopes %>%
  filter(!is.na(egfr_slope_2yr), !is.na(lithium_status)) %>%
  mutate(
    lithium_status = factor(
      lithium_status,
      levels = c("not on lithium", "starting lithium",
                 "on lithium",     "stopping lithium"),
      labels = c("Not on lithium", "Starting lithium",
                 "On lithium",     "Stopping lithium")
    )
  )

# n per group for labels at the bottom of the plot
slope_boxplot_n <- slope_boxplot_data %>%
  group_by(lithium_status) %>%
  summarise(n = n(), .groups = "drop")

figure_4 <- ggplot(
  slope_boxplot_data,
  aes(x = lithium_status, y = egfr_slope_2yr, fill = lithium_status)
) +
  # Horizontal reference line at slope = 0
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_boxplot(
    outlier.shape = NA,
    width         = 0.55,
    linewidth     = 0.5,
    alpha         = 0.85
  ) +
  scale_fill_manual(values = c(
    "Not on lithium"   = lancet_cols[8],
    "Starting lithium" = lancet_cols[3],
    "On lithium"       = lancet_cols[1],
    "Stopping lithium" = lancet_cols[2]
  )) +
  scale_y_continuous(
    limits = c(-10, 10),
    breaks = seq(-10, 10, by = 2),
    # Prefix positive values with "+" for clarity
    labels = function(x) ifelse(x > 0, paste0("+", x), x)
  ) +
  # n labels at the bottom of each box
  geom_text(
    data        = slope_boxplot_n,
    aes(x = lithium_status, y = -9.5,
        label = paste0("n = ", comma(n))),
    inherit.aes = FALSE,
    size        = 3.8,
    colour      = "grey30"
  ) +
  labs(
    title    = "eGFR Slope by Lithium Status",
    subtitle = paste0("Two-year rolling eGFR slope at each measurement;",
                      " slopes beyond \u00b110 ml/min/1.73m\u00b2/year excluded"),
    x        = "Lithium status",
    y        = "eGFR slope (ml/min/1.73m\u00b2/year)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 15, colour = "black"),
    plot.subtitle   = element_text(size = 11, colour = "grey30", margin = margin(b = 8)),
    axis.title      = element_text(face = "bold", colour = "black"),
    axis.text       = element_text(colour = "black"),
    axis.text.x     = element_text(size = 12),
    legend.position = "none",
    plot.margin     = margin(16, 16, 16, 16)
  )

figure_4


# =============================================================================
# Figure 5: CKD Stage Distribution at Lithium Discontinuation (Table)
# =============================================================================
# Reformats the ckd_stage column from stopping_egfr_patient with readable
# labels and builds a styled flextable for inclusion in the dissertation.

ckd_table_data <- stopping_egfr_patient %>%
  mutate(
    ckd_stage = factor(
      ckd_stage,
      levels = c("Normal(>90)", "Mild(60-90)", "Moderate(30-60)", "Severe(<30)"),
      labels = c(
        "Normal (eGFR \u226590)",
        "Mild CKD (eGFR 60\u201390)",
        "Moderate CKD (eGFR 30\u201360)",
        "Severe CKD (eGFR <30)"
      )
    )
  ) %>%
  count(ckd_stage, .drop = FALSE) %>%
  mutate(
    Percentage = paste0(round(n / sum(n) * 100, 1), "%")
  ) %>%
  rename(
    "CKD Stage at Discontinuation" = ckd_stage,
    "N"                            = n,
    "Percentage"                   = Percentage
  )

figure_5 <- flextable(ckd_table_data) %>%
  set_caption(
    caption = "Table 1. CKD stage distribution at time of lithium discontinuation"
  ) %>%
  bold(part = "header") %>%
  # Lancet navy header background
  bg(part = "header", bg = "#00468B") %>%
  color(part = "header", color = "white") %>%
  # Alternating row shading for readability
  bg(i = seq(2, nrow(ckd_table_data), 2), bg = "#F0F4FA", part = "body") %>%
  border_outer(part = "all",  border = officer::fp_border(color = "#AAAAAA", width = 0.5)) %>%
  border_inner_h(part = "body", border = officer::fp_border(color = "#DDDDDD", width = 0.5)) %>%
  align(j = c("N", "Percentage"), align = "center", part = "all") %>%
  width(j = "CKD Stage at Discontinuation", width = 2.8) %>%
  width(j = "N",          width = 0.7) %>%
  width(j = "Percentage", width = 1.0) %>%
  add_footer_lines(
    "Percentages may not sum to exactly 100% due to rounding."
  ) %>%
  color(part = "footer", color = "grey40") %>%
  fontsize(part = "footer", size = 9) %>%
  theme_vanilla() %>%
  autofit()

figure_5


# =============================================================================
# Figure 6: Distribution of eGFR at Lithium Discontinuation (Histogram)
# =============================================================================
# Vertical dashed lines at eGFR 30, 45, and 60 mark the three CKD clinical
# thresholds discussed in the results. The eGFR 45 threshold is included
# because it is referenced in the results text (CKD stage 3b boundary).

egfr_stopping_stats <- stopping_egfr_patient %>%
  summarise(
    n           = n(),
    mean_egfr   = mean(egfr,    na.rm = TRUE),
    median_egfr = median(egfr,  na.rm = TRUE),
    sd_egfr     = sd(egfr,      na.rm = TRUE),
    min_egfr    = min(egfr,     na.rm = TRUE),
    max_egfr    = max(egfr,     na.rm = TRUE)
  )

figure_6 <- ggplot(stopping_egfr_patient, aes(x = egfr)) +
  geom_histogram(
    binwidth  = 5,
    fill      = lancet_cols[1],
    colour    = "white",
    linewidth = 0.3
  ) +
  # CKD clinical threshold lines
  geom_vline(xintercept = 60, linetype = "dashed",
             colour = lancet_cols[4], linewidth = 0.9) +
  geom_vline(xintercept = 45, linetype = "dashed",
             colour = lancet_cols[5], linewidth = 0.9) +
  geom_vline(xintercept = 30, linetype = "dashed",
             colour = lancet_cols[2], linewidth = 0.9) +
  # Threshold annotations — staggered vertically to avoid overlap
  annotate("text", x = 60, y = Inf,
           label    = "eGFR 60\n(CKD stage 3a)",
           vjust    = 1.3, hjust = -0.07, size = 3.3,
           colour   = lancet_cols[4], fontface = "bold") +
  annotate("text", x = 45, y = Inf,
           label    = "eGFR 45\n(CKD stage 3b)",
           vjust    = 3.4, hjust = -0.07, size = 3.3,
           colour   = lancet_cols[5], fontface = "bold") +
  annotate("text", x = 30, y = Inf,
           label    = "eGFR 30\n(CKD stage 4)",
           vjust    = 5.5, hjust = -0.07, size = 3.3,
           colour   = lancet_cols[2], fontface = "bold") +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_x_continuous(breaks = seq(0, 150, by = 15)) +
  labs(
    title    = "Distribution of eGFR at Time of Lithium Discontinuation",
    subtitle = sprintf(
      "n = %d  |  Mean = %.1f (SD = %.1f)  |  Median = %.1f  |  Range: %.1f\u2013%.1f ml/min/1.73m\u00b2",
      egfr_stopping_stats$n,
      egfr_stopping_stats$mean_egfr,
      egfr_stopping_stats$sd_egfr,
      egfr_stopping_stats$median_egfr,
      egfr_stopping_stats$min_egfr,
      egfr_stopping_stats$max_egfr
    ),
    x = "eGFR at discontinuation (ml/min/1.73m\u00b2)",
    y = "Number of patients"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 15, colour = "black"),
    plot.subtitle = element_text(size = 12, colour = "grey30", margin = margin(b = 8)),
    axis.title    = element_text(face = "bold", colour = "black"),
    axis.text     = element_text(colour = "black"),
    plot.margin   = margin(16, 16, 16, 16)
  )

figure_6


# =============================================================================
# Figure 7: Kaplan-Meier Plot — Time to Lithium Discontinuation by eGFR Group
# =============================================================================
# Uses km_fit and km_df produced in Script 8 (8_Survival_Analysis.R).
# Lancet colours: green = reference (>=60), blue = <60, red = <30.

figure_7 <- ggsurvplot(
  km_fit,
  data              = km_df,
  risk.table        = TRUE,
  pval              = TRUE,
  pval.size         = 4,
  pval.coord        = c(0.3, 0.12),
  xlim              = c(0, 10),
  break.time.by     = 2,
  xlab              = "Years from eGFR threshold crossing",
  ylab              = "Probability of remaining on lithium",
  title             = "Time to Lithium Discontinuation by eGFR Threshold Group",
  legend.title      = "eGFR group",
  legend.labs       = c("eGFR \u226560 (reference)", "eGFR <60, never <30", "eGFR <30"),
  palette           = c(lancet_cols[3], lancet_cols[4], lancet_cols[2]),
  size              = 1.0,
  censor            = FALSE,
  risk.table.height = 0.28,
  risk.table.y.text = FALSE,
  risk.table.title  = "Number at risk",
  ggtheme           = theme_classic(base_size = 13),
  fontsize          = 4,
  tables.theme      = theme_classic(base_size = 11),
  risk.table.col    = "strata"
)

# Apply bold title and axis labels consistently with other figures
figure_7$plot <- figure_7$plot +
  theme(
    plot.title   = element_text(face = "bold", size = 14, colour = "black"),
    axis.title   = element_text(face = "bold", colour = "black"),
    axis.text    = element_text(colour = "black"),
    legend.title = element_text(face = "bold")
  )

figure_7


# =============================================================================
# Figure 8: Reasons for Lithium Discontinuation (Horizontal Bar Chart)
# =============================================================================
# Counts are hardcoded from the completed manual review of 100 TrakCare records
# (see Script 9). Categories were pre-defined before review (see methods).
# Bars are ordered by ascending percentage. Labels show both n and percentage.

discontinuation_reasons <- data.frame(
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
    percent = n / sum(n) * 100,
    # Factor levels in ascending order so bars read bottom-to-top on the chart
    reason  = factor(reason, levels = reason[order(percent)])
  )

figure_8 <- ggplot(
  discontinuation_reasons,
  aes(x = reason, y = percent, fill = reason)
) +
  geom_bar(stat = "identity", width = 0.65) +
  geom_text(
    aes(label = sprintf("n = %d (%.0f%%)", n, percent)),
    hjust  = -0.1,
    size   = 4,
    colour = "black"
  ) +
  scale_fill_manual(values = c(
    "Did not stop lithium"                  = lancet_cols[8],
    "Unknown"                               = lancet_cols[8],
    "Non-renal side effects"                = lancet_cols[4],
    "Renal side effects"                    = lancet_cols[2],
    "Unsuitable medication"                 = lancet_cols[3],
    "Lithium toxicity\n(unspecified cause)" = lancet_cols[5]
  )) +
  scale_y_continuous(
    limits = c(0, max(discontinuation_reasons$percent) + 12),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_flip() +
  labs(
    title    = "Documented Reasons for Lithium Discontinuation",
    subtitle = "Exploratory manual review of 100 patients classified as stopping lithium",
    x        = NULL,
    y        = "Percentage of reviewed patients (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 15, colour = "black"),
    plot.subtitle   = element_text(size = 12, colour = "grey30", margin = margin(b = 8)),
    axis.title.x    = element_text(face = "bold", colour = "black"),
    axis.text       = element_text(colour = "black"),
    axis.text.y     = element_text(size = 12),
    legend.position = "none",
    plot.margin     = margin(16, 24, 16, 16)
  )

figure_8


# =============================================================================
# End of Script 10
# =============================================================================
