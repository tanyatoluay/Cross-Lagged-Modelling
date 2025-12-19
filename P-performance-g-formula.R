###############################################################################
# Parametric g-formula — Final Performance Evaluation (Publication-Ready)
#
# Purpose:
#   1. Load completed g-formula simulation outputs
#   2. Extract counterfactual mean outcomes under static interventions
#   3. Compute interventional contrasts (0.25 − 0)
#   4. Produce publication-ready tables and figures
#   5. Save all outputs under methods/g-formula/performance-g-formula/
#
# Assumptions:
#   - g-formula outputs are stored as .rds files under:
#       methods/g-formula/final/s{scenario}/{direction}/treat{value}.rds
#   - Each .rds file contains gfoRmula::gformula() output with columns:
#       "Interv.", "g-form mean", "Mean lower 95% CI", "Mean upper 95% CI"
#   - Interv. == 1 corresponds to the intervention regime
#
###############################################################################

## ---------------------------------------------------------------------------
## 0. Libraries
## ---------------------------------------------------------------------------

library(data.table)
library(tidyverse)

## ---------------------------------------------------------------------------
## 1. Colour palette (fixed; encodes causal direction)
## ---------------------------------------------------------------------------
## NOTE:
## Colour encodes causal direction (forward vs reverse),
## not statistical method. This is intentional and documented in the thesis.

morris_cols <- c(
  "forward" = "#0072B2",  # blue
  "reverse" = "#D55E00"   # vermillion
)

morris_fill <- c(
  "forward" = "#56B4E9",  # light blue
  "reverse" = "#E69F00"   # orange
)

## ---------------------------------------------------------------------------
## 2. Publication theme (used consistently across figures)
## ---------------------------------------------------------------------------

theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(size = 0.3, colour = "grey85"),
    
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10),
    
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    
    plot.title = element_text(
      size = 12, face = "bold", hjust = 0
    )
  )

## ---------------------------------------------------------------------------
## 3. Output directory structure
## ---------------------------------------------------------------------------

perf_dir <- file.path("methods", "g-formula", "performance-g-formula")

dir.create(perf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(perf_dir, "tables"),  showWarnings = FALSE)
dir.create(file.path(perf_dir, "figures"), showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 4. Load all g-formula outputs
## ---------------------------------------------------------------------------

base_dir <- "methods/g-formula/final"

files <- list.files(
  base_dir,
  pattern = "treat.*\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0L) {
  stop("No g-formula result files found. Check directory structure.")
}

## Combine all simulation outputs into a single table
gf_all <- rbindlist(lapply(files, readRDS), fill = TRUE)

## ---------------------------------------------------------------------------
## 5. Extract intervention rows and harmonise columns
## ---------------------------------------------------------------------------
## Each g-formula output contains:
##   Interv. == 0 : natural course
##   Interv. == 1 : intervention
##
## We retain only Interv. == 1 (counterfactual under intervention).

gf <- gf_all %>%
  as_tibble() %>%
  filter(`Interv.` == 1) %>%
  transmute(
    scenario  = factor(scenario),
    direction = factor(direction, levels = c("forward", "reverse")),
    treat     = factor(treat, levels = c("0", "0.25")),
    rep_id,
    mean      = `g-form mean`,
    lcl       = `Mean lower 95% CI`,
    ucl       = `Mean upper 95% CI`
  )

## ---------------------------------------------------------------------------
## 6. Performance table: counterfactual mean outcomes
## ---------------------------------------------------------------------------
## Aggregated across Monte Carlo replicates

gf_summary <- gf %>%
  group_by(scenario, direction, treat) %>%
  summarise(
    mean_est = mean(mean),
    sd_est   = sd(mean),
    lcl      = mean(lcl),
    ucl      = mean(ucl),
    .groups  = "drop"
  )

write.csv(
  gf_summary,
  file.path(perf_dir, "tables", "gformula_counterfactual_means.csv"),
  row.names = FALSE
)

## ---------------------------------------------------------------------------
## 7. Performance table: interventional contrasts
## ---------------------------------------------------------------------------
## Contrast definition:
##   mean(0.25) − mean(0)

gf_contrast <- gf %>%
  group_by(scenario, direction, treat) %>%
  summarise(mean_est = mean(mean), .groups = "drop") %>%
  pivot_wider(
    names_from   = treat,
    values_from  = mean_est,
    names_prefix = "treat_"
  ) %>%
  mutate(effect = treat_0.25 - treat_0)

write.csv(
  gf_contrast,
  file.path(perf_dir, "tables", "gformula_interventional_effects.csv"),
  row.names = FALSE
)

###############################################################################
# Parametric g-formula — Final Performance Evaluation (Revised Visuals)
#
# Visual principles applied consistently:
#   - Causal direction shown via facet_wrap(), not colour competition
#   - Colours retained but used as fill only
#   - Reduced overplotting and clutter
#   - Clear within-panel comparisons
#   - Appendix-safe but publication-quality
###############################################################################

## ---------------------------------------------------------------------------
## Publication theme (shared across all figures)
## ---------------------------------------------------------------------------

theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(size = 0.3, colour = "grey85"),
    
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10),
    
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    
    legend.position = "none",
    
    plot.title = element_text(
      size = 12, face = "bold", hjust = 0
    )
  )

## ---------------------------------------------------------------------------
## Figure 1: Counterfactual means under static interventions
## ---------------------------------------------------------------------------
## Interpretation:
##   Within each panel (direction), compare treatment levels.
##   Across panels, compare causal direction.

p_means <- ggplot(
  gf_summary,
  aes(
    x = treat,
    y = mean_est,
    ymin = lcl,
    ymax = ucl,
    fill = direction
  )
) +
  geom_errorbar(
    width = 0.08,
    linewidth = 0.6
  ) +
  geom_point(
    size = 2.8,
    shape = 21,
    colour = "black"
  ) +
  facet_grid(direction ~ scenario, scales = "free_y") +
  scale_fill_manual(values = morris_fill) +
  labs(
    x = "Static intervention level",
    y = "Counterfactual mean outcome",
    title = "Counterfactual mean outcomes under static interventions"
  ) +
  theme_pub

ggsave(
  filename = file.path(perf_dir, "figures", "gformula_counterfactual_means.png"),
  plot = p_means,
  width = 9,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

## ---------------------------------------------------------------------------
## Figure 2: Interventional effects (main causal estimand)
## ---------------------------------------------------------------------------
## Contrast definition:
##   mean(0.25) − mean(0)
##
## Bars shown only for effect magnitude; direction separated by facets.

p_effect <- ggplot(
  gf_contrast,
  aes(
    x = scenario,
    y = effect,
    fill = direction
  )
) +
  geom_col(width = 0.6) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey40"
  ) +
  facet_wrap(~ direction, nrow = 1) +
  scale_fill_manual(values = morris_fill) +
  labs(
    x = "Scenario",
    y = "Mean difference (0.25 − 0)",
    title = "Interventional effects estimated via the parametric g-formula"
  ) +
  theme_pub

ggsave(
  filename = file.path(perf_dir, "figures", "gformula_interventional_effects.png"),
  plot = p_effect,
  width = 8,
  height = 4,
  dpi = 300,
  bg = "white"
)

## ---------------------------------------------------------------------------
## Figure 3: Monte Carlo variability of g-formula estimates (Appendix)
## ---------------------------------------------------------------------------
## Design rationale:
##   - Distributional shape shown via violin
##   - Central tendency via thin boxplot
##   - Outliers suppressed (Monte Carlo tails are not inferential)
##   - Direction separated by facets for cognitive clarity

p_mc <- ggplot(
  gf,
  aes(
    x = treat,
    y = mean,
    fill = direction
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.8,
    linewidth = 0.3
  ) +
  geom_boxplot(
    width = 0.12,
    linewidth = 0.4,
    outlier.shape = NA
  ) +
  facet_grid(direction ~ scenario, scales = "free_y") +
  scale_fill_manual(values = morris_fill) +
  labs(
    x = "Static intervention level",
    y = "Counterfactual mean (replicate level)",
    title = "Monte Carlo variability of g-formula estimates"
  ) +
  theme_pub

ggsave(
  filename = file.path(perf_dir, "figures", "gformula_replicate_distributions.png"),
  plot = p_mc,
  width = 9,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

## ---------------------------------------------------------------------------
## End of revised visualisation section
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## 11. Final message
## ---------------------------------------------------------------------------

cat(
  "\nParametric g-formula performance evaluation completed.\n",
  "Results saved under:\n  ", perf_dir, "\n",
  sep = ""
)
