# =====================================================================
# DESCRIPTIVE STATISTICS + POOLED DATASET CREATION FOR SAVED DGM DATA
# ---------------------------------------------------------------------
# Purpose:
#   Reads all saved .rds files created by the simulation pipeline,
#   produces descriptive statistics for each time point, and creates
#   pooled datasets that can be used later for g-formula analyses.
#
# Expected input files:
#   data/dgm1_rep1.rds, data/dgm1_rep2.rds, ..., data/dgm3_rep2000.rds
#
# Main outputs:
#   meta/descriptive/descriptive_stats_by_rep_time.csv
#   meta/descriptive/descriptive_stats_by_rep_time.rds
#
#   meta/descriptive/descriptive_stats_across_reps.csv
#   meta/descriptive/descriptive_stats_across_reps.rds
#
#   meta/descriptive/correlations_by_rep_time.csv
#   meta/descriptive/correlations_by_rep_time.rds
#
#   meta/descriptive/correlations_across_reps.csv
#   meta/descriptive/correlations_across_reps.rds
#
#   meta/descriptive/pooled_person_time_scenario1.csv.gz
#   meta/descriptive/pooled_person_time_scenario2.csv.gz
#   meta/descriptive/pooled_person_time_scenario3.csv.gz
#
#   meta/descriptive/pooled_baseline_scenario1.csv.gz
#   meta/descriptive/pooled_baseline_scenario2.csv.gz
#   meta/descriptive/pooled_baseline_scenario3.csv.gz
#
#   meta/descriptive/output_documentation.txt
#   meta/descriptive/data_dictionary.csv
#
# Notes:
#   - Replicate-level descriptives are exact for each saved dataset.
#   - Across-replicate descriptives are summaries of replicate-level
#     statistics, not pooled subject-level summaries across all files.
#   - The pooled person-time datasets stack all replicates within each
#     scenario into one long dataset and add a unique pooled_id.
# =====================================================================

# =====================================================================
# 0. WORKING DIRECTORY + LIBRARIES
# =====================================================================
setwd("/data/cephfs-1/home/users/tato10_c/work/causal_inference_panel_data")

library(data.table)
library(tidyverse)

# =====================================================================
# 1. DIRECTORIES
# =====================================================================
save_dir <- "data"
meta_dir <- "meta"
desc_dir <- file.path(meta_dir, "descriptive")

dir.create(desc_dir, showWarnings = FALSE, recursive = TRUE)

# =====================================================================
# 2. FIND INPUT FILES
# =====================================================================
files <- list.files(
  path = save_dir,
  pattern = "^dgm[0-9]+_rep[0-9]+\\.rds$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No input files found in 'data/'.")
}

# =====================================================================
# 3. HELPER FUNCTIONS
# =====================================================================

# -------------------------------
# Extract scenario and rep_id from filename
# Example filename: dgm2_rep157.rds
# -------------------------------
parse_file_info <- function(path) {
  nm <- basename(path)
  m <- stringr::str_match(nm, "^dgm([0-9]+)_rep([0-9]+)\\.rds$")
  
  tibble(
    file = path,
    scenario = as.integer(m[, 2]),
    rep_id   = as.integer(m[, 3])
  )
}

# -------------------------------
# Descriptive statistics for one numeric vector
# -------------------------------
get_desc <- function(x) {
  tibble(
    n      = sum(!is.na(x)),
    mean   = mean(x, na.rm = TRUE),
    sd     = sd(x, na.rm = TRUE),
    min    = min(x, na.rm = TRUE),
    p25    = unname(quantile(x, 0.25, na.rm = TRUE, type = 7)),
    median = median(x, na.rm = TRUE),
    p75    = unname(quantile(x, 0.75, na.rm = TRUE, type = 7)),
    max    = max(x, na.rm = TRUE)
  )
}

# -------------------------------
# Read one file and return:
#   1) long data for descriptive stats
#   2) descriptive stats by time/variable
#   3) correlations by time
#   4) full dataset with scenario/rep_id/pool identifiers
# -------------------------------
process_one_file <- function(path, scenario_from_name, rep_id_from_name) {
  
  dat <- readRDS(path) %>%
    as_tibble()
  
  # Standardize fields from filename so file naming is the source of truth
  dat <- dat %>%
    mutate(
      scenario = scenario_from_name,
      rep_id   = rep_id_from_name
    ) %>%
    arrange(id, time)
  
  # Unique subject ID after pooling across replicates
  dat <- dat %>%
    mutate(
      pooled_id = paste0("s", scenario, "_r", rep_id, "_i", id)
    )
  
  # Long format for descriptive summaries
  dat_long <- dat %>%
    select(pooled_id, id, rep_id, scenario, time, c, X, Y) %>%
    pivot_longer(
      cols = c(c, X, Y),
      names_to = "variable",
      values_to = "value"
    )
  
  # Descriptive statistics by replicate, time, variable
  desc <- dat_long %>%
    group_by(scenario, rep_id, time, variable) %>%
    group_modify(~ get_desc(.x$value)) %>%
    ungroup()
  
  # Correlations by replicate and time
  cors <- dat %>%
    group_by(scenario, rep_id, time) %>%
    summarise(
      cor_c_X = cor(c, X),
      cor_c_Y = cor(c, Y),
      cor_X_Y = cor(X, Y),
      .groups = "drop"
    )
  
  list(
    desc = desc,
    cors = cors,
    dat  = dat
  )
}

# -------------------------------
# Write documentation text file
# -------------------------------
write_output_documentation <- function(outfile) {
  lines <- c(
    "DESCRIPTIVE OUTPUT DOCUMENTATION",
    "================================",
    "",
    "This folder contains descriptive summaries and pooled datasets",
    "generated from simulation output files stored in the 'data/' folder.",
    "",
    "INPUT DATA STRUCTURE",
    "--------------------",
    "Each input file has one row per person-time observation and includes:",
    "  id        : subject identifier within a replicate",
    "  time      : observed time point",
    "  c         : time-varying confounder",
    "  X         : time-varying exposure",
    "  Y         : time-varying outcome",
    "  scenario  : data-generating scenario",
    "",
    "IMPORTANT OUTPUT FILES",
    "----------------------",
    "",
    "1) descriptive_stats_by_rep_time.csv",
    "   Replicate-level descriptive statistics by scenario, rep_id, time, and variable.",
    "   One row per scenario x rep_id x time x variable.",
    "   This is the exact within-dataset summary for each saved simulation replicate.",
    "",
    "2) descriptive_stats_across_reps.csv",
    "   Summary of replicate-level descriptives across replications.",
    "   One row per scenario x time x variable.",
    "   Example: mean_of_mean is the average of replicate-specific means.",
    "   These are not pooled individual-level quantiles across all files.",
    "",
    "3) correlations_by_rep_time.csv",
    "   Replicate-level correlations among c, X, and Y for each time point.",
    "   One row per scenario x rep_id x time.",
    "",
    "4) correlations_across_reps.csv",
    "   Summary of replicate-level correlations across replications.",
    "   One row per scenario x time.",
    "",
    "5) pooled_person_time_scenario*.csv.gz",
    "   Pooled long-format dataset for each scenario.",
    "   All replicates are stacked together.",
    "   This file is usually the most directly useful starting point for",
    "   longitudinal modeling and later g-formula implementation.",
    "",
    "   Variables in pooled person-time datasets:",
    "     pooled_id  : unique subject identifier across all replicates",
    "     id         : original subject identifier within replicate",
    "     rep_id     : replicate number",
    "     scenario   : scenario number",
    "     time       : observed time point",
    "     c          : time-varying confounder",
    "     X          : time-varying exposure",
    "     Y          : time-varying outcome",
    "",
    "6) pooled_baseline_scenario*.csv.gz",
    "   Baseline-only subset from each pooled scenario dataset.",
    "   Contains time = 0 observations only.",
    "   Useful for baseline tables, initial distribution checks, and any",
    "   models that need baseline covariates only.",
    "",
    "INTERPRETATION NOTES",
    "--------------------",
    "Replicate-level descriptives are exact summaries within one simulation run.",
    "Across-replicate summaries quantify how those replicate-level summaries vary",
    "across simulations.",
    "",
    "Because all replicates were generated with the same sample size (N = 1000),",
    "the average of replicate-specific means corresponds to the pooled mean across",
    "replicates. However, the average of replicate-specific quantiles is not the",
    "same thing as pooled individual-level quantiles from all records combined.",
    "",
    "RECOMMENDED USE FOR G-FORMULA PREPARATION",
    "-----------------------------------------",
    "For g-formula work, start from the pooled_person_time_scenario*.csv.gz files.",
    "These retain the full person-time structure and replicate indicator.",
    "Depending on your modeling plan, you may later create lagged variables such as:",
    "  c_lag1, X_lag1, Y_lag1",
    "within each pooled_id trajectory."
  )
  
  writeLines(lines, outfile)
}

# -------------------------------
# Data dictionary
# -------------------------------
make_data_dictionary <- function() {
  tibble(
    file_group = c(
      rep("descriptive_stats_by_rep_time", 12),
      rep("descriptive_stats_across_reps", 14),
      rep("correlations_by_rep_time", 6),
      rep("correlations_across_reps", 8),
      rep("pooled_person_time_scenario*", 8),
      rep("pooled_baseline_scenario*", 8)
    ),
    variable = c(
      "scenario", "rep_id", "time", "variable", "n", "mean", "sd", "min", "p25", "median", "p75", "max",
      "scenario", "time", "variable", "n_reps", "mean_of_mean", "sd_of_mean", "mean_of_sd", "sd_of_sd",
      "mean_of_min", "mean_of_p25", "mean_of_median", "mean_of_p75", "mean_of_max", "min_mean", "max_mean", "min_median",
      "scenario", "rep_id", "time", "cor_c_X", "cor_c_Y", "cor_X_Y",
      "scenario", "time", "n_reps", "mean_cor_c_X", "sd_cor_c_X", "mean_cor_c_Y", "sd_cor_c_Y", "mean_cor_X_Y",
      "pooled_id", "id", "rep_id", "scenario", "time", "c", "X", "Y",
      "pooled_id", "id", "rep_id", "scenario", "time", "c", "X", "Y"
    ),
    description = c(
      "Simulation scenario number",
      "Replication index",
      "Observed time point",
      "Variable name: c, X, or Y",
      "Number of non-missing observations",
      "Arithmetic mean",
      "Standard deviation",
      "Minimum",
      "25th percentile",
      "Median",
      "75th percentile",
      "Maximum",
      
      "Simulation scenario number",
      "Observed time point",
      "Variable name: c, X, or Y",
      "Number of replications included",
      "Mean of replicate-specific means",
      "Standard deviation of replicate-specific means",
      "Mean of replicate-specific standard deviations",
      "Standard deviation of replicate-specific standard deviations",
      "Mean of replicate-specific minima",
      "Mean of replicate-specific 25th percentiles",
      "Mean of replicate-specific medians",
      "Mean of replicate-specific 75th percentiles",
      "Mean of replicate-specific maxima",
      "Minimum among replicate-specific means",
      "Maximum among replicate-specific means",
      "Minimum among replicate-specific medians",
      
      "Simulation scenario number",
      "Replication index",
      "Observed time point",
      "Correlation between c and X",
      "Correlation between c and Y",
      "Correlation between X and Y",
      
      "Simulation scenario number",
      "Observed time point",
      "Number of replications included",
      "Mean of replicate-specific cor(c, X)",
      "Standard deviation of replicate-specific cor(c, X)",
      "Mean of replicate-specific cor(c, Y)",
      "Standard deviation of replicate-specific cor(c, Y)",
      "Mean of replicate-specific cor(X, Y)",
      
      "Unique subject ID across all replicates within scenario",
      "Original subject ID within one replicate",
      "Replication index",
      "Simulation scenario number",
      "Observed time point",
      "Time-varying confounder",
      "Time-varying exposure",
      "Time-varying outcome",
      
      "Unique subject ID across all replicates within scenario",
      "Original subject ID within one replicate",
      "Replication index",
      "Simulation scenario number",
      "Observed time point, always 0 in baseline file",
      "Baseline value of time-varying confounder",
      "Baseline value of time-varying exposure",
      "Baseline value of time-varying outcome"
    )
  )
}

# =====================================================================
# 4. READ FILE INFO
# =====================================================================
file_info <- bind_rows(lapply(files, parse_file_info)) %>%
  arrange(scenario, rep_id)

# =====================================================================
# 5. PROCESS ALL FILES
# =====================================================================
results <- vector("list", nrow(file_info))

for (i in seq_len(nrow(file_info))) {
  fi <- file_info[i, ]
  
  if (i %% 100 == 0) {
    message("Processing ", i, " / ", nrow(file_info))
  }
  
  results[[i]] <- process_one_file(
    path = fi$file,
    scenario_from_name = fi$scenario,
    rep_id_from_name = fi$rep_id
  )
}

# Combine all outputs
desc_by_rep_time <- bind_rows(lapply(results, `[[`, "desc")) %>%
  arrange(scenario, rep_id, variable, time)

cors_by_rep_time <- bind_rows(lapply(results, `[[`, "cors")) %>%
  arrange(scenario, rep_id, time)

all_dat <- bind_rows(lapply(results, `[[`, "dat")) %>%
  arrange(scenario, rep_id, id, time)

# =====================================================================
# 6. ACROSS-REPLICATE SUMMARIES
# =====================================================================
desc_across_reps <- desc_by_rep_time %>%
  group_by(scenario, time, variable) %>%
  summarise(
    n_reps         = n_distinct(rep_id),
    mean_of_mean   = mean(mean),
    sd_of_mean     = sd(mean),
    mean_of_sd     = mean(sd),
    sd_of_sd       = sd(sd),
    mean_of_min    = mean(min),
    mean_of_p25    = mean(p25),
    mean_of_median = mean(median),
    mean_of_p75    = mean(p75),
    mean_of_max    = mean(max),
    min_mean       = min(mean),
    max_mean       = max(mean),
    min_median     = min(median),
    max_median     = max(median),
    .groups = "drop"
  ) %>%
  arrange(scenario, variable, time)

cors_across_reps <- cors_by_rep_time %>%
  group_by(scenario, time) %>%
  summarise(
    n_reps       = n_distinct(rep_id),
    mean_cor_c_X = mean(cor_c_X),
    sd_cor_c_X   = sd(cor_c_X),
    mean_cor_c_Y = mean(cor_c_Y),
    sd_cor_c_Y   = sd(cor_c_Y),
    mean_cor_X_Y = mean(cor_X_Y),
    sd_cor_X_Y   = sd(cor_X_Y),
    .groups = "drop"
  ) %>%
  arrange(scenario, time)

# =====================================================================
# 7. POOLED DATASETS PER SCENARIO
# =====================================================================
pooled_s1 <- all_dat %>%
  filter(scenario == 1) %>%
  select(pooled_id, id, rep_id, scenario, time, c, X, Y)

pooled_s2 <- all_dat %>%
  filter(scenario == 2) %>%
  select(pooled_id, id, rep_id, scenario, time, c, X, Y)

pooled_s3 <- all_dat %>%
  filter(scenario == 3) %>%
  select(pooled_id, id, rep_id, scenario, time, c, X, Y)

baseline_s1 <- pooled_s1 %>%
  filter(time == 0)

baseline_s2 <- pooled_s2 %>%
  filter(time == 0)

baseline_s3 <- pooled_s3 %>%
  filter(time == 0)

# =====================================================================
# 8. SAVE OUTPUTS
# =====================================================================

# -------------------------------
# Replicate-level descriptive statistics
# -------------------------------
write_csv(
  desc_by_rep_time,
  file.path(desc_dir, "descriptive_stats_by_rep_time.csv")
)
saveRDS(
  desc_by_rep_time,
  file.path(desc_dir, "descriptive_stats_by_rep_time.rds")
)

# -------------------------------
# Across-replicate summaries
# -------------------------------
write_csv(
  desc_across_reps,
  file.path(desc_dir, "descriptive_stats_across_reps.csv")
)
saveRDS(
  desc_across_reps,
  file.path(desc_dir, "descriptive_stats_across_reps.rds")
)

# -------------------------------
# Correlations
# -------------------------------
write_csv(
  cors_by_rep_time,
  file.path(desc_dir, "correlations_by_rep_time.csv")
)
saveRDS(
  cors_by_rep_time,
  file.path(desc_dir, "correlations_by_rep_time.rds")
)

write_csv(
  cors_across_reps,
  file.path(desc_dir, "correlations_across_reps.csv")
)
saveRDS(
  cors_across_reps,
  file.path(desc_dir, "correlations_across_reps.rds")
)

# -------------------------------
# Pooled person-time datasets
# Gzipped CSVs keep file size more reasonable on HPC
# -------------------------------
write_csv(
  pooled_s1,
  file.path(desc_dir, "pooled_person_time_scenario1.csv.gz")
)
write_csv(
  pooled_s2,
  file.path(desc_dir, "pooled_person_time_scenario2.csv.gz")
)
write_csv(
  pooled_s3,
  file.path(desc_dir, "pooled_person_time_scenario3.csv.gz")
)

saveRDS(
  pooled_s1,
  file.path(desc_dir, "pooled_person_time_scenario1.rds")
)
saveRDS(
  pooled_s2,
  file.path(desc_dir, "pooled_person_time_scenario2.rds")
)
saveRDS(
  pooled_s3,
  file.path(desc_dir, "pooled_person_time_scenario3.rds")
)

# -------------------------------
# Pooled baseline datasets
# -------------------------------
write_csv(
  baseline_s1,
  file.path(desc_dir, "pooled_baseline_scenario1.csv.gz")
)
write_csv(
  baseline_s2,
  file.path(desc_dir, "pooled_baseline_scenario2.csv.gz")
)
write_csv(
  baseline_s3,
  file.path(desc_dir, "pooled_baseline_scenario3.csv.gz")
)

saveRDS(
  baseline_s1,
  file.path(desc_dir, "pooled_baseline_scenario1.rds")
)
saveRDS(
  baseline_s2,
  file.path(desc_dir, "pooled_baseline_scenario2.rds")
)
saveRDS(
  baseline_s3,
  file.path(desc_dir, "pooled_baseline_scenario3.rds")
)

# -------------------------------
# Documentation
# -------------------------------
write_output_documentation(
  file.path(desc_dir, "output_documentation.txt")
)

data_dictionary <- make_data_dictionary()
write_csv(
  data_dictionary,
  file.path(desc_dir, "data_dictionary.csv")
)

# -------------------------------
# Small machine-readable summary
# -------------------------------
summary_overview <- tibble(
  quantity = c(
    "n_input_files",
    "n_rows_all_dat",
    "n_unique_subjects_all_dat",
    "n_rows_pooled_s1",
    "n_rows_pooled_s2",
    "n_rows_pooled_s3",
    "n_rows_baseline_s1",
    "n_rows_baseline_s2",
    "n_rows_baseline_s3"
  ),
  value = c(
    length(files),
    nrow(all_dat),
    n_distinct(all_dat$pooled_id),
    nrow(pooled_s1),
    nrow(pooled_s2),
    nrow(pooled_s3),
    nrow(baseline_s1),
    nrow(baseline_s2),
    nrow(baseline_s3)
  )
)

write_csv(
  summary_overview,
  file.path(desc_dir, "summary_overview.csv")
)

message("All descriptive outputs and pooled datasets saved to: ", desc_dir)
