###############################################################################
# SEM MONTE CARLO: RI-CLPM WITH TIME-VARYING CONFOUNDER (95% CI)
# FULL PIPELINE IMPLEMENTATION (fixed beta extraction + more stable variances)
#
# Fixes vs your previous run:
#   1) beta_xy / beta_yx appear 3× each (t1->t2, t2->t3, t3->t4) = 6 rows.
#      We collapse by label and take slice(1) to match your old pipeline style.
#   2) Adds equality constraints on within-person latent variances to reduce
#      Heywood cases (negative lv variances).
#
# Outputs:
#   - fits:        methods/sem_riclpm/fits/sem_riclpm_fit_s{S}_rep{R}.rds
#   - results all: methods/sem_riclpm/sem_riclpm_results_all.{rds,csv}
#   - performance: methods/sem_riclpm/performance/sem_riclpm_performance_by_scenario.{rds,csv}
###############################################################################

setwd("/data/cephfs-1/home/users/tato10_c/work/causal_inference_panel_data")

library(lavaan)
library(tidyverse)
library(data.table)
library(readr)

# -----------------------------------------------------------------------------
# Directories (pipeline style)
# -----------------------------------------------------------------------------
data_dir    <- file.path(getwd(), "data")
method_dir  <- file.path(getwd(), "methods", "sem_riclpm")
fits_dir    <- file.path(method_dir, "fits")
perf_dir    <- file.path(method_dir, "performance")

dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fits_dir,   showWarnings = FALSE)
dir.create(perf_dir,   showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
nsim     <- 2000L
ci_level <- 0.95

# -----------------------------------------------------------------------------
# True estimands (same as your other scripts)
# -----------------------------------------------------------------------------
true_estimands <- function(scenario) {
  tibble(
    ACE_XY = ifelse(scenario >= 2, 0.8, 0),
    ACE_YX = ifelse(scenario == 3, 0.5, 0)
  )
}

# -----------------------------------------------------------------------------
# Load dataset
# -----------------------------------------------------------------------------
get_dataset <- function(scenario, rep_id) {
  f <- file.path(data_dir, sprintf("dgm%i_rep%i.rds", scenario, rep_id))
  if (!file.exists(f)) stop("Dataset missing: ", f)
  readRDS(f)
}

# -----------------------------------------------------------------------------
# RI-CLPM model (T=4) + time-varying confounder c_t
# - Cross-lags on within-person deviations (xw/yw)
# - Confounder AR estimated (rho_c)
# - Added equality constraints on within-person latent variances (v_xw, v_yw)
#   to reduce Heywood cases.
# -----------------------------------------------------------------------------
sem_model_riclpm <- "
  # (A) Random intercepts (between-person)
  RI_X =~ 1*X_1 + 1*X_2 + 1*X_3 + 1*X_4
  RI_Y =~ 1*Y_1 + 1*Y_2 + 1*Y_3 + 1*Y_4
  RI_X ~~ RI_Y
  RI_X ~~ RI_X
  RI_Y ~~ RI_Y

  # (B) Within-person deviations
  xw1 =~ 1*X_1
  xw2 =~ 1*X_2
  xw3 =~ 1*X_3
  xw4 =~ 1*X_4

  yw1 =~ 1*Y_1
  yw2 =~ 1*Y_2
  yw3 =~ 1*Y_3
  yw4 =~ 1*Y_4

  # (C) Exact decomposition: X_t = RI_X + xw_t, Y_t = RI_Y + yw_t
  X_1 ~~ 0*X_1
  X_2 ~~ 0*X_2
  X_3 ~~ 0*X_3
  X_4 ~~ 0*X_4

  Y_1 ~~ 0*Y_1
  Y_2 ~~ 0*Y_2
  Y_3 ~~ 0*Y_3
  Y_4 ~~ 0*Y_4

  # (D) Orthogonality (RI uncorrelated with within deviations)
  RI_X ~~ 0*xw1 + 0*xw2 + 0*xw3 + 0*xw4
  RI_Y ~~ 0*yw1 + 0*yw2 + 0*yw3 + 0*yw4
  RI_X ~~ 0*yw1 + 0*yw2 + 0*yw3 + 0*yw4
  RI_Y ~~ 0*xw1 + 0*xw2 + 0*xw3 + 0*xw4

  # (E) Within-wave covariances (within-person synchronous association)
  xw1 ~~ yw1
  xw2 ~~ yw2
  xw3 ~~ yw3
  xw4 ~~ yw4

  # (F) CHANGED: constrain within-person latent variances equal across waves
  # WHY: reduces Heywood cases (negative lv variances) in RI-CLPM fits.
  xw1 ~~ v_xw*xw1
  xw2 ~~ v_xw*xw2
  xw3 ~~ v_xw*xw3
  xw4 ~~ v_xw*xw4

  yw1 ~~ v_yw*yw1
  yw2 ~~ v_yw*yw2
  yw3 ~~ v_yw*yw3
  yw4 ~~ v_yw*yw4

  # (G) Confounder AR(1) (estimated, NOT forced)
  c_2 ~ rho_c*c_1
  c_3 ~ rho_c*c_2
  c_4 ~ rho_c*c_3

  # (H) Baseline within deviations depend on baseline confounder
  xw1 ~ gamma_cx0*c_1
  yw1 ~ gamma_cy0*c_1

  # (I) Within-person dynamics + time-varying confounding
  xw2 ~ phi_x*xw1 + beta_yx*yw1 + gamma_cx*c_1
  yw2 ~ phi_y*yw1 + beta_xy*xw1 + gamma_cy*c_1

  xw3 ~ phi_x*xw2 + beta_yx*yw2 + gamma_cx*c_2
  yw3 ~ phi_y*yw2 + beta_xy*xw2 + gamma_cy*c_2

  xw4 ~ phi_x*xw3 + beta_yx*yw3 + gamma_cx*c_3
  yw4 ~ phi_y*yw3 + beta_xy*xw3 + gamma_cy*c_3

  # (J) Identify within-person deviation means as 0
  xw1 ~ 0*1
  xw2 ~ 0*1
  xw3 ~ 0*1
  xw4 ~ 0*1
  yw1 ~ 0*1
  yw2 ~ 0*1
  yw3 ~ 0*1
  yw4 ~ 0*1
"

# -----------------------------------------------------------------------------
# Morris-style MC metrics (same as your master script)
# -----------------------------------------------------------------------------
mc_se <- function(estimates, theta, ci_low, ci_high, se_mod) {
  nsim <- length(estimates)
  
  bias <- mean(estimates) - theta
  bias_mcse <- sqrt(var(estimates) / nsim)
  
  empSE <- sd(estimates)
  empSE_mcse <- empSE * sqrt(1 / (2 * (nsim - 1)))
  
  MSE <- mean((estimates - theta)^2)
  MSE_mcse <- sqrt(sum(((estimates - theta)^2 - MSE)^2) / (nsim * (nsim - 1)))
  
  cover <- mean(ci_low <= theta & ci_high >= theta)
  cover_mcse <- sqrt(cover * (1 - cover) / nsim)
  
  s2 <- se_mod^2
  m_s2 <- mean(s2, na.rm = TRUE)
  ModSE <- sqrt(m_s2)
  
  var_m <- var(s2, na.rm = TRUE) / nsim
  ModSE_mcse <- if (m_s2 > 0) sqrt(var_m / (4 * m_s2)) else NA_real_
  
  theta_hat <- mean(estimates)
  cover_be  <- mean(ci_low <= theta_hat & ci_high >= theta_hat)
  cover_be_mcse <- sqrt(cover_be * (1 - cover_be) / nsim)
  
  rel_ModSE <- ModSE / empSE - 1
  rel_ModSE_mcse <- sqrt(
    (ModSE_mcse / empSE)^2 +
      (ModSE * empSE_mcse / empSE^2)^2
  )
  
  data.frame(
    bias, bias_mcse,
    empSE, empSE_mcse,
    MSE, MSE_mcse,
    cover, cover_mcse,
    ModSE, ModSE_mcse,
    cover_be, cover_be_mcse,
    rel_ModSE, rel_ModSE_mcse
  )
}

# -----------------------------------------------------------------------------
# Fit one replicate (scenario, rep)
# -----------------------------------------------------------------------------
fit_sem_riclpm <- function(scenario, rep_id) {
  
  dat <- get_dataset(scenario, rep_id)
  setDT(dat)
  
  df_wide <- dat %>%
    select(id, time, X, Y, c) %>%
    pivot_wider(
      names_from  = time,
      values_from = c(X, Y, c),
      names_sep   = "_"
    )
  
  fit <- sem(
    model = sem_model_riclpm,
    data  = df_wide,
    fixed.x = FALSE,
    missing = "ML",
    meanstructure = TRUE
  )
  
  saveRDS(
    fit,
    file.path(fits_dir, sprintf("sem_riclpm_fit_s%i_rep%i.rds", scenario, rep_id))
  )
  
  est_table <- parameterEstimates(
    fit,
    standardized = FALSE,
    ci = TRUE,
    level = ci_level
  )
  
  # CHANGED: collapse repeated labeled rows (beta_* appears at 3 transitions)
  # WHY: lavaan reports separate rows for each equation even if labels match.
  # We take the FIRST occurrence (t1->t2), matching your earlier SEM pipeline style.
  est_xy_yx <- est_table %>%
    filter(label %in% c("beta_xy", "beta_yx")) %>%
    group_by(label) %>%
    slice(1) %>%          # CHANGED
    ungroup() %>%
    select(label, est, se, ci.lower, ci.upper)
  
  if (nrow(est_xy_yx) != 2) {
    stop("Expected 2 rows after collapsing repeated labels; got ", nrow(est_xy_yx),
         " (scenario ", scenario, ", rep ", rep_id, ")")
  }
  
  truth <- true_estimands(scenario)
  
  tibble(
    scenario = scenario,
    rep_id   = rep_id,
    
    beta_xy_hat  = est_xy_yx$est[est_xy_yx$label == "beta_xy"],
    beta_xy_se   = est_xy_yx$se[est_xy_yx$label == "beta_xy"],
    beta_xy_low  = est_xy_yx$ci.lower[est_xy_yx$label == "beta_xy"],
    beta_xy_high = est_xy_yx$ci.upper[est_xy_yx$label == "beta_xy"],
    beta_xy_true = truth$ACE_XY,
    
    beta_yx_hat  = est_xy_yx$est[est_xy_yx$label == "beta_yx"],
    beta_yx_se   = est_xy_yx$se[est_xy_yx$label == "beta_yx"],
    beta_yx_low  = est_xy_yx$ci.lower[est_xy_yx$label == "beta_yx"],
    beta_yx_high = est_xy_yx$ci.upper[est_xy_yx$label == "beta_yx"],
    beta_yx_true = truth$ACE_YX
  )
}

# -----------------------------------------------------------------------------
# Run all scenarios
# -----------------------------------------------------------------------------
set.seed(20251026)

res_list <- vector("list", 3)

for (s in 1:3) {
  message("Running RI-CLPM — Scenario ", s)
  out <- vector("list", nsim)
  
  for (r in 1:nsim) {
    if (r %% 50 == 0) message("  Replicate ", r, "/", nsim)
    out[[r]] <- fit_sem_riclpm(s, r)
  }
  
  res_list[[s]] <- bind_rows(out)
}

sem_riclpm_results <- bind_rows(res_list)

# -----------------------------------------------------------------------------
# Save replicate-level results
# -----------------------------------------------------------------------------
saveRDS(sem_riclpm_results, file.path(method_dir, "sem_riclpm_results_all.rds"))
readr::write_csv(sem_riclpm_results, file.path(method_dir, "sem_riclpm_results_all.csv"))

# -----------------------------------------------------------------------------
# Performance (scenario × effect) in Morris-style
# -----------------------------------------------------------------------------
sem_long <- bind_rows(
  sem_riclpm_results %>%
    transmute(
      method = "SEM_RI_CLPM",
      scenario, rep_id,
      effect = "XY",
      est = beta_xy_hat,
      ci_low = beta_xy_low,
      ci_high = beta_xy_high,
      theta = beta_xy_true,
      se_mod = beta_xy_se
    ),
  sem_riclpm_results %>%
    transmute(
      method = "SEM_RI_CLPM",
      scenario, rep_id,
      effect = "YX",
      est = beta_yx_hat,
      ci_low = beta_yx_low,
      ci_high = beta_yx_high,
      theta = beta_yx_true,
      se_mod = beta_yx_se
    )
)

sem_dt <- as.data.table(sem_long)

sem_riclpm_perf <- sem_dt[, mc_se(est, unique(theta), ci_low, ci_high, se_mod),
                          by = .(scenario, effect)]

sem_riclpm_perf[, `:=`(
  RMSE = sqrt(MSE),
  RMSE_mcse = ifelse(MSE > 0, 0.5 * MSE_mcse / sqrt(MSE), NA_real_)
)]

sem_riclpm_perf <- as_tibble(sem_riclpm_perf) %>%
  mutate(
    scenario_label = recode(
      as.character(scenario),
      "1" = "S1: No causal effect",
      "2" = "S2: X → Y",
      "3" = "S3: X ↔ Y"
    )
  ) %>%
  select(
    scenario, scenario_label, effect,
    bias, bias_mcse,
    empSE, empSE_mcse,
    MSE, MSE_mcse, RMSE, RMSE_mcse,
    cover, cover_mcse,
    ModSE, ModSE_mcse,
    cover_be, cover_be_mcse,
    rel_ModSE, rel_ModSE_mcse
  )

saveRDS(sem_riclpm_perf, file.path(perf_dir, "sem_riclpm_performance_by_scenario.rds"))
readr::write_csv(sem_riclpm_perf, file.path(perf_dir, "sem_riclpm_performance_by_scenario.csv"))

message("FINISHED RI-CLPM. Outputs in: ", method_dir)
