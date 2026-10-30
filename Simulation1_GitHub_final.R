# ======================================================================
# SIMULATION 1
# Monte Carlo study of miscalibration, Net Benefit, Weighted Net Benefit,
# prior sensitivity, and comparison with local frequentist calibration error
# ======================================================================

# ----------------------------------------------------------------------
# Packages and reproducibility
# ----------------------------------------------------------------------
library(tidyverse)
library(pROC)
library(furrr)
library(progressr)
library(ggpubr)

set.seed(123)

plan(multisession, workers = max(1, parallel::detectCores() - 1))
handlers(global = TRUE)
handlers("txtprogressbar")

# ----------------------------------------------------------------------
# 1. Simulation settings
# ----------------------------------------------------------------------
target_prev <- 0.20
target_auc <- 0.75

n_values <- c(500L, 1000L, 10000L)
n_sim <- 1000L

thresholds <- c(0.10, 0.20, 0.30, 0.40, 0.50)
reported_thresholds <- c(0.10, 0.30, 0.50)
bin_width <- 0.10

# Prior mean is fixed at the target/global prevalence (0.20).
# Prior concentration m = alpha + beta.
# Complete sensitivity analysis: m = 2, 4, ..., 50.
prior_m_values <- seq(2, 50, by = 2)
reference_prior_m <- 10
paper_prior_m_values <- c(2, 10, 50)

scenarios <- c(
  "Underestimation",
  "Overestimation",
  "Under-Over",
  "Over-Under"
)

output_dir <- "simulation1_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------
# 2. Safe summary helpers
# ----------------------------------------------------------------------
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(quantile(x, probs = prob, names = FALSE, type = 7))
}

# ----------------------------------------------------------------------
# 3. Data-generating mechanism
# ----------------------------------------------------------------------
generate_lp <- function(n, target_auc, target_prev) {
  sigma_lp <- sqrt(2) * qnorm(target_auc)
  lp_raw <- rnorm(n, mean = 0, sd = sigma_lp)

  intercept <- uniroot(
    function(a) mean(plogis(a + lp_raw)) - target_prev,
    interval = c(-15, 15)
  )$root

  prob <- plogis(intercept + lp_raw)
  y <- rbinom(n = n, size = 1, prob = prob)

  list(prob = prob, y = y)
}

# ----------------------------------------------------------------------
# 4. Miscalibration transformations
# ----------------------------------------------------------------------
shift_val <- 0.5
center_lp <- qlogis(0.3)
tilt_val <- 1.4

miscal_transform <- function(p, scenario) {
  switch(
    scenario,
    "Underestimation" = plogis(-shift_val + qlogis(p)),
    "Overestimation" = plogis(+shift_val + qlogis(p)),
    "Under-Over" = plogis(
      center_lp + tilt_val * (qlogis(p) - center_lp)
    ),
    "Over-Under" = plogis(
      center_lp + (1 / tilt_val) * (qlogis(p) - center_lp)
    ),
    stop("Unknown scenario: ", scenario)
  )
}

# ----------------------------------------------------------------------
# 5. Net Benefit
# ----------------------------------------------------------------------
net_benefit <- function(y, p, t) {
  tp <- mean(p >= t & y == 1)
  fp <- mean(p >= t & y == 0)

  tp - fp * t / (1 - t)
}

# ----------------------------------------------------------------------
# 6. Local posterior quantities and calibration errors
# ----------------------------------------------------------------------
compute_local_metrics <- function(y, p, t, prior_m) {
  idx <- which(
    p >= (t - bin_width / 2) &
      p < (t + bin_width / 2)
  )

  p0 <- target_prev
  prior_alpha <- prior_m * p0
  prior_beta <- prior_m * (1 - p0)

  # Empty local bins are treated as non-estimable, not as zero error.
  if (length(idx) == 0L) {
    return(
      tibble(
        prior_m = prior_m,
        prior_alpha = prior_alpha,
        prior_beta = prior_beta,
        n_bin = 0L,
        k = 0L,
        pbar = NA_real_,
        empirical_rate = NA_real_,
        post_mean = NA_real_,
        prior_weight = NA_real_,
        signed_prior_pull = NA_real_,
        abs_prior_pull = NA_real_,
        local_freq_error = NA_real_,
        LCE = NA_real_,
        SDpost = NA_real_
      )
    )
  }

  n_bin <- length(idx)
  k <- sum(y[idx])
  pbar <- mean(p[idx])
  empirical_rate <- k / n_bin

  a_post <- prior_alpha + k
  b_post <- prior_beta + (n_bin - k)

  post_mean <- a_post / (a_post + b_post)
  post_sd <- sqrt(
    a_post * b_post /
      ((a_post + b_post)^2 * (a_post + b_post + 1))
  )

  prior_weight <- prior_m / (n_bin + prior_m)
  signed_prior_pull <- post_mean - empirical_rate
  abs_prior_pull <- abs(signed_prior_pull)

  # Local frequentist calibration error: bin-specific absolute calibration
  # error underlying ECE, evaluated in the same local bin used for LCE.
  local_freq_error <- abs(pbar - empirical_rate)

  tibble(
    prior_m = prior_m,
    prior_alpha = prior_alpha,
    prior_beta = prior_beta,
    n_bin = n_bin,
    k = k,
    pbar = pbar,
    empirical_rate = empirical_rate,
    post_mean = post_mean,
    prior_weight = prior_weight,
    signed_prior_pull = signed_prior_pull,
    abs_prior_pull = abs_prior_pull,
    local_freq_error = local_freq_error,
    LCE = abs(pbar - post_mean),
    SDpost = post_sd
  )
}

# ----------------------------------------------------------------------
# 7. Simulate one dataset
# ----------------------------------------------------------------------
simulate_dataset <- function(n, scenario) {
  dat <- generate_lp(
    n = n,
    target_auc = target_auc,
    target_prev = target_prev
  )

  y <- dat$y
  p <- miscal_transform(dat$prob, scenario)

  AUC <- if (length(unique(y)) == 2L) {
    as.numeric(
      pROC::auc(
        response = y,
        predictor = p,
        quiet = TRUE
      )
    )
  } else {
    NA_real_
  }

  threshold_tbl <- map_dfr(
    thresholds,
    function(t) {
      NB <- net_benefit(y = y, p = p, t = t)

      map_dfr(
        prior_m_values,
        function(m) {
          compute_local_metrics(
            y = y,
            p = p,
            t = t,
            prior_m = m
          ) %>%
            mutate(
              threshold = t,
              NB = NB,
              WNB = NB / (1 + LCE + SDpost)
            )
        }
      )
    }
  )

  list(auc = AUC, metrics = threshold_tbl)
}

# ----------------------------------------------------------------------
# 8. Full Monte Carlo design
# ----------------------------------------------------------------------
simulation_grid <- tidyr::crossing(
  n = n_values,
  scenario = scenarios,
  sim_id = seq_len(n_sim)
)

# ----------------------------------------------------------------------
# 9. Run simulations
# ----------------------------------------------------------------------
with_progress({
  p_progress <- progressor(steps = nrow(simulation_grid))

  results <- future_pmap_dfr(
    simulation_grid,
    function(n, scenario, sim_id) {
      sim <- simulate_dataset(n = n, scenario = scenario)

      p_progress(
        sprintf(
          "n = %s | %s | replicate %s/%s",
          n,
          scenario,
          sim_id,
          n_sim
        )
      )

      tibble(
        n = n,
        scenario = scenario,
        sim_id = sim_id,
        AUC = sim$auc,
        metrics = list(sim$metrics)
      )
    },
    .options = furrr_options(
      seed = TRUE,
      scheduling = 2
    )
  )
})

plan(sequential)

all_results <- results %>%
  unnest(metrics)

# ----------------------------------------------------------------------
# 10. AUC summary
# ----------------------------------------------------------------------
auc_summary <- results %>%
  group_by(n, scenario) %>%
  summarise(
    AUC_mean = safe_mean(AUC),
    AUC_low = safe_quantile(AUC, 0.025),
    AUC_high = safe_quantile(AUC, 0.975),
    .groups = "drop"
  ) %>%
  arrange(n, factor(scenario, levels = scenarios))

write.csv(
  auc_summary,
  file.path(output_dir, "Simulation1_AUC_summary.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 11. Paired WNB comparison relative to reference m = 10
# ----------------------------------------------------------------------
reference_results <- all_results %>%
  filter(prior_m == reference_prior_m) %>%
  select(
    n,
    scenario,
    sim_id,
    threshold,
    WNB_reference = WNB
  )

sensitivity_results <- all_results %>%
  left_join(
    reference_results,
    by = c("n", "scenario", "sim_id", "threshold")
  ) %>%
  mutate(
    delta_WNB = WNB - WNB_reference,
    abs_delta_WNB = abs(delta_WNB),
    prior_weight_percent = 100 * prior_weight,
    prior_pull_pp = 100 * signed_prior_pull,
    abs_prior_pull_pp = 100 * abs_prior_pull
  )

# ----------------------------------------------------------------------
# 12. Complete prior-sensitivity table: m = 2,4,...,50
# ----------------------------------------------------------------------
prior_sensitivity_table <- sensitivity_results %>%
  group_by(n, scenario, threshold, prior_m) %>%
  summarise(
    mean_n_bin = safe_mean(n_bin),
    empty_bin_percent = 100 * mean(n_bin == 0),
    mean_empirical_rate = safe_mean(empirical_rate),
    mean_posterior_rate = safe_mean(post_mean),
    mean_prior_weight_percent = safe_mean(prior_weight_percent),
    mean_abs_prior_pull_pp = safe_mean(abs_prior_pull_pp),
    mean_LCE = safe_mean(LCE),
    mean_SDpost = safe_mean(SDpost),
    mean_WNB = safe_mean(WNB),
    WNB_low = safe_quantile(WNB, 0.025),
    WNB_high = safe_quantile(WNB, 0.975),
    mean_delta_WNB = safe_mean(delta_WNB),
    mean_abs_delta_WNB = safe_mean(abs_delta_WNB),
    .groups = "drop"
  ) %>%
  arrange(
    n,
    factor(scenario, levels = scenarios),
    threshold,
    prior_m
  )

write.csv(
  prior_sensitivity_table,
  file.path(output_dir, "Simulation1_prior_sensitivity_complete_results.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 13. Selected prior-sensitivity table for manuscript (m = 2,10,50)
# ----------------------------------------------------------------------
prior_sensitivity_paper <- prior_sensitivity_table %>%
  filter(
    prior_m %in% paper_prior_m_values,
    threshold %in% reported_thresholds
  ) %>%
  select(
    n,
    scenario,
    threshold,
    prior_m,
    mean_n_bin,
    empty_bin_percent,
    mean_prior_weight_percent,
    mean_abs_prior_pull_pp,
    mean_LCE,
    mean_SDpost,
    mean_WNB,
    mean_delta_WNB
  ) %>%
  arrange(
    n,
    factor(scenario, levels = scenarios),
    threshold,
    prior_m
  )

write.csv(
  prior_sensitivity_paper,
  file.path(output_dir, "Simulation1_prior_sensitivity_paper_selected.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 14. Main Simulation 1 results: reference prior m = 10
# ----------------------------------------------------------------------
main_results <- sensitivity_results %>%
  filter(prior_m == reference_prior_m) %>%
  group_by(n, scenario, threshold) %>%
  summarise(
    mean_n_bin = safe_mean(n_bin),
    empty_bin_percent = 100 * mean(n_bin == 0),
    NB_mean = safe_mean(NB),
    NB_low = safe_quantile(NB, 0.025),
    NB_high = safe_quantile(NB, 0.975),
    WNB_mean = safe_mean(WNB),
    WNB_low = safe_quantile(WNB, 0.025),
    WNB_high = safe_quantile(WNB, 0.975),
    LCE_mean = safe_mean(LCE),
    SDpost_mean = safe_mean(SDpost),
    .groups = "drop"
  ) %>%
  left_join(auc_summary, by = c("n", "scenario")) %>%
  mutate(
    relative_attenuation_percent = if_else(
      NB_mean > 0,
      100 * (NB_mean - WNB_mean) / NB_mean,
      NA_real_
    )
  ) %>%
  arrange(
    n,
    factor(scenario, levels = scenarios),
    threshold
  )

write.csv(
  main_results,
  file.path(output_dir, "Simulation1_main_results_reference_m10.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 15. LCE vs local frequentist calibration error (Table S5)
# ----------------------------------------------------------------------
# This comparison is restricted to the reference prior m = 10 and uses
# exactly the same local bin for both measures.
lce_vs_freq_complete <- sensitivity_results %>%
  filter(prior_m == reference_prior_m) %>%
  mutate(
    difference_LCE_minus_freq = LCE - local_freq_error
  ) %>%
  group_by(n, scenario, threshold) %>%
  summarise(
    mean_n_bin = safe_mean(n_bin),
    empty_bin_percent = 100 * mean(n_bin == 0),
    mean_predicted_probability = safe_mean(pbar),
    mean_empirical_event_rate = safe_mean(empirical_rate),
    mean_local_freq_error = safe_mean(local_freq_error),
    mean_LCE = safe_mean(LCE),
    mean_difference_LCE_minus_freq = safe_mean(difference_LCE_minus_freq),
    .groups = "drop"
  ) %>%
  arrange(
    n,
    factor(scenario, levels = scenarios),
    threshold
  )

write.csv(
  lce_vs_freq_complete,
  file.path(output_dir, "Simulation1_LCE_vs_local_frequentist_error_complete.csv"),
  row.names = FALSE
)

lce_vs_freq_paper <- lce_vs_freq_complete %>%
  filter(threshold %in% reported_thresholds) %>%
  select(
    n,
    scenario,
    threshold,
    mean_n_bin,
    empty_bin_percent,
    mean_local_freq_error,
    mean_LCE,
    mean_difference_LCE_minus_freq
  )

write.csv(
  lce_vs_freq_paper,
  file.path(output_dir, "Simulation1_LCE_vs_local_frequentist_error_paper.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------
# 16. Figure 1: miscalibration scenarios
# ----------------------------------------------------------------------
calibration_curve_data <- tidyr::crossing(
  true_probability = seq(0.001, 0.999, length.out = 2000),
  scenario = scenarios,
  n = n_values
) %>%
  mutate(
    predicted_probability = map2_dbl(
      true_probability,
      scenario,
      ~ miscal_transform(.x, .y)
    ),
    scenario = factor(scenario, levels = scenarios),
    n = factor(n, levels = n_values, labels = paste0("n = ", n_values))
  ) %>%
  filter(predicted_probability >= 0.08, predicted_probability <= 0.50)

miscalibration_plot <- ggplot(
  calibration_curve_data,
  aes(
    x = predicted_probability,
    y = true_probability,
    color = scenario,
    group = scenario
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_line(linewidth = 1.1) +
  facet_grid(scenario ~ n) +
  scale_color_manual(values = c(
    "Underestimation" = "#1B9E77",
    "Overestimation" = "#D95F02",
    "Under-Over" = "#7570B3",
    "Over-Under" = "#E7298A"
  )) +
  scale_x_continuous(
    breaks = thresholds,
    limits = c(0.08, 0.52)
  ) +
  coord_cartesian(ylim = c(0.04, 0.64)) +
  labs(
    x = "Predicted probability",
    y = "Observed event rate"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12)
  )

ggsave(
  filename = file.path(output_dir, "Figure_miscalibration_scenarios.pdf"),
  plot = miscalibration_plot,
  width = 12,
  height = 10
)

# ----------------------------------------------------------------------
# 17. Figure 2: LCE and posterior SD across m = 2,4,...,50
# ----------------------------------------------------------------------
plot_data <- prior_sensitivity_table %>%
  mutate(
    n = factor(n, levels = n_values, labels = paste0("n = ", n_values)),
    scenario = factor(scenario, levels = scenarios)
  )

lce_plot <- ggplot(
  plot_data,
  aes(
    x = threshold,
    y = mean_LCE,
    color = prior_m,
    group = prior_m
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_grid(scenario ~ n) +
  scale_x_continuous(breaks = thresholds) +
  scale_color_viridis_c(option = "C", end = 0.95) +
  labs(
    x = "Decision threshold",
    y = "Local Calibration Error",
    color = "Prior concentration (m)"
  ) +
  theme_bw(base_size = 12)

sd_plot <- ggplot(
  plot_data,
  aes(
    x = threshold,
    y = mean_SDpost,
    color = prior_m,
    group = prior_m
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_grid(scenario ~ n) +
  scale_x_continuous(breaks = thresholds) +
  scale_color_viridis_c(option = "C", end = 0.95) +
  labs(
    x = "Decision threshold",
    y = "Posterior SD",
    color = "Prior concentration (m)"
  ) +
  theme_bw(base_size = 12)

fig_prior_components <- ggarrange(
  lce_plot,
  sd_plot,
  labels = c("a", "b"),
  ncol = 2,
  nrow = 1,
  align = "h",
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  filename = file.path(output_dir, "Figure_prior_components_m2_to_m50.pdf"),
  plot = fig_prior_components,
  width = 14,
  height = 9
)

# ----------------------------------------------------------------------
# 18. Figure 3: paired change in WNB relative to m = 10
# ----------------------------------------------------------------------
delta_plot_data <- prior_sensitivity_table %>%
  filter(prior_m != reference_prior_m) %>%
  mutate(
    n = factor(n, levels = n_values, labels = paste0("n = ", n_values)),
    scenario = factor(scenario, levels = scenarios)
  )

delta_wnb_plot <- ggplot(
  delta_plot_data,
  aes(
    x = threshold,
    y = mean_delta_WNB,
    color = prior_m,
    group = prior_m
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_grid(scenario ~ n) +
  scale_x_continuous(breaks = thresholds) +
  scale_color_viridis_c(option = "C", end = 0.95) +
  labs(
    x = "Decision threshold",
    y = expression(Delta * " WNB relative to m = 10"),
    color = "Prior concentration (m)"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(output_dir, "Figure_delta_WNB_vs_m10_m2_to_m50.pdf"),
  plot = delta_wnb_plot,
  width = 12,
  height = 10
)

# ----------------------------------------------------------------------
# 19. Console checks
# ----------------------------------------------------------------------
cat("\nSimulation 1 completed.\n")
cat("Output directory:", normalizePath(output_dir), "\n\n")

cat("Main results (reference m = 10):\n")
print(main_results, n = Inf)

cat("\nSelected prior-sensitivity results (m = 2, 10, 50):\n")
print(prior_sensitivity_paper, n = Inf)

cat("\nLCE vs local frequentist calibration error (reported thresholds):\n")
print(lce_vs_freq_paper, n = Inf)

cat("\nMaximum mean absolute change in WNB relative to m = 10:\n")
print(
  prior_sensitivity_table %>%
    filter(prior_m != reference_prior_m) %>%
    arrange(desc(mean_abs_delta_WNB)) %>%
    slice_head(n = 10)
)
