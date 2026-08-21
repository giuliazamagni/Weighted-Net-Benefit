# ======================================================================
# SIMULATION 1
# Monte Carlo study of miscalibration, Net Benefit, and Weighted Net Benefit
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
bin_width <- 0.10

# Prior sensitivity: m = alpha + beta, with prior mean fixed at 0.20.
# m = 2  -> Beta(0.4, 1.6)
# m = 10 -> Beta(2, 8)
# m = 50 -> Beta(10, 40)
prior_m_values <- c(2, 10, 50)
reference_prior_m <- 10

scenarios <- c(
  "Underestimation",
  "Overestimation",
  "Under-Over",
  "Over-Under"
)

# ----------------------------------------------------------------------
# 2. Data-generating mechanism
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
# 3. Miscalibration transformations
# ----------------------------------------------------------------------
shift_val <- 0.5
center_lp <- qlogis(0.3)
tilt_val <- 1.4

miscal_transform <- function(p, scenario) {
  switch(
    scenario,
    "Underestimation" = plogis(-shift_val + qlogis(p)),
    "Overestimation" = plogis(+shift_val + qlogis(p)),
    "Under-Over" = plogis(center_lp + tilt_val * (qlogis(p) - center_lp)),
    "Over-Under" = plogis(center_lp + (1 / tilt_val) * (qlogis(p) - center_lp))
  )
}

# ----------------------------------------------------------------------
# 4. Net Benefit
# ----------------------------------------------------------------------
net_benefit <- function(y, p, t) {
  tp <- mean(p >= t & y == 1)
  fp <- mean(p >= t & y == 0)
  tp - fp * t / (1 - t)
}

# ----------------------------------------------------------------------
# 5. Local Calibration Error + posterior SD
# ----------------------------------------------------------------------
compute_lce_sd <- function(y, p, t, prior_m) {
  idx <- which(p >= t - bin_width / 2 & p < t + bin_width / 2)

  prior_alpha <- prior_m * target_prev
  prior_beta <- prior_m * (1 - target_prev)

  if (length(idx) == 0) {
    return(tibble(
      prior_m = prior_m,
      prior_alpha = prior_alpha,
      prior_beta = prior_beta,
      n_bin = 0L,
      k = 0L,
      pbar = NA_real_,
      post_mean = NA_real_,
      LCE = 0,
      SD = 0
    ))
  }

  n_bin <- length(idx)
  k <- sum(y[idx])
  pbar <- mean(p[idx])

  a_post <- prior_alpha + k
  b_post <- prior_beta + (n_bin - k)
  post_mean <- a_post / (a_post + b_post)
  post_sd <- sqrt(
    a_post * b_post /
      ((a_post + b_post)^2 * (a_post + b_post + 1))
  )

  tibble(
    prior_m = prior_m,
    prior_alpha = prior_alpha,
    prior_beta = prior_beta,
    n_bin = n_bin,
    k = k,
    pbar = pbar,
    post_mean = post_mean,
    LCE = abs(pbar - post_mean),
    SD = post_sd
  )
}

# ----------------------------------------------------------------------
# 6. Simulate one dataset
# ----------------------------------------------------------------------
simulate_dataset <- function(n, scenario) {
  dat <- generate_lp(n = n, target_auc = target_auc, target_prev = target_prev)
  y <- dat$y
  p <- miscal_transform(p = dat$prob, scenario = scenario)

  AUC <- if (length(unique(y)) == 2) {
    as.numeric(pROC::auc(response = y, predictor = p, quiet = TRUE))
  } else {
    NA_real_
  }

  nb_tbl <- map_dfr(thresholds, function(t) {
    NB <- net_benefit(y = y, p = p, t = t)

    map_dfr(prior_m_values, function(m) {
      compute_lce_sd(y = y, p = p, t = t, prior_m = m) %>%
        mutate(
          threshold = t,
          NB = NB,
          WNB = NB / (1 + LCE + SD)
        )
    })
  })

  list(auc = AUC, nb = nb_tbl)
}

# ----------------------------------------------------------------------
# 7. Full Monte Carlo design
# ----------------------------------------------------------------------
simulation_grid <- tidyr::crossing(
  n = n_values,
  scenario = scenarios,
  sim_id = seq_len(n_sim)
)

# ----------------------------------------------------------------------
# 8. Run simulations
# ----------------------------------------------------------------------
with_progress({
  p_progress <- progressor(steps = nrow(simulation_grid))

  results <- future_pmap_dfr(
    simulation_grid,
    function(n, scenario, sim_id) {
      sim <- simulate_dataset(n = n, scenario = scenario)

      p_progress(sprintf(
        "n = %s | %s | replicate %s/%s",
        n, scenario, sim_id, n_sim
      ))

      tibble(
        n = n,
        scenario = scenario,
        sim_id = sim_id,
        AUC = sim$auc,
        nb = list(sim$nb)
      )
    },
    .options = furrr_options(seed = TRUE, scheduling = 2)
  )
})

# ----------------------------------------------------------------------
# 9. Unnest NB / WNB results
# ----------------------------------------------------------------------
nb_all <- results %>% unnest(nb)

# ----------------------------------------------------------------------
# 10. AUC summary
# Mean + empirical 95% simulation interval
# ----------------------------------------------------------------------
auc_summary <- results %>%
  group_by(n, scenario) %>%
  summarise(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_low = as.numeric(quantile(AUC, probs = 0.025, na.rm = TRUE)),
    AUC_high = as.numeric(quantile(AUC, probs = 0.975, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(n, factor(scenario, levels = scenarios))

print(auc_summary, n = Inf)

# ----------------------------------------------------------------------
# 11. NB / WNB summary
# Means + empirical 95% simulation intervals
# ----------------------------------------------------------------------
nb_summary_table <- nb_all %>%
  group_by(n, scenario, threshold, prior_m) %>%
  summarise(
    NB_mean = mean(NB, na.rm = TRUE),
    NB_low = as.numeric(quantile(NB, probs = 0.025, na.rm = TRUE)),
    NB_high = as.numeric(quantile(NB, probs = 0.975, na.rm = TRUE)),
    WNB_mean = mean(WNB, na.rm = TRUE),
    WNB_low = as.numeric(quantile(WNB, probs = 0.025, na.rm = TRUE)),
    WNB_high = as.numeric(quantile(WNB, probs = 0.975, na.rm = TRUE)),
    LCE_mean = mean(LCE, na.rm = TRUE),
    SD_mean = mean(SD, na.rm = TRUE),
    mean_n_bin = mean(n_bin, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n, factor(scenario, levels = scenarios), threshold, prior_m)

print(nb_summary_table, n = Inf)

# ----------------------------------------------------------------------
# 12. Main results: reference prior m = 10
# ----------------------------------------------------------------------
main_results <- nb_summary_table %>%
  filter(prior_m == reference_prior_m) %>%
  left_join(auc_summary, by = c("n", "scenario")) %>%
  mutate(
    relative_attenuation_percent = if_else(
      NB_mean > 0,
      100 * (NB_mean - WNB_mean) / NB_mean,
      NA_real_
    )
  ) %>%
  select(
    n, scenario, threshold,
    AUC_mean, AUC_low, AUC_high,
    NB_mean, NB_low, NB_high,
    WNB_mean, WNB_low, WNB_high,
    relative_attenuation_percent,
    LCE_mean, SD_mean, mean_n_bin
  ) %>%
  arrange(n, factor(scenario, levels = scenarios), threshold)

print(main_results, n = Inf)

# ----------------------------------------------------------------------
# 13. Formatted table for manuscript
# ----------------------------------------------------------------------
main_results_formatted <- main_results %>%
  mutate(
    AUC = sprintf("%.3f [%.3f–%.3f]", AUC_mean, AUC_low, AUC_high),
    NB = sprintf("%.4f [%.4f–%.4f]", NB_mean, NB_low, NB_high),
    WNB = sprintf("%.4f [%.4f–%.4f]", WNB_mean, WNB_low, WNB_high),
    `Relative attenuation (%)` = if_else(
      is.na(relative_attenuation_percent),
      NA_character_,
      sprintf("%.2f", relative_attenuation_percent)
    ),
    LCE = sprintf("%.4f", LCE_mean),
    SDpost = sprintf("%.4f", SD_mean),
    `Mean n in local bin` = sprintf("%.1f", mean_n_bin)
  ) %>%
  select(
    n,
    Scenario = scenario,
    Threshold = threshold,
    AUC, NB, WNB,
    `Relative attenuation (%)`,
    LCE, SDpost,
    `Mean n in local bin`
  )

print(main_results_formatted, n = Inf)

# ----------------------------------------------------------------------
# 14. Prior-sensitivity summaries
# ----------------------------------------------------------------------
wnb_reference_comparison <- nb_summary_table %>%
  select(n, scenario, threshold, prior_m, WNB_mean) %>%
  pivot_wider(
    names_from = prior_m,
    values_from = WNB_mean,
    names_prefix = "m_"
  ) %>%
  mutate(
    diff_m2_vs_m10 = m_2 - m_10,
    diff_m50_vs_m10 = m_50 - m_10,
    abs_diff_m2_vs_m10 = abs(m_2 - m_10),
    abs_diff_m50_vs_m10 = abs(m_50 - m_10)
  )

sensitivity_summary <- wnb_reference_comparison %>%
  group_by(n, scenario) %>%
  summarise(
    max_abs_diff_m2_vs_m10 = max(abs_diff_m2_vs_m10, na.rm = TRUE),
    threshold_max_diff_m2 = threshold[which.max(abs_diff_m2_vs_m10)],
    max_abs_diff_m50_vs_m10 = max(abs_diff_m50_vs_m10, na.rm = TRUE),
    threshold_max_diff_m50 = threshold[which.max(abs_diff_m50_vs_m10)],
    .groups = "drop"
  )

penalty_summary <- nb_summary_table %>%
  mutate(
    relative_attenuation_percent = if_else(
      NB_mean > 0,
      100 * (NB_mean - WNB_mean) / NB_mean,
      NA_real_
    )
  )

print(wnb_reference_comparison, n = Inf)
print(sensitivity_summary, n = Inf)
print(penalty_summary, n = Inf)

# ----------------------------------------------------------------------
# 15. Plot data
# ----------------------------------------------------------------------
plot_data <- nb_summary_table %>%
  mutate(
    prior_m = factor(
      prior_m,
      levels = c(2, 10, 50),
      labels = c("m = 2", "m = 10", "m = 50")
    ),
    n = factor(n, levels = n_values, labels = paste0("n = ", n_values)),
    scenario = factor(scenario, levels = scenarios)
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
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.6) +
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

print(miscalibration_plot)

# ----------------------------------------------------------------------
# 17. Figure 2: LCE and posterior SD sensitivity
# ----------------------------------------------------------------------
lce_sensitivity_plot <- ggplot(
  plot_data,
  aes(x = threshold, y = LCE_mean, color = prior_m, group = prior_m)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_grid(scenario ~ n) +
  scale_x_continuous(breaks = thresholds) +
  labs(
    x = "Decision threshold",
    y = "Local Calibration Error",
    color = "Prior concentration"
  ) +
  theme_bw(base_size = 14)

sd_sensitivity_plot <- ggplot(
  plot_data,
  aes(x = threshold, y = SD_mean, color = prior_m, group = prior_m)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_grid(scenario ~ n) +
  scale_x_continuous(breaks = thresholds) +
  labs(
    x = "Decision threshold",
    y = "Posterior SD",
    color = "Prior concentration"
  ) +
  theme_bw(base_size = 14)

fig_prior_components <- ggarrange(
  lce_sensitivity_plot,
  sd_sensitivity_plot,
  labels = c("a", "b"),
  ncol = 2,
  nrow = 1,
  align = "h",
  common.legend = TRUE,
  legend = "bottom"
)

print(lce_sensitivity_plot)
print(sd_sensitivity_plot)
print(fig_prior_components)

# ----------------------------------------------------------------------
# 18. Figure 3: change in WNB relative to reference m = 10
# ----------------------------------------------------------------------
delta_plot_data <- wnb_reference_comparison %>%
  select(
    n, scenario, threshold,
    diff_m2_vs_m10, diff_m50_vs_m10
  ) %>%
  pivot_longer(
    cols = c(diff_m2_vs_m10, diff_m50_vs_m10),
    names_to = "comparison",
    values_to = "delta_WNB"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      "diff_m2_vs_m10" = "m = 2 vs m = 10",
      "diff_m50_vs_m10" = "m = 50 vs m = 10"
    ),
    comparison = factor(
      comparison,
      levels = c("m = 2 vs m = 10", "m = 50 vs m = 10")
    ),
    n = factor(n, levels = n_values, labels = paste0("n = ", n_values)),
    scenario = factor(scenario, levels = scenarios)
  )

delta_wnb_plot <- ggplot(
  delta_plot_data,
  aes(x = threshold, y = delta_WNB, color = comparison, group = comparison)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_grid(scenario ~ n) +
  scale_x_continuous(breaks = thresholds) +
  labs(
    x = "Decision threshold",
    y = expression(Delta * " WNB relative to m = 10"),
    color = ""
  ) +
  theme_bw(base_size = 14)

print(delta_wnb_plot)

# ----------------------------------------------------------------------
# 19. Figure 4: relative attenuation of NB
# ----------------------------------------------------------------------
penalty_plot_data <- penalty_summary %>%
  mutate(
    prior_m = factor(
      prior_m,
      levels = c(2, 10, 50),
      labels = c("m = 2", "m = 10", "m = 50")
    ),
    n = factor(n, levels = n_values, labels = paste0("n = ", n_values)),
    scenario = factor(scenario, levels = scenarios)
  )

penalty_plot <- ggplot(
  penalty_plot_data,
  aes(
    x = threshold,
    y = relative_attenuation_percent,
    color = prior_m,
    group = prior_m
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_grid(scenario ~ n) +
  scale_x_continuous(breaks = thresholds) +
  labs(
    x = "Decision threshold",
    y = "Relative attenuation of NB (%)",
    color = "Prior concentration"
  ) +
  theme_bw(base_size = 14)

print(penalty_plot)
