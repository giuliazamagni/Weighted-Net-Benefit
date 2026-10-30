# ==============================================================================
# SIMULATION 2
# Logistic-like vs RF-like model
#
# ------------------------------------------------------------------------------
# 1. Packages and reproducibility
# ------------------------------------------------------------------------------

library(tidyverse)
library(pROC)
library(furrr)
library(progressr)
library(ggpubr)

set.seed(789)

plan(
  multisession,
  workers = max(1, parallel::detectCores() - 1)
)

handlers(global = TRUE)
handlers("txtprogressbar")


# ------------------------------------------------------------------------------
# 2. Simulation settings
# ------------------------------------------------------------------------------

n <- 4000L
target_prev <- 0.15

thresholds <- seq(0.10, 0.50, by = 0.10)
bin_half_width <- 0.05

# Reference prior:
#   p0 = 0.15
#   m  = 10
#   Beta(1.5, 8.5)
reference_prior_m <- 10

# Prior-sensitivity analysis:
#   m = 2, 4, ..., 50
#   prior mean fixed at 0.15
prior_m_values <- seq(2, 50, by = 2)

# Number of bootstrap replicates used for uncertainty summaries
B <- 500L

output_dir <- "simulation2_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(
    quantile(
      x,
      probs = prob,
      names = FALSE,
      type = 7
    )
  )
}

net_benefit <- function(y, p, threshold) {

  tp <- mean(p >= threshold & y == 1)
  fp <- mean(p >= threshold & y == 0)

  tp - fp * threshold / (1 - threshold)
}


# ------------------------------------------------------------------------------
# 4. Local posterior quantities
# ------------------------------------------------------------------------------

compute_local_metrics <- function(
    y,
    p,
    threshold,
    prior_m
) {

  lower <- max(0, threshold - bin_half_width)
  upper <- min(1, threshold + bin_half_width)

  idx <- which(
    p >= lower &
      p < upper
  )

  n_bin <- length(idx)

  prior_alpha <- prior_m * target_prev
  prior_beta <- prior_m * (1 - target_prev)

  # A local calibration estimate is not available when the bin is empty.
  # Return NA rather than treating calibration error or uncertainty as zero.
  if (n_bin == 0L) {

    return(
      tibble(
        threshold = threshold,
        prior_m = prior_m,
        prior_alpha = prior_alpha,
        prior_beta = prior_beta,
        n_bin = 0L,
        k = 0L,
        pbar = NA_real_,
        empirical_rate = NA_real_,
        posterior_rate = NA_real_,
        prior_weight = NA_real_,
        signed_prior_pull = NA_real_,
        abs_prior_pull = NA_real_,
        LCE = NA_real_,
        SDpost = NA_real_
      )
    )
  }

  k <- sum(y[idx])
  pbar <- mean(p[idx])

  empirical_rate <- k / n_bin

  alpha_post <- prior_alpha + k
  beta_post <- prior_beta + n_bin - k

  posterior_rate <- alpha_post / (alpha_post + beta_post)

  SDpost <- sqrt(
    alpha_post * beta_post /
      (
        (alpha_post + beta_post)^2 *
          (alpha_post + beta_post + 1)
      )
  )

  prior_weight <- prior_m / (n_bin + prior_m)

  signed_prior_pull <- posterior_rate - empirical_rate
  abs_prior_pull <- abs(signed_prior_pull)

  LCE <- abs(pbar - posterior_rate)

  tibble(
    threshold = threshold,
    prior_m = prior_m,
    prior_alpha = prior_alpha,
    prior_beta = prior_beta,
    n_bin = n_bin,
    k = k,
    pbar = pbar,
    empirical_rate = empirical_rate,
    posterior_rate = posterior_rate,
    prior_weight = prior_weight,
    signed_prior_pull = signed_prior_pull,
    abs_prior_pull = abs_prior_pull,
    LCE = LCE,
    SDpost = SDpost
  )
}


# ------------------------------------------------------------------------------
# 5. Calculate NB, LCE, SDpost and WNB for one prediction vector
# ------------------------------------------------------------------------------

calculate_model_results <- function(
    y,
    p,
    model_label,
    prior_values = prior_m_values
) {

  map_dfr(
    thresholds,
    function(threshold) {

      NB <- net_benefit(
        y = y,
        p = p,
        threshold = threshold
      )

      map_dfr(
        prior_values,
        function(m) {

          compute_local_metrics(
            y = y,
            p = p,
            threshold = threshold,
            prior_m = m
          ) %>%
            mutate(
              model = model_label,
              NB = NB,

              # WNB definition used in the manuscript
              WNB = NB / (1 + LCE + SDpost)
            )
        }
      )
    }
  )
}


# ------------------------------------------------------------------------------
# 6. Simulate nonlinear data-generating process
# ------------------------------------------------------------------------------

X <- matrix(
  rnorm(n * 5),
  nrow = n,
  ncol = 5
)

colnames(X) <- paste0("x", 1:5)

lp_true <-
  2.5 * sin(X[, 1]) +
  1.8 * (X[, 2]^2) -
  0.9 * X[, 3] +
  0.5 * X[, 4] * X[, 5]

find_intercept <- function(offset) {
  mean(
    plogis(
      lp_true + offset
    )
  ) - target_prev
}

intercept <- uniroot(
  find_intercept,
  c(-10, 10)
)$root

lp_true <- lp_true + intercept
p_true <- plogis(lp_true)

y <- rbinom(
  n = n,
  size = 1,
  prob = p_true
)

dat <- tibble(
  y = y,
  x1 = X[, 1],
  x2 = X[, 2],
  x3 = X[, 3],
  x4 = X[, 4],
  x5 = X[, 5],
  p_true = p_true
)


# ------------------------------------------------------------------------------
# 7. Logistic-like model
# ------------------------------------------------------------------------------

model_logit <- glm(
  y ~ x1 + x2 + x3,
  data = dat,
  family = binomial()
)

# The perturbation is generated once and stored with each observation.
# During bootstrap resampling it therefore travels with the resampled subject.
dat$logit_noise <- rnorm(
  n,
  mean = 0,
  sd = 0.015
)

dat$p_logit <-
  predict(
    model_logit,
    newdata = dat,
    type = "response"
  ) +
  dat$logit_noise

dat$p_logit <- pmin(
  pmax(
    dat$p_logit,
    0.001
  ),
  0.5
)


# ------------------------------------------------------------------------------
# 8. RF-like model
# ------------------------------------------------------------------------------

# This is an RF-like prediction profile rather than a fitted random-forest model.
# It is constructed from the true nonlinear risk plus noise and then compressed
# on the logit scale to induce the intended S-shaped miscalibration pattern.

p_rf_base <-
  dat$p_true +
  rnorm(
    n,
    mean = 0,
    sd = 0.03
  )

p_rf_base <- pmin(
  pmax(
    p_rf_base,
    0.001
  ),
  0.999
)

lp_compressed <-
  0.5 *
  qlogis(
    p_rf_base
  )

dat$p_rf <-
  plogis(
    lp_compressed
  )

dat$p_rf <- pmin(
  pmax(
    dat$p_rf,
    0.001
  ),
  0.999
)


# ------------------------------------------------------------------------------
# 9. Original-data results for all prior concentrations
# ------------------------------------------------------------------------------

AUC_original <- tibble(
  model = c("Logistic", "RF"),
  AUC_original = c(
    as.numeric(
      pROC::auc(
        dat$y,
        dat$p_logit,
        quiet = TRUE
      )
    ),
    as.numeric(
      pROC::auc(
        dat$y,
        dat$p_rf,
        quiet = TRUE
      )
    )
  )
)

results_all <- bind_rows(
  calculate_model_results(
    y = dat$y,
    p = dat$p_logit,
    model_label = "Logistic",
    prior_values = prior_m_values
  ),
  calculate_model_results(
    y = dat$y,
    p = dat$p_rf,
    model_label = "RF",
    prior_values = prior_m_values
  )
)

cat(
  "\nObserved event prevalence:",
  round(mean(dat$y), 4),
  "\n\n"
)

print(AUC_original)
print(results_all, n = Inf)


# ------------------------------------------------------------------------------
# 10. Prior sensitivity relative to m = 10
# ------------------------------------------------------------------------------

reference_WNB <- results_all %>%
  filter(
    prior_m == reference_prior_m
  ) %>%
  select(
    model,
    threshold,
    WNB_reference = WNB
  )

prior_sensitivity_complete <- results_all %>%
  left_join(
    reference_WNB,
    by = c(
      "model",
      "threshold"
    )
  ) %>%
  mutate(
    delta_WNB = WNB - WNB_reference,
    abs_delta_WNB = abs(delta_WNB),

    prior_weight_percent = 100 * prior_weight,

    prior_pull_pp = 100 * signed_prior_pull,
    abs_prior_pull_pp = 100 * abs_prior_pull
  ) %>%
  select(
    model,
    threshold,
    prior_m,
    n_bin,
    pbar,
    empirical_rate,
    posterior_rate,
    prior_weight_percent,
    prior_pull_pp,
    abs_prior_pull_pp,
    LCE,
    SDpost,
    NB,
    WNB,
    delta_WNB,
    abs_delta_WNB
  ) %>%
  arrange(
    factor(
      model,
      levels = c("Logistic", "RF")
    ),
    threshold,
    prior_m
  )

write.csv(
  prior_sensitivity_complete,
  file.path(
    output_dir,
    "Simulation2_prior_sensitivity_complete_results.csv"
  ),
  row.names = FALSE
)


# Selected values for manuscript/supplementary table: m = 2, 10, 50
prior_sensitivity_selected <- prior_sensitivity_complete %>%
  filter(
    prior_m %in% c(2, 10, 50)
  ) %>%
  arrange(
    factor(
      model,
      levels = c("Logistic", "RF")
    ),
    threshold,
    prior_m
  )

write.csv(
  prior_sensitivity_selected,
  file.path(
    output_dir,
    "Simulation2_prior_sensitivity_selected_m2_m10_m50.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 11. Bootstrap under the reference prior m = 10
# ------------------------------------------------------------------------------

bootstrap_once <- function(b) {

  idx <- sample.int(
    nrow(dat),
    size = nrow(dat),
    replace = TRUE
  )

  boot_dat <- dat[idx, ]

  # Logistic-like model:
  # refit the logistic regression in the bootstrap sample, then add the
  # observation-level perturbation carried with the resampled observations.
  boot_logit <- glm(
    y ~ x1 + x2 + x3,
    data = boot_dat,
    family = binomial()
  )

  p_logit_boot <-
    predict(
      boot_logit,
      newdata = boot_dat,
      type = "response"
    ) +
    boot_dat$logit_noise

  p_logit_boot <- pmin(
    pmax(
      p_logit_boot,
      0.001
    ),
    0.5
  )

  # RF-like model:
  # there is no separately fitted RF object in this illustrative simulation.
  # The pre-generated RF-like predictions are therefore resampled together
  # with the observations.
  p_rf_boot <- boot_dat$p_rf

  auc_logit <- if (length(unique(boot_dat$y)) == 2L) {
    as.numeric(
      pROC::auc(
        boot_dat$y,
        p_logit_boot,
        quiet = TRUE
      )
    )
  } else {
    NA_real_
  }

  auc_rf <- if (length(unique(boot_dat$y)) == 2L) {
    as.numeric(
      pROC::auc(
        boot_dat$y,
        p_rf_boot,
        quiet = TRUE
      )
    )
  } else {
    NA_real_
  }

  auc_tbl <- tibble(
    bootstrap_id = b,
    model = c("Logistic", "RF"),
    AUC = c(
      auc_logit,
      auc_rf
    )
  )

  nb_tbl <- bind_rows(
    calculate_model_results(
      y = boot_dat$y,
      p = p_logit_boot,
      model_label = "Logistic",
      prior_values = reference_prior_m
    ),
    calculate_model_results(
      y = boot_dat$y,
      p = p_rf_boot,
      model_label = "RF",
      prior_values = reference_prior_m
    )
  ) %>%
    mutate(
      bootstrap_id = b
    )

  list(
    auc = auc_tbl,
    nb = nb_tbl
  )
}


cat(
  "\nRunning",
  B,
  "bootstrap replicates under the reference prior m = 10...\n\n"
)

with_progress({

  p_progress <- progressor(
    steps = B
  )

  bootstrap_results <- future_map(
    seq_len(B),
    function(b) {

      out <- bootstrap_once(b)

      p_progress(
        sprintf(
          "bootstrap %d/%d",
          b,
          B
        )
      )

      out
    },
    .options = furrr_options(
      seed = TRUE,
      scheduling = 2
    )
  )
})

plan(sequential)

auc_boot_all <- map_dfr(
  bootstrap_results,
  "auc"
)

nb_boot_all <- map_dfr(
  bootstrap_results,
  "nb"
)


# ------------------------------------------------------------------------------
# 12. Bootstrap summaries
#
# Intervals are percentile-based 95% bootstrap intervals.
# ------------------------------------------------------------------------------

auc_bootstrap_summary <- auc_boot_all %>%
  group_by(model) %>%
  summarise(
    AUC_mean = safe_mean(AUC),
    AUC_low = safe_quantile(AUC, 0.025),
    AUC_high = safe_quantile(AUC, 0.975),
    .groups = "drop"
  )

main_bootstrap_summary <- nb_boot_all %>%
  group_by(
    model,
    threshold
  ) %>%
  summarise(
    NB_mean = safe_mean(NB),
    NB_low = safe_quantile(NB, 0.025),
    NB_high = safe_quantile(NB, 0.975),

    WNB_mean = safe_mean(WNB),
    WNB_low = safe_quantile(WNB, 0.025),
    WNB_high = safe_quantile(WNB, 0.975),

    LCE_mean = safe_mean(LCE),
    SDpost_mean = safe_mean(SDpost),
    mean_n_bin = safe_mean(n_bin),

    valid_WNB_boot = sum(is.finite(WNB)),

    .groups = "drop"
  ) %>%
  left_join(
    auc_bootstrap_summary,
    by = "model"
  ) %>%
  mutate(
    relative_attenuation_percent = if_else(
      NB_mean > 0,
      100 * (NB_mean - WNB_mean) / NB_mean,
      NA_real_
    )
  ) %>%
  arrange(
    factor(
      model,
      levels = c("Logistic", "RF")
    ),
    threshold
  )

write.csv(
  main_bootstrap_summary,
  file.path(
    output_dir,
    "Simulation2_main_bootstrap_summary_m10.csv"
  ),
  row.names = FALSE
)


# Formatted manuscript-facing summary
main_bootstrap_formatted <- main_bootstrap_summary %>%
  mutate(
    AUC = sprintf(
      "%.3f [%.3f–%.3f]",
      AUC_mean,
      AUC_low,
      AUC_high
    ),

    NB = sprintf(
      "%.4f [%.4f–%.4f]",
      NB_mean,
      NB_low,
      NB_high
    ),

    WNB = sprintf(
      "%.4f [%.4f–%.4f]",
      WNB_mean,
      WNB_low,
      WNB_high
    ),

    `Relative attenuation (%)` = if_else(
      is.na(relative_attenuation_percent),
      NA_character_,
      sprintf(
        "%.2f",
        relative_attenuation_percent
      )
    ),

    LCE = sprintf(
      "%.4f",
      LCE_mean
    ),

    SDpost = sprintf(
      "%.4f",
      SDpost_mean
    ),

    `Mean n in local bin` = sprintf(
      "%.1f",
      mean_n_bin
    )
  ) %>%
  select(
    Model = model,
    Threshold = threshold,
    AUC,
    NB,
    WNB,
    `Relative attenuation (%)`,
    LCE,
    SDpost,
    `Mean n in local bin`,
    valid_WNB_boot
  )

write.csv(
  main_bootstrap_formatted,
  file.path(
    output_dir,
    "Simulation2_main_bootstrap_summary_m10_formatted.csv"
  ),
  row.names = FALSE
)

print(
  main_bootstrap_summary,
  n = Inf
)

print(
  main_bootstrap_formatted,
  n = Inf
)


# ------------------------------------------------------------------------------
# 13. Figure 1: prior sensitivity of LCE and SDpost
# ------------------------------------------------------------------------------

plot_data <- prior_sensitivity_complete %>%
  mutate(
    model = factor(
      model,
      levels = c(
        "Logistic",
        "RF"
      )
    )
  )

p_lce <- ggplot(
  plot_data,
  aes(
    x = threshold,
    y = LCE,
    color = prior_m,
    group = prior_m
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 1.4
  ) +
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  scale_x_continuous(
    breaks = thresholds
  ) +
  scale_color_viridis_c(
    option = "C",
    end = 0.95
  ) +
  labs(
    x = "Decision threshold",
    y = "Local Calibration Error (LCE)",
    color = "Prior concentration (m)"
  ) +
  theme_bw(
    base_size = 12
  )

p_sd <- ggplot(
  plot_data,
  aes(
    x = threshold,
    y = SDpost,
    color = prior_m,
    group = prior_m
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 1.4
  ) +
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  scale_x_continuous(
    breaks = thresholds
  ) +
  scale_color_viridis_c(
    option = "C",
    end = 0.95
  ) +
  labs(
    x = "Decision threshold",
    y = "Posterior SD",
    color = "Prior concentration (m)"
  ) +
  theme_bw(
    base_size = 12
  )

fig_prior_components <- ggarrange(
  p_lce,
  p_sd,
  labels = c(
    "a",
    "b"
  ),
  ncol = 2,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom",
  align = "h"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_prior_components_m2_to_m50_simulation2.pdf"
  ),
  plot = fig_prior_components,
  width = 12,
  height = 6
)


# ------------------------------------------------------------------------------
# 14. Figure 2: delta WNB relative to m = 10
# ------------------------------------------------------------------------------

p_delta_WNB <- prior_sensitivity_complete %>%
  filter(
    prior_m != reference_prior_m
  ) %>%
  mutate(
    model = factor(
      model,
      levels = c(
        "Logistic",
        "RF"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = threshold,
      y = delta_WNB,
      color = prior_m,
      group = prior_m
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 1.4
  ) +
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  scale_x_continuous(
    breaks = thresholds
  ) +
  scale_color_viridis_c(
    option = "C",
    end = 0.95
  ) +
  labs(
    x = "Decision threshold",
    y = expression(
      Delta *
        "WNB relative to m = 10"
    ),
    color = "Prior concentration (m)"
  ) +
  theme_bw(
    base_size = 12
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure_delta_WNB_vs_m10_m2_to_m50_simulation2.pdf"
  ),
  plot = p_delta_WNB,
  width = 10,
  height = 5.5
)


# ------------------------------------------------------------------------------
# 15. Optional figure: bootstrap mean NB versus WNB under m = 10
# ------------------------------------------------------------------------------

plot_nb_wnb <- main_bootstrap_summary %>%
  select(
    model,
    threshold,
    NB_mean,
    WNB_mean
  ) %>%
  pivot_longer(
    cols = c(
      NB_mean,
      WNB_mean
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      NB_mean = "NB",
      WNB_mean = "WNB"
    )
  )

observed_prev <- mean(dat$y)

df_treat <- tibble(
  threshold = thresholds,
  `Treat all` =
    observed_prev -
    (1 - observed_prev) *
    threshold /
    (1 - threshold),

  `Treat none` = 0
) %>%
  pivot_longer(
    cols = c(
      `Treat all`,
      `Treat none`
    ),
    names_to = "strategy",
    values_to = "value"
  )

p_nb_wnb_boot <- ggplot() +
  geom_line(
    data = plot_nb_wnb,
    aes(
      x = threshold,
      y = value,
      color = model,
      linetype = metric,
      group = interaction(
        model,
        metric
      )
    ),
    linewidth = 1.1
  ) +
  geom_line(
    data = df_treat,
    aes(
      x = threshold,
      y = value,
      color = strategy,
      linetype = strategy,
      group = strategy
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c(
      Logistic = "orange",
      RF = "green4",
      `Treat all` = "red4",
      `Treat none` = "gray40"
    ),
    breaks = c(
      "Logistic",
      "RF",
      "Treat all",
      "Treat none"
    )
  ) +
  scale_linetype_manual(
    values = c(
      NB = "solid",
      WNB = "dashed",
      `Treat all` = "dotdash",
      `Treat none` = "dotdash"
    ),
    breaks = c(
      "NB",
      "WNB",
      "Treat all",
      "Treat none"
    )
  ) +
  scale_x_continuous(
    breaks = thresholds,
    limits = c(
      0.08,
      0.52
    )
  ) +
  labs(
    x = "Decision threshold",
    y = "Net Benefit",
    color = NULL,
    linetype = NULL,
    title = "Net Benefit (solid) vs Weighted Net Benefit (dashed)"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(
      hjust = 0.5
    ),
    legend.box = "horizontal"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure_NB_vs_WNB_bootstrap_mean_m10_simulation2.pdf"
  ),
  plot = p_nb_wnb_boot,
  width = 9,
  height = 6
)


# ------------------------------------------------------------------------------
# 16. Reproducibility checks
# ------------------------------------------------------------------------------

cat(
  "\nPrior concentrations evaluated:",
  paste(
    prior_m_values,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Reference prior concentration:",
  reference_prior_m,
  "\n"
)

cat(
  "Reference prior:",
  sprintf(
    "Beta(%.1f, %.1f)",
    reference_prior_m * target_prev,
    reference_prior_m * (1 - target_prev)
  ),
  "\n"
)

cat(
  "Observed event prevalence:",
  round(
    mean(dat$y),
    4
  ),
  "\n"
)

cat(
  "Bootstrap replicates:",
  B,
  "\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)

cat(
  "\nSimulation 2 completed successfully.\n"
)
