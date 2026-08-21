# ==============================================================================
# SIMULATION 2: Logistic-like vs RF-like model
# ==============================================================================
# Single simulated dataset (n = 4000)
# Bootstrap percentile 95% CIs for AUC, Net Benefit (NB), and
# Weighted Net Benefit (WNB)
# Bootstrap replicates: 500
# Prior sensitivity: m = 2, 10, 50; prior mean fixed at 0.15 (i.e., global prevalence)
# ==============================================================================

# Packages ---------------------------------------------------------------------
library(tidyverse)
library(pROC)
library(furrr)
library(progressr)
library(ggpubr)

set.seed(789)

# Parallel processing ----------------------------------------------------------
plan(multisession, workers = max(1, parallel::detectCores() - 1))
handlers(global = TRUE)
handlers("txtprogressbar")

# Settings ---------------------------------------------------------------------
n <- 4000L
target_prev <- 0.15
prior_m_values <- c(2, 10, 50)
reference_prior_m <- 10
thr_centers <- seq(0.1, 0.5, by = 0.1)
bin_half_width <- 0.05
B <- 500L

# Helper functions -------------------------------------------------------------
net_benefit <- function(y, p, thr) {
  tp <- mean(p >= thr & y == 1)
  fp <- mean(p >= thr & y == 0)
  tp - fp * thr / (1 - thr)
}

compute_lce_sd_bin <- function(y, p, lower, upper, prior_m) {
  idx <- which(p >= lower & p < upper)
  n_bin <- length(idx)
  prior_alpha <- prior_m * target_prev
  prior_beta <- prior_m * (1 - target_prev)

  if (n_bin == 0) {
    return(tibble(
      prior_m = prior_m,
      prior_alpha = prior_alpha,
      prior_beta = prior_beta,
      n_bin = 0L,
      k = 0L,
      pbar = NA_real_,
      post_mean = NA_real_,
      prior_weight = 1,
      LCE = NA_real_,
      post_sd = NA_real_
    ))
  }

  k <- sum(y[idx])
  pbar <- mean(p[idx])
  alpha_post <- prior_alpha + k
  beta_post <- prior_beta + n_bin - k
  post_mean <- alpha_post / (alpha_post + beta_post)
  post_sd <- sqrt(
    alpha_post * beta_post /
      ((alpha_post + beta_post)^2 * (alpha_post + beta_post + 1))
  )
  prior_weight <- prior_m / (n_bin + prior_m)

  tibble(
    prior_m = prior_m,
    prior_alpha = prior_alpha,
    prior_beta = prior_beta,
    n_bin = n_bin,
    k = k,
    pbar = pbar,
    post_mean = post_mean,
    prior_weight = prior_weight,
    LCE = abs(pbar - post_mean),
    post_sd = post_sd
  )
}

# Threshold-centered bins ------------------------------------------------------
bins <- tibble(
  t = thr_centers,
  lower = thr_centers - bin_half_width,
  upper = thr_centers + bin_half_width
)

calc_curve <- function(y, p, label) {
  map_dfr(seq_len(nrow(bins)), function(i) {
    thr <- bins$t[i]
    lower <- bins$lower[i]
    upper <- bins$upper[i]
    NB <- net_benefit(y = y, p = p, thr = thr)

    map_dfr(prior_m_values, function(m) {
      compute_lce_sd_bin(
        y = y,
        p = p,
        lower = lower,
        upper = upper,
        prior_m = m
      ) %>%
        mutate(
          t = thr,
          NB = NB,
          WNB = NB / (1 + LCE + post_sd),
          model = label
        )
    })
  })
}

# Simulate nonlinear data-generating process ----------------------------------
X <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
colnames(X) <- paste0("x", 1:5)

lp_true <-
  2.5 * sin(X[, 1]) +
  1.8 * (X[, 2]^2) -
  0.9 * X[, 3] +
  0.5 * X[, 4] * X[, 5]

find_intercept <- function(offset) {
  mean(plogis(lp_true + offset)) - target_prev
}

intercept <- uniroot(find_intercept, c(-10, 10))$root
lp_true <- lp_true + intercept
p_true <- plogis(lp_true)
y <- rbinom(n = n, size = 1, prob = p_true)

dat <- tibble(
  y = y,
  x1 = X[, 1],
  x2 = X[, 2],
  x3 = X[, 3],
  x4 = X[, 4],
  x5 = X[, 5],
  p_true = p_true
)

# Logistic-like model ----------------------------------------------------------
model_logit <- glm(y ~ x1 + x2 + x3, data = dat, family = binomial())

# Noise is generated once and resampled with observations during bootstrap.
dat$logit_noise <- rnorm(n, mean = 0, sd = 0.015)
dat$p_logit <- predict(model_logit, newdata = dat, type = "response") + dat$logit_noise
dat$p_logit <- pmin(pmax(dat$p_logit, 0.001), 0.5)

# RF-like model ----------------------------------------------------------------
p_rf_base <- dat$p_true + rnorm(n, mean = 0, sd = 0.03)
p_rf_base <- pmin(pmax(p_rf_base, 0.001), 0.999)
lp_compressed <- 0.5 * qlogis(p_rf_base)
dat$p_rf <- plogis(lp_compressed)
dat$p_rf <- pmin(pmax(dat$p_rf, 0.001), 0.999)

# Original point estimates -----------------------------------------------------
AUC_original <- tibble(
  model = c("Logistic", "RF"),
  AUC = c(
    as.numeric(pROC::auc(dat$y, dat$p_logit, quiet = TRUE)),
    as.numeric(pROC::auc(dat$y, dat$p_rf, quiet = TRUE))
  )
)

results_original <- bind_rows(
  calc_curve(dat$y, dat$p_logit, "Logistic"),
  calc_curve(dat$y, dat$p_rf, "RF")
)

cat("\nObserved prevalence:", round(mean(dat$y), 3), "\n")
print(AUC_original)
print(bins)
print(results_original, n = Inf)

# Bootstrap --------------------------------------------------------------------
bootstrap_once <- function(b) {
  idx <- sample.int(nrow(dat), size = nrow(dat), replace = TRUE)
  boot_dat <- dat[idx, ]

  # Refit logistic model in the bootstrap sample.
  boot_logit <- glm(y ~ x1 + x2 + x3, data = boot_dat, family = binomial())
  p_logit_boot <-
    predict(boot_logit, newdata = boot_dat, type = "response") +
    boot_dat$logit_noise
  p_logit_boot <- pmin(pmax(p_logit_boot, 0.001), 0.5)

  # No fitted RF exists in this simulation; RF-like predictions are resampled.
  p_rf_boot <- boot_dat$p_rf

  auc_logit <- if (length(unique(boot_dat$y)) == 2) {
    as.numeric(pROC::auc(boot_dat$y, p_logit_boot, quiet = TRUE))
  } else {
    NA_real_
  }

  auc_rf <- if (length(unique(boot_dat$y)) == 2) {
    as.numeric(pROC::auc(boot_dat$y, p_rf_boot, quiet = TRUE))
  } else {
    NA_real_
  }

  auc_tbl <- tibble(
    bootstrap_id = b,
    model = c("Logistic", "RF"),
    AUC = c(auc_logit, auc_rf)
  )

  nb_tbl <- bind_rows(
    calc_curve(boot_dat$y, p_logit_boot, "Logistic"),
    calc_curve(boot_dat$y, p_rf_boot, "RF")
  ) %>%
    mutate(bootstrap_id = b)

  list(auc = auc_tbl, nb = nb_tbl)
}

cat("\nRunning", B, "bootstrap replicates...\n\n")

with_progress({
  p_progress <- progressor(steps = B)
  bootstrap_results <- future_map(
    seq_len(B),
    function(b) {
      out <- bootstrap_once(b)
      p_progress(sprintf("bootstrap %d/%d", b, B))
      out
    },
    .options = furrr_options(seed = TRUE, scheduling = 2)
  )
})

auc_boot_all <- map_dfr(bootstrap_results, "auc")
nb_boot_all <- map_dfr(bootstrap_results, "nb")

# Bootstrap percentile 95% CIs -------------------------------------------------
auc_boot_ci <- auc_boot_all %>%
  group_by(model) %>%
  summarise(
    AUC_low = as.numeric(quantile(AUC, probs = 0.025, na.rm = TRUE)),
    AUC_high = as.numeric(quantile(AUC, probs = 0.975, na.rm = TRUE)),
    .groups = "drop"
  )

auc_summary <- AUC_original %>%
  rename(AUC_estimate = AUC) %>%
  left_join(auc_boot_ci, by = "model")

nb_boot_ci <- nb_boot_all %>%
  group_by(model, t, prior_m) %>%
  summarise(
    NB_low = as.numeric(quantile(NB, probs = 0.025, na.rm = TRUE)),
    NB_high = as.numeric(quantile(NB, probs = 0.975, na.rm = TRUE)),
    WNB_low = as.numeric(quantile(WNB, probs = 0.025, na.rm = TRUE)),
    WNB_high = as.numeric(quantile(WNB, probs = 0.975, na.rm = TRUE)),
    valid_WNB_boot = sum(!is.na(WNB)),
    .groups = "drop"
  )

results_with_ci <- results_original %>%
  left_join(nb_boot_ci, by = c("model", "t", "prior_m"))

# Reference prior (m = 10) -----------------------------------------------------
results_reference <- results_with_ci %>%
  filter(prior_m == reference_prior_m) %>%
  left_join(auc_summary, by = "model") %>%
  mutate(
    relative_attenuation_percent = if_else(
      NB > 0,
      100 * (NB - WNB) / NB,
      NA_real_
    )
  ) %>%
  select(
    model,
    t,
    AUC_estimate,
    AUC_low,
    AUC_high,
    NB,
    NB_low,
    NB_high,
    WNB,
    WNB_low,
    WNB_high,
    relative_attenuation_percent,
    LCE,
    post_sd,
    n_bin,
    valid_WNB_boot
  )

results_reference_formatted <- results_reference %>%
  mutate(
    AUC = sprintf("%.3f [%.3f–%.3f]", AUC_estimate, AUC_low, AUC_high),
    NB = sprintf("%.4f [%.4f–%.4f]", NB, NB_low, NB_high),
    WNB = sprintf("%.4f [%.4f–%.4f]", WNB, WNB_low, WNB_high),
    `Relative attenuation (%)` = if_else(
      is.na(relative_attenuation_percent),
      NA_character_,
      sprintf("%.2f", relative_attenuation_percent)
    )
  ) %>%
  select(
    Model = model,
    Threshold = t,
    AUC,
    NB,
    WNB,
    `Relative attenuation (%)`
  )

print(auc_summary)
print(results_with_ci, n = Inf)
print(results_reference, n = Inf)
print(results_reference_formatted, n = Inf)

# Prior sensitivity ------------------------------------------------------------
wnb_reference_comparison <- results_original %>%
  select(model, t, prior_m, WNB) %>%
  pivot_wider(names_from = prior_m, values_from = WNB, names_prefix = "m_") %>%
  mutate(
    diff_m2_vs_m10 = m_2 - m_10,
    diff_m50_vs_m10 = m_50 - m_10,
    abs_diff_m2_vs_m10 = abs(diff_m2_vs_m10),
    abs_diff_m50_vs_m10 = abs(diff_m50_vs_m10)
  )

sensitivity_summary <- wnb_reference_comparison %>%
  group_by(model) %>%
  summarise(
    max_abs_diff_m2_vs_m10 = max(abs_diff_m2_vs_m10, na.rm = TRUE),
    threshold_max_m2 = t[which.max(abs_diff_m2_vs_m10)],
    max_abs_diff_m50_vs_m10 = max(abs_diff_m50_vs_m10, na.rm = TRUE),
    threshold_max_m50 = t[which.max(abs_diff_m50_vs_m10)],
    .groups = "drop"
  )

print(sensitivity_summary)

# Prior-sensitivity figures ----------------------------------------------------
plot_data <- results_original %>%
  mutate(
    prior_m = factor(
      prior_m,
      levels = c(2, 10, 50),
      labels = c("m = 2", "m = 10", "m = 50")
    ),
    model = factor(model, levels = c("Logistic", "RF"))
  )

p_wnb_sensitivity <- ggplot(
  plot_data,
  aes(x = t, y = WNB, color = prior_m, group = prior_m)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_wrap(~model, nrow = 1) +
  scale_x_continuous(breaks = thr_centers) +
  labs(
    x = "Decision threshold",
    y = "Weighted Net Benefit",
    color = "Prior concentration"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

p_lce_sensitivity <- ggplot(
  plot_data,
  aes(x = t, y = LCE, color = prior_m, group = prior_m)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_wrap(~model, nrow = 1) +
  scale_x_continuous(breaks = thr_centers) +
  labs(
    x = "Decision threshold",
    y = "Local Calibration Error",
    color = "Prior concentration"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

p_sd_sensitivity <- ggplot(
  plot_data,
  aes(x = t, y = post_sd, color = prior_m, group = prior_m)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_wrap(~model, nrow = 1) +
  scale_x_continuous(breaks = thr_centers) +
  labs(
    x = "Decision threshold",
    y = "Posterior SD",
    color = "Prior concentration"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

fig_components_sensitivity <- ggarrange(
  p_lce_sensitivity,
  p_sd_sensitivity,
  labels = c("a", "b"),
  ncol = 2,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom",
  align = "h"
)

delta_plot_data <- wnb_reference_comparison %>%
  select(model, t, diff_m2_vs_m10, diff_m50_vs_m10) %>%
  pivot_longer(
    cols = c(diff_m2_vs_m10, diff_m50_vs_m10),
    names_to = "comparison",
    values_to = "delta_WNB"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      diff_m2_vs_m10 = "m = 2 vs m = 10",
      diff_m50_vs_m10 = "m = 50 vs m = 10"
    )
  )

p_delta_wnb <- ggplot(
  delta_plot_data,
  aes(x = t, y = delta_WNB, color = comparison, group = comparison)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_wrap(~model, nrow = 1) +
  scale_x_continuous(breaks = thr_centers) +
  labs(
    x = "Decision threshold",
    y = expression(Delta * " WNB relative to m = 10"),
    color = ""
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

# Mean NB/WNB across bootstrap replicates --------------------------------------
plot_nb_wnb <- nb_boot_all %>%
  filter(prior_m == reference_prior_m) %>%
  group_by(model, t) %>%
  summarise(
    NB = mean(NB, na.rm = TRUE),
    WNB = mean(WNB, na.rm = TRUE),
    .groups = "drop"
  )

prev <- mean(dat$y)
df_treat <- tibble(
  t = thr_centers,
  `Treat all` = prev - (1 - prev) * t / (1 - t),
  `Treat none` = 0
) %>%
  pivot_longer(
    cols = c(`Treat all`, `Treat none`),
    names_to = "strategy",
    values_to = "NB"
  )

plot_long <- plot_nb_wnb %>%
  pivot_longer(cols = c(NB, WNB), names_to = "metric", values_to = "value")

p_nb_wnb_boot <- ggplot() +
  geom_line(
    data = plot_long,
    aes(
      x = t,
      y = value,
      color = model,
      linetype = metric,
      group = interaction(model, metric)
    ),
    linewidth = 1.1
  ) +
  geom_line(
    data = df_treat,
    aes(x = t, y = NB, color = strategy, linetype = strategy, group = strategy),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c(
      Logistic = "orange",
      RF = "green4",
      `Treat all` = "red4",
      `Treat none` = "gray40"
    ),
    breaks = c("Logistic", "RF", "Treat all", "Treat none")
  ) +
  scale_linetype_manual(
    values = c(
      NB = "solid",
      WNB = "dashed",
      `Treat all` = "dotdash",
      `Treat none` = "dotdash"
    ),
    breaks = c("NB", "WNB", "Treat all", "Treat none")
  ) +
  scale_x_continuous(breaks = thr_centers, limits = c(0.08, 0.52)) +
  coord_cartesian(ylim = c(-0.025, 0.105)) +
  labs(
    x = "Decision threshold",
    y = "Net Benefit",
    color = NULL,
    linetype = NULL,
    title = "Net Benefit (solid) vs Weighted Net Benefit (dashed)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    legend.box = "horizontal"
  )

print(plot_nb_wnb)
print(p_nb_wnb_boot)

# Mean LCE and posterior SD across bootstrap replicates ------------------------
plot_components_boot <- nb_boot_all %>%
  filter(prior_m == reference_prior_m) %>%
  group_by(model, t) %>%
  summarise(
    LCE = mean(LCE, na.rm = TRUE),
    post_sd = mean(post_sd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = factor(model, levels = c("Logistic", "RF")))

p_lce_boot <- ggplot(
  plot_components_boot,
  aes(x = t, y = LCE, color = model, group = model)
) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c(Logistic = "orange", RF = "green4")) +
  scale_x_continuous(breaks = thr_centers, limits = c(0.08, 0.52)) +
  coord_cartesian(ylim = c(0, 0.16)) +
  labs(x = "Predicted probability", y = "Local Calibration Error (LCE)", color = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", legend.title = element_blank())

p_sd_boot <- ggplot(
  plot_components_boot,
  aes(x = t, y = post_sd, color = model, group = model)
) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c(Logistic = "orange", RF = "green4")) +
  scale_x_continuous(breaks = thr_centers, limits = c(0.08, 0.52)) +
  coord_cartesian(ylim = c(0, 0.16)) +
  labs(x = "Predicted probability", y = "Posterior SD", color = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", legend.title = element_blank())

fig_components_boot <- ggarrange(
  p_lce_boot,
  p_sd_boot,
  ncol = 2,
  nrow = 1,
  common.legend = FALSE,
  align = "h"
)

print(plot_components_boot)
print(fig_components_boot)
