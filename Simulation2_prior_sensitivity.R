# ==============================================================================
# SIMULATION 2: Logistic-like vs RF-like model
# Prior-sensitivity analysis
#
# Prior concentration:
#   m = 2, 4, 6, ..., 50
# Reference:
#   m = 10
# ==============================================================================

library(tidyverse)
library(pROC)
library(ggpubr)

set.seed(789)

# ------------------------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------------------------

n <- 4000L
target_prev <- 0.15

prior_m_values <- seq(2, 50, by = 2)
reference_prior_m <- 10

thresholds <- seq(0.10, 0.50, by = 0.10)
bin_half_width <- 0.05

output_dir <- "simulation2_prior_sensitivity_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Net Benefit
# ------------------------------------------------------------------------------

net_benefit <- function(y, p, threshold) {
  tp <- mean(p >= threshold & y == 1)
  fp <- mean(p >= threshold & y == 0)

  tp - fp * threshold / (1 - threshold)
}

# ------------------------------------------------------------------------------
# 3. Local posterior quantities
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

  if (n_bin == 0L) {
    return(
      tibble(
        threshold = threshold,
        prior_m = prior_m,
        n_bin = 0L,
        empirical_rate = NA_real_,
        posterior_rate = NA_real_,
        prior_weight = NA_real_,
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

  abs_prior_pull <- abs(
    posterior_rate - empirical_rate
  )

  LCE <- abs(
    pbar - posterior_rate
  )

  tibble(
    threshold = threshold,
    prior_m = prior_m,
    n_bin = n_bin,
    empirical_rate = empirical_rate,
    posterior_rate = posterior_rate,
    prior_weight = prior_weight,
    abs_prior_pull = abs_prior_pull,
    LCE = LCE,
    SDpost = SDpost
  )
}

# ------------------------------------------------------------------------------
# 4. Simulate nonlinear data-generating process
# ------------------------------------------------------------------------------

X <- matrix(
  rnorm(n * 5),
  nrow = n,
  ncol = 5
)

colnames(X) <- paste0(
  "x",
  1:5
)

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
# 5. Logistic-like model
# ------------------------------------------------------------------------------

model_logit <- glm(
  y ~ x1 + x2 + x3,
  data = dat,
  family = binomial()
)

dat$p_logit <-
  predict(
    model_logit,
    newdata = dat,
    type = "response"
  ) +
  rnorm(
    n,
    mean = 0,
    sd = 0.015
  )

dat$p_logit <- pmin(
  pmax(
    dat$p_logit,
    0.001
  ),
  0.5
)

# ------------------------------------------------------------------------------
# 6. RF-like model
# ------------------------------------------------------------------------------

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
# 7. Calculate LCE, SDpost and WNB for all prior concentrations
# ------------------------------------------------------------------------------

calculate_model_results <- function(
    y,
    p,
    model_label
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
        prior_m_values,
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

              # WNB definition used in the paper
              WNB = NB / (
                1 +
                  LCE +
                  SDpost
              )
            )
        }
      )
    }
  )
}

results_all <- bind_rows(
  calculate_model_results(
    y = dat$y,
    p = dat$p_logit,
    model_label = "Logistic"
  ),
  calculate_model_results(
    y = dat$y,
    p = dat$p_rf,
    model_label = "RF"
  )
)

# ------------------------------------------------------------------------------
# 8. Difference in WNB relative to m = 10
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

results_complete <- results_all %>%
  left_join(
    reference_WNB,
    by = c(
      "model",
      "threshold"
    )
  ) %>%
  mutate(
    delta_WNB =
      WNB -
      WNB_reference,

    abs_delta_WNB =
      abs(
        delta_WNB
      ),

    prior_weight_percent =
      100 *
      prior_weight,

    abs_prior_pull_pp =
      100 *
      abs_prior_pull
  ) %>%
  select(
    model,
    threshold,
    prior_m,
    n_bin,
    empirical_rate,
    posterior_rate,
    prior_weight_percent,
    abs_prior_pull_pp,
    LCE,
    SDpost,
    WNB,
    delta_WNB,
    abs_delta_WNB
  ) %>%
  arrange(
    factor(
      model,
      levels = c(
        "Logistic",
        "RF"
      )
    ),
    threshold,
    prior_m
  )

# ------------------------------------------------------------------------------
# 9. Export complete sensitivity table
# ------------------------------------------------------------------------------

write.csv(
  results_complete,
  file.path(
    output_dir,
    "Simulation2_prior_sensitivity_complete_results.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Figure 1: LCE and posterior SD
# ------------------------------------------------------------------------------

plot_data <- results_complete %>%
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
# 11. Figure 2: WNB difference relative to m = 10
# ------------------------------------------------------------------------------

p_delta_WNB <- results_complete %>%
  filter(
    prior_m != reference_prior_m
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
