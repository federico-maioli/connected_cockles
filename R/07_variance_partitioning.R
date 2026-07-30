# Variance partitioning for the cockle biomass (Tweedie) SDMs.
#
# Nakagawa & Schielzeth (2013) marginal / conditional R2, all on the link scale:
#   var_f  variance of the fixed-effects linear predictor
#   var_l  variance of the spatial random field (0 when the field is off)
#   var_d  Tweedie observation-level (distribution) variance
# marginal R2    = var_f / (var_f + var_l + var_d)             fixed effects
# conditional R2 = (var_f + var_l) / (var_f + var_l + var_d)   fixed + field
# Total variance is split into: spatial field, environment, survey, connectivity,
# unexplained. The fixed (marginal) share is divided among the three fixed blocks
# by each block's covariance with the fixed predictor, so the blocks sum exactly
# to the marginal R2.
#
# Run for every connectivity metric, twice:
#   space_env_conn  with the spatial field
#   env_conn        without the spatial field
# Without the field, whatever the field was absorbing has to go somewhere, so the
# comparison shows how much of each block's share depends on it being there.
#
# Fits are reused from 06 so the partitioning matches the AIC table and the
# coefficient figures exactly (same data, same mesh, same estimates).
#
# Writes data/derived/variance_partition.rds; the figures and the LaTeX table are
# built from it by fig3_variance.R and tab2_variance.R.

library(tidyverse)
library(here)
library(sdmTMB)

env_cols <- c("depth_std", "temp_std", "oxy_std", "sal_std", "shear_max_std")

# split one fitted model's total variance into components
partition_fit <- function(fit, conn_term) {
  # with a field, predict() splits est into est_non_rf + est_rf; with the field
  # off it returns est only, which is already the fixed-effects predictor
  pr <- predict(fit)
  var_f <- if (is.null(pr$est_non_rf)) var(pr$est) else var(pr$est_non_rf)
  var_l <- if (is.null(pr$est_rf)) 0 else var(pr$est_rf)
  rp <- tidy(fit, effects = "ran_pars")
  phi <- rp$estimate[rp$term == "phi"]
  p_tw <- rp$estimate[rp$term == "tweedie_p"]
  var_d <- mean(log1p(phi * exp(pr$est)^(p_tw - 2))) # Tweedie observation variance
  tot <- var_f + var_l + var_d

  # each fixed block's share = cov(block, fixed predictor) / tot; the three sum
  # exactly to the marginal (fixed) R2 (block covariances absorbed by this split)
  coefs <- tidy(fit, effects = "fixed") |>
    select(term, estimate) |>
    deframe()
  X <- model.matrix(delete.response(terms(formula(fit))), data = fit$data)
  lp_term <- sweep(X, 2, coefs[colnames(X)], "*")
  b_env <- rowSums(lp_term[, env_cols, drop = FALSE])
  b_survey <- lp_term[, "surveyStock"]
  b_conn <- lp_term[, conn_term]
  eta_fixed <- b_env + b_survey + b_conn

  tibble(
    component = c("Spatial field", "Environment", "Survey (gear)", "Connectivity", "Unexplained"),
    share = c(
      var_l / tot,
      cov(b_env, eta_fixed) / tot,
      cov(b_survey, eta_fixed) / tot,
      cov(b_conn, eta_fixed) / tot,
      var_d / tot
    ),
    marginal_r2 = var_f / tot,
    conditional_r2 = (var_f + var_l) / tot
  )
}

# 01 Load fits ----
# only converged models are partitioned
fits <- readRDS(here("data", "derived", "sdm_fits.rds"))
ok_ids <- readRDS(here("data", "derived", "sdm_model_comparison.rds")) |>
  filter(converged) |>
  transmute(id = paste(response, model, sep = "_")) |>
  pull(id)

# order the metrics as they should appear downstream
metrics <- c(
  "log_biomass_in_strength", "presence_in_strength", "deg_in", "in_strength",
  "eigen_centrality", "closeness_centrality"
)

# 02 Partition every metric, with and without the spatial field ----
parts <- expand_grid(
  metric = metrics,
  structure = c("space_env_conn", "env_conn")
) |>
  mutate(id = paste("biomass", structure, metric, sep = "_")) |>
  filter(id %in% ok_ids) |>
  mutate(part = map2(id, metric, \(i, m) partition_fit(fits[[i]], paste0(m, "_std")))) |>
  unnest(part)

# 03 Report ----
parts |>
  filter(component == "Connectivity") |>
  select(metric, structure, connectivity_share = share, marginal_r2, conditional_r2) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(metric, structure) |>
  print(n = Inf)

# 04 Save ----
saveRDS(parts, here("data", "derived", "variance_partition.rds"))
