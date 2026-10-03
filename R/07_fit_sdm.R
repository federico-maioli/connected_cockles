# Fit every SDM for cockle biomass as a delta-gamma (hurdle) model, and
# record each one. A delta-gamma model has two components fitted together:
#   1 presence     binomial (logit)  - where cockles occur
#   2 biomass      gamma (log)       - how much biomass, where they occur
# each with its own coefficients and its own spatial field, so presence and
# conditional biomass come from one model and no separate presence-absence
# model is needed. (Chosen over Tweedie and delta-lognormal by AIC, see
# R/exploration/compare_delta_families.R.)
#
# For every model below:
#   1. fit it - the formula is written out in full at the call site, nothing
#      is assembled from a table
#   2. record it - check convergence, get AIC, and split the variance of each
#      component separately (Nakagawa & Schielzeth R2; the fixed-effect share
#      further split into Environment / Connectivity blocks HMSC-style - see
#      nakagawa_sdmtmb() in R/helpers.R). The biomass component is partitioned
#      over the positive observations only, the ones it describes.
#   3. save it - the fitted model and its variance table go to their own
#      files, so nothing needs re-fitting to look at either one later
#
# MAIN - in-strength is the connectivity metric. Six structures, crossing
# whether the spatial field is on and whether environment/connectivity are
# included:
#   space            spatial field only
#   space_env        env + spatial field
#   env              env                     (no field)
#   space_conn       connectivity + spatial field
#   space_env_conn   env + connectivity + spatial field
#   env_conn         env + connectivity      (no field)
# env = depth + depth^2 + temp + oxy + sal + shear_max (depth enters as a
# quadratic: presence peaks at intermediate depth, see
# R/exploration/check_quadratic_terms.R).
#
# SENSITIVITY - two separate checks:
#   - the three connectivity structures (space_conn, space_env_conn, env_conn)
#     repeated for the other four metrics (in-degree, in-closeness,
#     eigenvector centrality, transitivity)
#   - the full model (space_env_conn / env_conn, in-strength), with one
#     environment term dropped at a time, spatial field on and off
#
# Every model includes survey (Stock, Stock2018, KSKV) and year as factors in
# both components: the surveys use different gear and protocols, and
# occurrence and biomass differ between years. Continuous covariates are z-scored; all models share the same
# complete-case dataset so their AIC values are comparable.
#
# Saved per model, in data/sdm/main/ or data/sdm/sensitivity/:
#   <model>.rds            the fitted sdmTMB object
#   <model>_variance.rds   its nakagawa_sdmtmb() table, one block per component
#                          (part = "presence" / "biomass")
# plus one model_comparison.rds per folder (AIC, marginal/conditional R2 of
# each component).

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)

source(here("R", "helpers.R"))

# 01 Load data ----
dat <- readRDS(here("data", "cockles", "derived", "cockles_env_conn.rds"))

# connectivity metrics used throughout: presence-weighted (source rows of the
# flow matrix weighted by the suitability model's presence probability, see
# 05_weight_connectivity.R). Switching weighting only means changing these lines
dat <- dat |>
  mutate(
    conn_in_degree = conn_in_degree_presence,
    conn_in_strength = conn_in_strength_presence,
    conn_in_closeness = conn_in_closeness_presence,
    conn_eigen = conn_eigen_presence,
    conn_transitivity = conn_transitivity_presence
  )

# 02 Standardise covariates ----
env_vars <- c("depth", "temp", "sal", "oxy", "shear_max")
conn_vars <- c("conn_in_degree", "conn_in_strength", "conn_in_closeness", "conn_eigen", "conn_transitivity")
dat <- dat |>
  mutate(across(all_of(c(env_vars, conn_vars)), ~ as.numeric(scale(.x)), .names = "{.col}_std"))

# 03 Predictor blocks for variance partitioning ----
# named so nakagawa_sdmtmb() can report each block's own share of the
# fixed-effect variance; defined once here and reused below so the same
# column names are not retyped at every record() call.
# Models with no predictors (space) pass no blocks at all - there is nothing
# to split. Models with only one kind of predictor (space_env, env,
# space_conn) still get a named block, so their fixed-effect variance is
# labelled Environment or Connectivity rather than a generic "fixed".
env_cols <- c("depth_std", "I(depth_std^2)", "temp_std", "sal_std", "oxy_std", "shear_max_std")
# survey (Stock grab 2021-2025 incl. HighRes = reference, Stock2018 suction
# dredge, KSKV dredge) enters every model, both components:
# the surveys use different gear and protocols, so catchability differs
survey_cols <- c("surveyStock2018", "surveyKSKV")
# year as a factor in every model too (2018 reference): annual differences in
# occurrence and biomass (recruitment, fishing), kept out of the spatial field
year_cols <- paste0("year", levels(dat$year)[-1])
blocks_survey <- list(Survey = survey_cols, Year = year_cols)

blocks_env <- list(Survey = survey_cols, Year = year_cols, Environment = env_cols)

blocks_conn_in_strength <- list(Survey = survey_cols, Year = year_cols, Connectivity = "conn_in_strength_std")
blocks_env_conn_in_strength <- list(Survey = survey_cols, Year = year_cols, Environment = env_cols, Connectivity = "conn_in_strength_std")

blocks_conn_in_degree <- list(Survey = survey_cols, Year = year_cols, Connectivity = "conn_in_degree_std")
blocks_env_conn_in_degree <- list(Survey = survey_cols, Year = year_cols, Environment = env_cols, Connectivity = "conn_in_degree_std")

blocks_conn_in_closeness <- list(Survey = survey_cols, Year = year_cols, Connectivity = "conn_in_closeness_std")
blocks_env_conn_in_closeness <- list(Survey = survey_cols, Year = year_cols, Environment = env_cols, Connectivity = "conn_in_closeness_std")

blocks_conn_eigen <- list(Survey = survey_cols, Year = year_cols, Connectivity = "conn_eigen_std")
blocks_env_conn_eigen <- list(Survey = survey_cols, Year = year_cols, Environment = env_cols, Connectivity = "conn_eigen_std")

blocks_conn_transitivity <- list(Survey = survey_cols, Year = year_cols, Connectivity = "conn_transitivity_std")
blocks_env_conn_transitivity <- list(Survey = survey_cols, Year = year_cols, Environment = env_cols, Connectivity = "conn_transitivity_std")

# leave-one-out: the full in-strength block list, minus one environment term
blocks_env_conn_drop_depth <- list(Survey = survey_cols, Year = year_cols, Environment = setdiff(env_cols, c("depth_std", "I(depth_std^2)")), Connectivity = "conn_in_strength_std")
blocks_env_conn_drop_temp <- list(Survey = survey_cols, Year = year_cols, Environment = setdiff(env_cols, "temp_std"), Connectivity = "conn_in_strength_std")
blocks_env_conn_drop_sal <- list(Survey = survey_cols, Year = year_cols, Environment = setdiff(env_cols, "sal_std"), Connectivity = "conn_in_strength_std")
blocks_env_conn_drop_oxy <- list(Survey = survey_cols, Year = year_cols, Environment = setdiff(env_cols, "oxy_std"), Connectivity = "conn_in_strength_std")
blocks_env_conn_drop_shear_max <- list(Survey = survey_cols, Year = year_cols, Environment = setdiff(env_cols, "shear_max_std"), Connectivity = "conn_in_strength_std")

# 04 Shared model dataset ----
# keep rows complete on every covariate used by any model below, so every
# fit shares the same n and AIC values are comparable across all of them
used_vars <- c(paste0(env_vars, "_std"), paste0(conn_vars, "_std"))
dat <- dat |>
  filter(
    !is.na(biomass), !is.na(present), !is.na(x_utm), !is.na(y_utm),
    if_all(all_of(used_vars), ~ !is.na(.x))
  )

# 05 Barrier mesh ----
# variable-resolution mesh with the coastline built in, from 03_build_mesh.R:
# fine in densely surveyed beds, coarse on the sparse grid and on land; land
# triangles are a barrier so correlation does not cross headlands
mesh <- readRDS(here("data", "mesh", "mesh.rds"))
barrier_mesh <- sdmTMBextra::add_barrier_mesh(
  make_mesh(dat, c("x_utm", "y_utm"), mesh = mesh$mesh),
  mesh$land_barrier,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)

# 06 Bookkeeping ----
# not a model-fitting helper - just avoids retyping "check convergence, get
# the AIC, partition the variance, save everything, add a row to the
# comparison table" after every sdmTMB() call below. Every model's formula and
# predictor blocks are still written out explicitly at the call site.
dir.create(here("data", "sdm", "main"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("data", "sdm", "sensitivity"), showWarnings = FALSE, recursive = TRUE)

record <- function(fit, model, spatial, folder, blocks = NULL) {
  converged <- isTRUE(suppressMessages(sanity(fit, silent = TRUE))$all_ok)
  saveRDS(fit, here("data", "sdm", folder, paste0(model, ".rds")))

  # one variance table per delta component: 1 presence, 2 biomass where present
  variance <- if (converged) {
    bind_rows(
      nakagawa_sdmtmb(fit, blocks = blocks, model = 1) |> mutate(part = "presence"),
      nakagawa_sdmtmb(fit, blocks = blocks, model = 2) |> mutate(part = "biomass")
    )
  } else {
    NULL
  }
  saveRDS(variance, here("data", "sdm", folder, paste0(model, "_variance.rds")))

  r2 <- function(part, col) if (converged) variance[[col]][variance$part == part][1] else NA_real_
  tibble(
    model = model, spatial = spatial, converged = converged,
    aic = if (converged) AIC(fit) else NA_real_,
    r2_marginal_presence = r2("presence", "r2_marginal"),
    r2_conditional_presence = r2("presence", "r2_conditional"),
    r2_marginal_biomass = r2("biomass", "r2_marginal"),
    r2_conditional_biomass = r2("biomass", "r2_conditional")
  )
}

comparison_main <- tibble()
comparison_sensitivity <- tibble()

# 07 Main models: in-strength ----
fit_space <- sdmTMB(
  biomass ~ survey + year,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_main <- bind_rows(comparison_main, record(fit_space, "space", "on", "main", blocks_survey))

fit_space_env <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_main <- bind_rows(comparison_main, record(fit_space_env, "space_env", "on", "main", blocks_env))

fit_env <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_main <- bind_rows(comparison_main, record(fit_env, "env", "off", "main", blocks_env))

fit_space_conn <- sdmTMB(
  biomass ~ survey + year + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_main <- bind_rows(comparison_main, record(fit_space_conn, "space_conn", "on", "main", blocks_conn_in_strength))

fit_space_env_conn <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_main <- bind_rows(comparison_main, record(fit_space_env_conn, "space_env_conn", "on", "main", blocks_env_conn_in_strength))

fit_env_conn <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_main <- bind_rows(comparison_main, record(fit_env_conn, "env_conn", "off", "main", blocks_env_conn_in_strength))

print(comparison_main |> mutate(delta_aic = aic - min(aic, na.rm = TRUE)) |> arrange(aic), n = Inf)

# 08 Sensitivity: other connectivity metrics ----
# same three connectivity structures as the main model, repeated once per
# metric (in-strength itself is already fitted above, not repeated here)

# in-degree
fit_space_conn_in_degree <- sdmTMB(
  biomass ~ survey + year + conn_in_degree_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_conn_in_degree, "space_conn_in_degree", "on", "sensitivity", blocks_conn_in_degree))

fit_space_env_conn_in_degree <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_degree_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_in_degree, "space_env_conn_in_degree", "on", "sensitivity", blocks_env_conn_in_degree))

fit_env_conn_in_degree <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_degree_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_env_conn_in_degree, "env_conn_in_degree", "off", "sensitivity", blocks_env_conn_in_degree))

# in-closeness
fit_space_conn_in_closeness <- sdmTMB(
  biomass ~ survey + year + conn_in_closeness_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_conn_in_closeness, "space_conn_in_closeness", "on", "sensitivity", blocks_conn_in_closeness))

fit_space_env_conn_in_closeness <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_closeness_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_in_closeness, "space_env_conn_in_closeness", "on", "sensitivity", blocks_env_conn_in_closeness))

fit_env_conn_in_closeness <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_closeness_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_env_conn_in_closeness, "env_conn_in_closeness", "off", "sensitivity", blocks_env_conn_in_closeness))

# eigenvector centrality
fit_space_conn_eigen <- sdmTMB(
  biomass ~ survey + year + conn_eigen_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_conn_eigen, "space_conn_eigen", "on", "sensitivity", blocks_conn_eigen))

fit_space_env_conn_eigen <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_eigen_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_eigen, "space_env_conn_eigen", "on", "sensitivity", blocks_env_conn_eigen))

fit_env_conn_eigen <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_eigen_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_env_conn_eigen, "env_conn_eigen", "off", "sensitivity", blocks_env_conn_eigen))

# transitivity
fit_space_conn_transitivity <- sdmTMB(
  biomass ~ survey + year + conn_transitivity_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_conn_transitivity, "space_conn_transitivity", "on", "sensitivity", blocks_conn_transitivity))

fit_space_env_conn_transitivity <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_transitivity_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_transitivity, "space_env_conn_transitivity", "on", "sensitivity", blocks_env_conn_transitivity))

fit_env_conn_transitivity <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_transitivity_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_env_conn_transitivity, "env_conn_transitivity", "off", "sensitivity", blocks_env_conn_transitivity))

# 09 Sensitivity: environment leave-one-out ----
# the full model (env + in-strength), minus one environment term at a time
# (dropping depth drops depth^2 too),
# spatial field on and off - does the connectivity coefficient depend on any
# single environment covariate?

# drop depth
fit_space_env_conn_drop_depth <- sdmTMB(
  biomass ~ survey + year + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_drop_depth, "space_env_conn_drop_depth", "on", "sensitivity", blocks_env_conn_drop_depth))

fit_env_conn_drop_depth <- sdmTMB(
  biomass ~ survey + year + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_env_conn_drop_depth, "env_conn_drop_depth", "off", "sensitivity", blocks_env_conn_drop_depth))

# drop temperature
fit_space_env_conn_drop_temp <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_drop_temp, "space_env_conn_drop_temp", "on", "sensitivity", blocks_env_conn_drop_temp))

fit_env_conn_drop_temp <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_env_conn_drop_temp, "env_conn_drop_temp", "off", "sensitivity", blocks_env_conn_drop_temp))

# drop salinity
fit_space_env_conn_drop_sal <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_drop_sal, "space_env_conn_drop_sal", "on", "sensitivity", blocks_env_conn_drop_sal))

fit_env_conn_drop_sal <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_env_conn_drop_sal, "env_conn_drop_sal", "off", "sensitivity", blocks_env_conn_drop_sal))

# drop oxygen
fit_space_env_conn_drop_oxy <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_drop_oxy, "space_env_conn_drop_oxy", "on", "sensitivity", blocks_env_conn_drop_oxy))

fit_env_conn_drop_oxy <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_env_conn_drop_oxy, "env_conn_drop_oxy", "off", "sensitivity", blocks_env_conn_drop_oxy))

# drop shear stress
fit_space_env_conn_drop_shear_max <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_drop_shear_max, "space_env_conn_drop_shear_max", "on", "sensitivity", blocks_env_conn_drop_shear_max))

fit_env_conn_drop_shear_max <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_env_conn_drop_shear_max, "env_conn_drop_shear_max", "off", "sensitivity", blocks_env_conn_drop_shear_max))


print(comparison_sensitivity, n = Inf)

# 10 Save comparison tables ----
saveRDS(comparison_main, here("data", "sdm", "main", "model_comparison.rds"))
saveRDS(comparison_sensitivity, here("data", "sdm", "sensitivity", "model_comparison.rds"))
