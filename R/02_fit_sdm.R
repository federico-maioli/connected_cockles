# Fit every SDM for cockle presence and biomass, and partition each one's
# variance. One script, three things happen for every model below:
#   1. fit it - the formula is written out in full at the call site, nothing
#      is assembled from a table
#   2. record it - check convergence, get AIC, and split its variance into
#      components (Nakagawa & Schielzeth R2; the fixed-effect share is further
#      split into Environment / Connectivity blocks HMSC-style - see
#      nakagawa_sdmtmb() in R/helpers.R for both)
#   3. save it - the fitted model and its variance table go to their own
#      files, so nothing needs re-fitting to look at either one later
#
# MAIN - in-strength is the connectivity metric (see figS5_space_confounding.R:
# it is the only one of the five whose coefficient survives dropping the
# spatial field). Six structures, crossing whether the spatial field is on and
# whether environment/connectivity are included:
#   space            spatial field only
#   space_env        env + spatial field
#   env              env                     (no field)
#   space_conn       connectivity + spatial field
#   space_env_conn   env + connectivity + spatial field
#   env_conn         env + connectivity      (no field)
# env = depth + temp + oxy + sal + shear_max.
#
# SENSITIVITY - two separate checks:
#   - the three connectivity structures (space_conn, space_env_conn, env_conn)
#     repeated for the other four metrics (in-degree, in-closeness,
#     eigenvector centrality, transitivity)
#   - the full model (space_env_conn / env_conn, in-strength), with one
#     environment term dropped at a time, spatial field on and off
#
# Every model is fitted for both responses:
#   present ~ ...   binomial (logit)   presence / absence
#   biomass ~ ...   Tweedie (log)      biomass
# Only the Stock survey remains after 01_clean_cockles.R, so there is no
# survey term. Continuous covariates are z-scored; all models share the same
# complete-case dataset so their AIC values are comparable.
#
# Everything is saved per model, in data/sdm/main/ or data/sdm/sensitivity/:
#   <response>_<model>.rds            the fitted sdmTMB object
#   <response>_<model>_variance.rds   its nakagawa_sdmtmb() table
# plus one model_comparison.rds per folder (AIC, marginal/conditional R2).

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)

source(here("R", "helpers.R"))

# 01 Load data ----
dat <- readRDS(here("data", "derived", "cockles_clean.rds"))

# 02 Standardise covariates ----
env_vars <- c("depth", "temp", "sal", "oxy", "shear_max")
conn_vars <- c("conn_in_degree", "conn_in_strength", "conn_in_closeness", "conn_eigen", "conn_transitivity")
dat <- dat |>
  mutate(across(all_of(c(env_vars, conn_vars)), ~ as.numeric(scale(.x)), .names = "{.col}_std"))

# 03 Predictor blocks for variance partitioning ----
# named so nakagawa_sdmtmb() can report each block's own share of the
# fixed-effect variance; defined once here and reused below so the same
# column names are not retyped at every one of the ~56 record() calls.
# Models with no predictors (space) pass no blocks at all - there is nothing
# to split. Models with only one kind of predictor (space_env, env,
# space_conn) still get a named block, so their fixed-effect variance is
# labelled Environment or Connectivity rather than a generic "fixed".
env_cols <- c("depth_std", "temp_std", "sal_std", "oxy_std", "shear_max_std")

blocks_env <- list(Environment = env_cols)

blocks_conn_in_strength <- list(Connectivity = "conn_in_strength_std")
blocks_env_conn_in_strength <- list(Environment = env_cols, Connectivity = "conn_in_strength_std")

blocks_conn_in_degree <- list(Connectivity = "conn_in_degree_std")
blocks_env_conn_in_degree <- list(Environment = env_cols, Connectivity = "conn_in_degree_std")

blocks_conn_in_closeness <- list(Connectivity = "conn_in_closeness_std")
blocks_env_conn_in_closeness <- list(Environment = env_cols, Connectivity = "conn_in_closeness_std")

blocks_conn_eigen <- list(Connectivity = "conn_eigen_std")
blocks_env_conn_eigen <- list(Environment = env_cols, Connectivity = "conn_eigen_std")

blocks_conn_transitivity <- list(Connectivity = "conn_transitivity_std")
blocks_env_conn_transitivity <- list(Environment = env_cols, Connectivity = "conn_transitivity_std")

# leave-one-out: the full in-strength block list, minus one environment term
blocks_env_conn_drop_depth <- list(Environment = setdiff(env_cols, "depth_std"), Connectivity = "conn_in_strength_std")
blocks_env_conn_drop_temp <- list(Environment = setdiff(env_cols, "temp_std"), Connectivity = "conn_in_strength_std")
blocks_env_conn_drop_sal <- list(Environment = setdiff(env_cols, "sal_std"), Connectivity = "conn_in_strength_std")
blocks_env_conn_drop_oxy <- list(Environment = setdiff(env_cols, "oxy_std"), Connectivity = "conn_in_strength_std")
blocks_env_conn_drop_shear_max <- list(Environment = setdiff(env_cols, "shear_max_std"), Connectivity = "conn_in_strength_std")

# 04 Shared model dataset ----
# keep rows complete on every covariate used by any model below, so every
# fit shares the same n and AIC values are comparable across all of them
used_vars <- c(paste0(env_vars, "_std"), paste0(conn_vars, "_std"))
dat <- dat |>
  filter(
    !is.na(biomass), !is.na(present), !is.na(x_utm), !is.na(y_utm),
    if_all(all_of(used_vars), ~ !is.na(.x))
  )

# 05 Coastline barrier mesh ----
barrier_mesh <- build_barrier_mesh(dat)
cat(
  "mesh vertices:", barrier_mesh$mesh$n,
  "| water triangles:", length(barrier_mesh$normal_triangles),
  "| land triangles:", length(barrier_mesh$barrier_triangles), "\n"
)

# 06 Bookkeeping ----
# not a model-fitting helper - just avoids retyping "check convergence, get
# the AIC, partition the variance, save everything, add a row to the
# comparison table" after every one of the ~56 sdmTMB() calls below. Every
# model's formula and predictor blocks are still written out explicitly at
# the call site.
dir.create(here("data", "sdm", "main"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("data", "sdm", "sensitivity"), showWarnings = FALSE, recursive = TRUE)

record <- function(fit, response, model, spatial, folder, blocks = NULL) {
  converged <- isTRUE(suppressMessages(sanity(fit, silent = TRUE))$all_ok)
  saveRDS(fit, here("data", "sdm", folder, paste0(response, "_", model, ".rds")))

  variance <- if (converged) nakagawa_sdmtmb(fit, blocks = blocks) else NULL
  saveRDS(variance, here("data", "sdm", folder, paste0(response, "_", model, "_variance.rds")))

  tibble(
    response = response, model = model, spatial = spatial, converged = converged,
    aic = if (converged) AIC(fit) else NA_real_,
    r2_marginal = if (converged) variance$r2_marginal[1] else NA_real_,
    r2_conditional = if (converged) variance$r2_conditional[1] else NA_real_
  )
}

comparison_main <- tibble()
comparison_sensitivity <- tibble()

# 07 Main models: in-strength, biomass ----
fit_biomass_space <- sdmTMB(
  biomass ~ 1,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_main <- bind_rows(comparison_main, record(fit_biomass_space, "biomass", "space", "on", "main"))

fit_biomass_space_env <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_main <- bind_rows(comparison_main, record(fit_biomass_space_env, "biomass", "space_env", "on", "main", blocks_env))

fit_biomass_env <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_main <- bind_rows(comparison_main, record(fit_biomass_env, "biomass", "env", "off", "main", blocks_env))

fit_biomass_space_conn <- sdmTMB(
  biomass ~ conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_main <- bind_rows(comparison_main, record(fit_biomass_space_conn, "biomass", "space_conn", "on", "main", blocks_conn_in_strength))

fit_biomass_space_env_conn <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_main <- bind_rows(comparison_main, record(fit_biomass_space_env_conn, "biomass", "space_env_conn", "on", "main", blocks_env_conn_in_strength))

fit_biomass_env_conn <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_main <- bind_rows(comparison_main, record(fit_biomass_env_conn, "biomass", "env_conn", "off", "main", blocks_env_conn_in_strength))

# 08 Main models: in-strength, presence ----
fit_present_space <- sdmTMB(
  present ~ 1,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_main <- bind_rows(comparison_main, record(fit_present_space, "present", "space", "on", "main"))

fit_present_space_env <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_main <- bind_rows(comparison_main, record(fit_present_space_env, "present", "space_env", "on", "main", blocks_env))

fit_present_env <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_main <- bind_rows(comparison_main, record(fit_present_env, "present", "env", "off", "main", blocks_env))

fit_present_space_conn <- sdmTMB(
  present ~ conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_main <- bind_rows(comparison_main, record(fit_present_space_conn, "present", "space_conn", "on", "main", blocks_conn_in_strength))

fit_present_space_env_conn <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_main <- bind_rows(comparison_main, record(fit_present_space_env_conn, "present", "space_env_conn", "on", "main", blocks_env_conn_in_strength))

fit_present_env_conn <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_main <- bind_rows(comparison_main, record(fit_present_env_conn, "present", "env_conn", "off", "main", blocks_env_conn_in_strength))

cat("\nMain models (in-strength):\n")
print(comparison_main |> mutate(delta_aic = aic - min(aic, na.rm = TRUE), .by = response) |> arrange(response, aic), n = Inf)

# 09 Sensitivity: other connectivity metrics, biomass ----
# same three connectivity structures as the main model, repeated once per
# metric (in-strength itself is already fitted above, not repeated here)

# in-degree
fit_biomass_space_conn_in_degree <- sdmTMB(
  biomass ~ conn_in_degree_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_conn_in_degree, "biomass", "space_conn_in_degree", "on", "sensitivity", blocks_conn_in_degree))

fit_biomass_space_env_conn_in_degree <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_degree_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_env_conn_in_degree, "biomass", "space_env_conn_in_degree", "on", "sensitivity", blocks_env_conn_in_degree))

fit_biomass_env_conn_in_degree <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_degree_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_env_conn_in_degree, "biomass", "env_conn_in_degree", "off", "sensitivity", blocks_env_conn_in_degree))

# in-closeness
fit_biomass_space_conn_in_closeness <- sdmTMB(
  biomass ~ conn_in_closeness_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_conn_in_closeness, "biomass", "space_conn_in_closeness", "on", "sensitivity", blocks_conn_in_closeness))

fit_biomass_space_env_conn_in_closeness <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_closeness_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_env_conn_in_closeness, "biomass", "space_env_conn_in_closeness", "on", "sensitivity", blocks_env_conn_in_closeness))

fit_biomass_env_conn_in_closeness <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_closeness_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_env_conn_in_closeness, "biomass", "env_conn_in_closeness", "off", "sensitivity", blocks_env_conn_in_closeness))

# eigenvector centrality
fit_biomass_space_conn_eigen <- sdmTMB(
  biomass ~ conn_eigen_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_conn_eigen, "biomass", "space_conn_eigen", "on", "sensitivity", blocks_conn_eigen))

fit_biomass_space_env_conn_eigen <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_eigen_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_env_conn_eigen, "biomass", "space_env_conn_eigen", "on", "sensitivity", blocks_env_conn_eigen))

fit_biomass_env_conn_eigen <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_eigen_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_env_conn_eigen, "biomass", "env_conn_eigen", "off", "sensitivity", blocks_env_conn_eigen))

# transitivity
fit_biomass_space_conn_transitivity <- sdmTMB(
  biomass ~ conn_transitivity_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_conn_transitivity, "biomass", "space_conn_transitivity", "on", "sensitivity", blocks_conn_transitivity))

fit_biomass_space_env_conn_transitivity <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_transitivity_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_env_conn_transitivity, "biomass", "space_env_conn_transitivity", "on", "sensitivity", blocks_env_conn_transitivity))

fit_biomass_env_conn_transitivity <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_transitivity_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_env_conn_transitivity, "biomass", "env_conn_transitivity", "off", "sensitivity", blocks_env_conn_transitivity))

# 10 Sensitivity: other connectivity metrics, presence ----

# in-degree
fit_present_space_conn_in_degree <- sdmTMB(
  present ~ conn_in_degree_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_conn_in_degree, "present", "space_conn_in_degree", "on", "sensitivity", blocks_conn_in_degree))

fit_present_space_env_conn_in_degree <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_degree_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_env_conn_in_degree, "present", "space_env_conn_in_degree", "on", "sensitivity", blocks_env_conn_in_degree))

fit_present_env_conn_in_degree <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_degree_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_env_conn_in_degree, "present", "env_conn_in_degree", "off", "sensitivity", blocks_env_conn_in_degree))

# in-closeness
fit_present_space_conn_in_closeness <- sdmTMB(
  present ~ conn_in_closeness_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_conn_in_closeness, "present", "space_conn_in_closeness", "on", "sensitivity", blocks_conn_in_closeness))

fit_present_space_env_conn_in_closeness <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_closeness_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_env_conn_in_closeness, "present", "space_env_conn_in_closeness", "on", "sensitivity", blocks_env_conn_in_closeness))

fit_present_env_conn_in_closeness <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_closeness_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_env_conn_in_closeness, "present", "env_conn_in_closeness", "off", "sensitivity", blocks_env_conn_in_closeness))

# eigenvector centrality
fit_present_space_conn_eigen <- sdmTMB(
  present ~ conn_eigen_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_conn_eigen, "present", "space_conn_eigen", "on", "sensitivity", blocks_conn_eigen))

fit_present_space_env_conn_eigen <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_eigen_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_env_conn_eigen, "present", "space_env_conn_eigen", "on", "sensitivity", blocks_env_conn_eigen))

fit_present_env_conn_eigen <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_eigen_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_env_conn_eigen, "present", "env_conn_eigen", "off", "sensitivity", blocks_env_conn_eigen))

# transitivity
fit_present_space_conn_transitivity <- sdmTMB(
  present ~ conn_transitivity_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_conn_transitivity, "present", "space_conn_transitivity", "on", "sensitivity", blocks_conn_transitivity))

fit_present_space_env_conn_transitivity <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_transitivity_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_env_conn_transitivity, "present", "space_env_conn_transitivity", "on", "sensitivity", blocks_env_conn_transitivity))

fit_present_env_conn_transitivity <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_transitivity_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_env_conn_transitivity, "present", "env_conn_transitivity", "off", "sensitivity", blocks_env_conn_transitivity))

# 11 Sensitivity: environment leave-one-out, biomass ----
# the full model (env + in-strength), minus one environment term at a time,
# spatial field on and off - does the connectivity coefficient depend on any
# single environment covariate?

# drop depth
fit_biomass_space_env_conn_drop_depth <- sdmTMB(
  biomass ~ temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_env_conn_drop_depth, "biomass", "space_env_conn_drop_depth", "on", "sensitivity", blocks_env_conn_drop_depth))

fit_biomass_env_conn_drop_depth <- sdmTMB(
  biomass ~ temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_env_conn_drop_depth, "biomass", "env_conn_drop_depth", "off", "sensitivity", blocks_env_conn_drop_depth))

# drop temperature
fit_biomass_space_env_conn_drop_temp <- sdmTMB(
  biomass ~ depth_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_env_conn_drop_temp, "biomass", "space_env_conn_drop_temp", "on", "sensitivity", blocks_env_conn_drop_temp))

fit_biomass_env_conn_drop_temp <- sdmTMB(
  biomass ~ depth_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_env_conn_drop_temp, "biomass", "env_conn_drop_temp", "off", "sensitivity", blocks_env_conn_drop_temp))

# drop salinity
fit_biomass_space_env_conn_drop_sal <- sdmTMB(
  biomass ~ depth_std + temp_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_env_conn_drop_sal, "biomass", "space_env_conn_drop_sal", "on", "sensitivity", blocks_env_conn_drop_sal))

fit_biomass_env_conn_drop_sal <- sdmTMB(
  biomass ~ depth_std + temp_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_env_conn_drop_sal, "biomass", "env_conn_drop_sal", "off", "sensitivity", blocks_env_conn_drop_sal))

# drop oxygen
fit_biomass_space_env_conn_drop_oxy <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_env_conn_drop_oxy, "biomass", "space_env_conn_drop_oxy", "on", "sensitivity", blocks_env_conn_drop_oxy))

fit_biomass_env_conn_drop_oxy <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_env_conn_drop_oxy, "biomass", "env_conn_drop_oxy", "off", "sensitivity", blocks_env_conn_drop_oxy))

# drop shear stress
fit_biomass_space_env_conn_drop_shear_max <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_space_env_conn_drop_shear_max, "biomass", "space_env_conn_drop_shear_max", "on", "sensitivity", blocks_env_conn_drop_shear_max))

fit_biomass_env_conn_drop_shear_max <- sdmTMB(
  biomass ~ depth_std + temp_std + sal_std + oxy_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = tweedie(link = "log")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_biomass_env_conn_drop_shear_max, "biomass", "env_conn_drop_shear_max", "off", "sensitivity", blocks_env_conn_drop_shear_max))

# 12 Sensitivity: environment leave-one-out, presence ----

# drop depth
fit_present_space_env_conn_drop_depth <- sdmTMB(
  present ~ temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_env_conn_drop_depth, "present", "space_env_conn_drop_depth", "on", "sensitivity", blocks_env_conn_drop_depth))

fit_present_env_conn_drop_depth <- sdmTMB(
  present ~ temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_env_conn_drop_depth, "present", "env_conn_drop_depth", "off", "sensitivity", blocks_env_conn_drop_depth))

# drop temperature
fit_present_space_env_conn_drop_temp <- sdmTMB(
  present ~ depth_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_env_conn_drop_temp, "present", "space_env_conn_drop_temp", "on", "sensitivity", blocks_env_conn_drop_temp))

fit_present_env_conn_drop_temp <- sdmTMB(
  present ~ depth_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_env_conn_drop_temp, "present", "env_conn_drop_temp", "off", "sensitivity", blocks_env_conn_drop_temp))

# drop salinity
fit_present_space_env_conn_drop_sal <- sdmTMB(
  present ~ depth_std + temp_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_env_conn_drop_sal, "present", "space_env_conn_drop_sal", "on", "sensitivity", blocks_env_conn_drop_sal))

fit_present_env_conn_drop_sal <- sdmTMB(
  present ~ depth_std + temp_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_env_conn_drop_sal, "present", "env_conn_drop_sal", "off", "sensitivity", blocks_env_conn_drop_sal))

# drop oxygen
fit_present_space_env_conn_drop_oxy <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_env_conn_drop_oxy, "present", "space_env_conn_drop_oxy", "on", "sensitivity", blocks_env_conn_drop_oxy))

fit_present_env_conn_drop_oxy <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_env_conn_drop_oxy, "present", "env_conn_drop_oxy", "off", "sensitivity", blocks_env_conn_drop_oxy))

# drop shear stress
fit_present_space_env_conn_drop_shear_max <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_space_env_conn_drop_shear_max, "present", "space_env_conn_drop_shear_max", "on", "sensitivity", blocks_env_conn_drop_shear_max))

fit_present_env_conn_drop_shear_max <- sdmTMB(
  present ~ depth_std + temp_std + sal_std + oxy_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "off", family = binomial(link = "logit")
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_present_env_conn_drop_shear_max, "present", "env_conn_drop_shear_max", "off", "sensitivity", blocks_env_conn_drop_shear_max))

cat("\nSensitivity models (other metrics + environment leave-one-out):\n")
print(comparison_sensitivity, n = Inf)

# 13 Save comparison tables ----
saveRDS(comparison_main, here("data", "sdm", "main", "model_comparison.rds"))
saveRDS(comparison_sensitivity, here("data", "sdm", "sensitivity", "model_comparison.rds"))
