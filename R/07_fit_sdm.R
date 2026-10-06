# Fit the cockle SDMs as delta-gamma (hurdle) models. A delta-gamma model has
# two components fitted together:
#   1 presence     binomial (logit)  - where cockles occur
#   2 biomass      gamma (log)       - how much biomass, where they occur
# each with its own coefficients and its own spatial field.
#
# Every model includes survey (Stock grab 2021-2025 incl. HighRes = reference,
# Stock2018 suction dredge, KSKV dredge) and year (2018 reference) as factors
# and a spatial random field, in both components: the surveys use different
# gear and protocols, and occurrence and biomass differ between years.
#
# MAIN - in-strength is the connectivity metric:
#   space            survey + year + spatial field
#   space_env        + environment
#   space_conn       + in-strength
#   space_env_conn   + environment + in-strength (full model)
# env = depth + depth^2 + temp + oxy + sal + shear_max (depth enters as a
# quadratic: presence peaks at intermediate depth).
#
# SENSITIVITY - the full model with in-strength replaced by each of the other
# four connectivity metrics (in-degree, in-closeness, eigenvector centrality,
# transitivity).
#
# Continuous covariates are z-scored; all models share the same complete-case
# dataset so their AIC values are comparable.
#
# Saved: data/sdm/main/<model>.rds and data/sdm/sensitivity/<model>.rds (the
# fitted sdmTMB objects), plus one model_comparison.rds per folder
# (convergence, AIC of the delta model and of each part).

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)

source(here("R", "helpers.R")) # fmesher compatibility shim, needed by add_barrier_mesh() below

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

# 03 Shared model dataset ----
# keep rows complete on every covariate used by any model below, so every
# fit shares the same n and AIC values are comparable across all of them
used_vars <- c(paste0(env_vars, "_std"), paste0(conn_vars, "_std"))
dat <- dat |>
  filter(
    !is.na(biomass), !is.na(present), !is.na(x_utm), !is.na(y_utm),
    if_all(all_of(used_vars), ~ !is.na(.x))
  )

# 04 Barrier mesh ----
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

# data and mesh for fitting the two parts separately (see record() below):
# presence is biomass > 0, as in the delta model (the `present` column, from
# density, differs at two stations); biomass where present uses only the
# stations with cockles, on the same mesh and barrier
dat <- mutate(dat, present_biomass = as.integer(biomass > 0))
dat_pos <- filter(dat, biomass > 0)
barrier_mesh_pos <- sdmTMBextra::add_barrier_mesh(
  make_mesh(dat_pos, c("x_utm", "y_utm"), mesh = mesh$mesh),
  mesh$land_barrier,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)

# 05 Bookkeeping ----
# not a model-fitting helper - just avoids retyping "check convergence, get
# the AIC, save the fit" after every sdmTMB() call below. Every model's
# formula is written out explicitly at the call site.
dir.create(here("data", "sdm", "main"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("data", "sdm", "sensitivity"), showWarnings = FALSE, recursive = TRUE)

# AIC is also reported for each part. The two parts of a delta-gamma model
# share no parameters, so they are refitted separately with the same
# predictors: presence (binomial, biomass > 0, all stations) and biomass where
# present (gamma, stations with cockles). Their log-likelihoods sum exactly to
# the delta model's, so the two AICs add up to its AIC
record <- function(fit, model, folder) {
  rhs <- delete.response(terms(fit$formula[[1]]))
  fit_presence <- sdmTMB(update(rhs, present_biomass ~ .), data = dat, mesh = barrier_mesh, spatial = "on", family = binomial())
  fit_biomass <- sdmTMB(update(rhs, biomass ~ .), data = dat_pos, mesh = barrier_mesh_pos, spatial = "on", family = Gamma(link = "log"))
  ok <- function(f) isTRUE(suppressMessages(sanity(f, silent = TRUE))$all_ok)
  converged <- ok(fit) && ok(fit_presence) && ok(fit_biomass)
  saveRDS(fit, here("data", "sdm", folder, paste0(model, ".rds")))
  tibble(
    model = model, converged = converged,
    aic = if (converged) AIC(fit) else NA_real_,
    aic_presence = if (converged) AIC(fit_presence) else NA_real_,
    aic_biomass = if (converged) AIC(fit_biomass) else NA_real_
  )
}

comparison_main <- tibble()
comparison_sensitivity <- tibble()

# 06 Main models: in-strength ----
fit_space <- sdmTMB(
  biomass ~ survey + year,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_main <- bind_rows(comparison_main, record(fit_space, "space", "main"))

fit_space_env <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_main <- bind_rows(comparison_main, record(fit_space_env, "space_env", "main"))

fit_space_conn <- sdmTMB(
  biomass ~ survey + year + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_main <- bind_rows(comparison_main, record(fit_space_conn, "space_conn", "main"))

fit_space_env_conn <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_main <- bind_rows(comparison_main, record(fit_space_env_conn, "space_env_conn", "main"))

print(comparison_main |> mutate(delta_aic = aic - min(aic, na.rm = TRUE)) |> arrange(aic))

# 07 Sensitivity: the full model with the other connectivity metrics ----
fit_space_env_conn_in_degree <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_degree_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_in_degree, "space_env_conn_in_degree", "sensitivity"))

fit_space_env_conn_in_closeness <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_closeness_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_in_closeness, "space_env_conn_in_closeness", "sensitivity"))

fit_space_env_conn_eigen <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_eigen_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_eigen, "space_env_conn_eigen", "sensitivity"))

fit_space_env_conn_transitivity <- sdmTMB(
  biomass ~ survey + year + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_transitivity_std,
  data = dat, mesh = barrier_mesh, spatial = "on", family = delta_gamma()
)
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(fit_space_env_conn_transitivity, "space_env_conn_transitivity", "sensitivity"))

print(comparison_sensitivity)

# 08 Save comparison tables ----
saveRDS(comparison_main, here("data", "sdm", "main", "model_comparison.rds"))
saveRDS(comparison_sensitivity, here("data", "sdm", "sensitivity", "model_comparison.rds"))
