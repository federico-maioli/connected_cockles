# Fit the cockle SDMs as two-part (hurdle) models, each part its own sdmTMB
# fit with its own coefficients and spatial field:
#   occurrence   binomial (cloglog)       - where cockles occur; every station
#   biomass      generalised gamma (log)  - how much biomass where they occur;
#                                           grab stations with cockles only
# The two parts share no parameters, so this is a standard delta model
# (cloglog presence link) in which some surveys inform occurrence only. The
# generalised gamma (Gamma and lognormal as special cases) fits the grab
# biomass better than the Gamma (AIC -15, cleaner DHARMa residuals;
# R/exploration/check_biomass_family.R).
#
# The surveys differ in what they can tell us:
#   Stock (grab, 0.1 m2, 2021-2025, incl. HighRes) - full efficiency, reliable
#     biomass; the reference survey
#   Stock2018 (hydraulic suction dredge, 0.75 m wide tows of ~100-250 m) -
#     full efficiency, but over a ~1,000x larger area, so its positives are
#     mostly sparse catches far below what a grab can register (a single
#     cockle in a grab is already 10 ind/m2): occurrence only
#   KSKV (oyster dredge, 0.965 m wide, tracks of ~100 or 400 m) - presences
#     true, absences possibly false (variable, unknown efficiency): occurrence
#     only, with log(dredged area) as an offset so short and long tows differ
# Occurrence includes survey (grab = reference) and year as factors: on the
# cloglog scale the survey terms are each gear's detection relative to one
# grab sample (P(detect) = 1 - exp(-q a d)). Biomass where present has year
# but no survey term (grab only), and previous-season fishing effort in every
# model: the swept-area ratio of the cockle fishery in the season leading into
# the survey year (02_clean_cockles.R), as log(1 + SAR), a control for
# depletion of fished beds. It is left out of occurrence, where it mainly
# reflects fishers targeting the beds.
#
# MAIN - in-strength is the connectivity metric:
#   space            survey + year (occurrence) / year + fishing (biomass)
#                    + spatial field
#   space_env        + environment
#   space_conn       + in-strength
#   space_env_conn   + environment + in-strength (full model)
# env = depth + depth^2 + temp + sal + oxy + shear_max, plus sediment type
# (Folk 7: mud, muddy sand, sand, mixed) in occurrence only: in the grab
# positives, sand has one station.
# Depth enters as a quadratic: occurrence peaks at intermediate depth.
#
# SENSITIVITY - the full model with in-strength replaced by each of the other
# four connectivity metrics (in-degree, in-closeness, eigenvector centrality,
# transitivity).
#
# Continuous covariates are z-scored over all stations; all models of a part
# share the same complete-case dataset, so their AIC values are comparable.
#
# Saved: data/sdm/main/<model>.rds and data/sdm/sensitivity/<model>.rds, each
# a list(occurrence, biomass) of fitted sdmTMB objects, plus one
# model_comparison.rds per folder (convergence, AIC of each part and total).

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
  mutate(
    across(all_of(c(env_vars, conn_vars)), ~ as.numeric(scale(.x)), .names = "{.col}_std"),
    # fishing effort is strongly skewed (most stations unfished)
    sar_prev_std = as.numeric(scale(log1p(sar_prev)))
  )

# 03 Model datasets ----
# keep rows complete on every covariate used by any model below, so every
# fit of a part shares the same n and AIC values are comparable
used_vars <- c(paste0(env_vars, "_std"), paste0(conn_vars, "_std"), "sediment", "sar_prev_std")
dat <- dat |>
  filter(
    !is.na(biomass), !is.na(present), !is.na(x_utm), !is.na(y_utm),
    if_all(all_of(used_vars), ~ !is.na(.x))
  ) |>
  # KSKV detection scales with the dredged area; other surveys have a fixed
  # sample size, absorbed by their survey term
  mutate(offset_area = if_else(survey == "KSKV", log(tow_area), 0))
stopifnot(!anyNA(dat$offset_area))

# biomass where present: grab stations with cockles (year levels without grab
# data dropped)
dat_pos <- dat |>
  filter(present == 1, survey == "Stock") |>
  mutate(year = droplevels(year))

count(dat, survey, present) |>
  print()
count(dat, sediment, present) |>
  print()

# 04 Barrier mesh ----
# variable-resolution mesh with the coastline built in, from 03_build_mesh.R:
# fine in densely surveyed beds, coarse on the sparse grid and on land; land
# triangles are a barrier so correlation does not cross headlands. The
# biomass part uses the same mesh with its own data points
mesh <- readRDS(here("data", "mesh", "mesh.rds"))
barrier_mesh <- sdmTMBextra::add_barrier_mesh(
  make_mesh(dat, c("x_utm", "y_utm"), mesh = mesh$mesh),
  mesh$land_barrier,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)
barrier_mesh_pos <- sdmTMBextra::add_barrier_mesh(
  make_mesh(dat_pos, c("x_utm", "y_utm"), mesh = mesh$mesh),
  mesh$land_barrier,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)

# 05 Fitting and bookkeeping ----
# the two parts with everything but the formula fixed; every model's formulas
# are written out explicitly at the call site
fit_occurrence <- function(formula) {
  sdmTMB(
    formula,
    data = dat, mesh = barrier_mesh, spatial = "on",
    family = binomial(link = "cloglog"), offset = dat$offset_area
  )
}
fit_biomass <- function(formula) {
  sdmTMB(
    formula,
    data = dat_pos, mesh = barrier_mesh_pos, spatial = "on",
    family = gengamma(link = "log")
  )
}

# check convergence, get the AIC of each part (the total is their sum: the
# parts share no parameters), save both fits
dir.create(here("data", "sdm", "main"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("data", "sdm", "sensitivity"), showWarnings = FALSE, recursive = TRUE)
record <- function(occurrence, biomass, model, folder) {
  ok <- function(f) isTRUE(suppressMessages(sanity(f, silent = TRUE))$all_ok)
  converged <- ok(occurrence) && ok(biomass)
  saveRDS(list(occurrence = occurrence, biomass = biomass), here("data", "sdm", folder, paste0(model, ".rds")))
  tibble(
    model = model, converged = converged,
    aic_presence = if (converged) AIC(occurrence) else NA_real_,
    aic_biomass = if (converged) AIC(biomass) else NA_real_,
    aic = aic_presence + aic_biomass
  )
}

comparison_main <- tibble()
comparison_sensitivity <- tibble()

# 06 Main models: in-strength ----
comparison_main <- bind_rows(comparison_main, record(
  fit_occurrence(present ~ survey + year),
  fit_biomass(biomass ~ year + sar_prev_std),
  "space", "main"
))

comparison_main <- bind_rows(comparison_main, record(
  fit_occurrence(present ~ survey + year + sediment + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std),
  fit_biomass(biomass ~ year + sar_prev_std + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std),
  "space_env", "main"
))

comparison_main <- bind_rows(comparison_main, record(
  fit_occurrence(present ~ survey + year + conn_in_strength_std),
  fit_biomass(biomass ~ year + sar_prev_std + conn_in_strength_std),
  "space_conn", "main"
))

comparison_main <- bind_rows(comparison_main, record(
  fit_occurrence(present ~ survey + year + sediment + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std),
  fit_biomass(biomass ~ year + sar_prev_std + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std),
  "space_env_conn", "main"
))

print(comparison_main |> mutate(delta_aic = aic - min(aic, na.rm = TRUE)) |> arrange(aic))

# 07 Sensitivity: the full model with the other connectivity metrics ----
comparison_sensitivity <- bind_rows(comparison_sensitivity, record(
  fit_occurrence(present ~ survey + year + sediment + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_degree_std),
  fit_biomass(biomass ~ year + sar_prev_std + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_degree_std),
  "space_env_conn_in_degree", "sensitivity"
))

comparison_sensitivity <- bind_rows(comparison_sensitivity, record(
  fit_occurrence(present ~ survey + year + sediment + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_closeness_std),
  fit_biomass(biomass ~ year + sar_prev_std + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_in_closeness_std),
  "space_env_conn_in_closeness", "sensitivity"
))

comparison_sensitivity <- bind_rows(comparison_sensitivity, record(
  fit_occurrence(present ~ survey + year + sediment + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_eigen_std),
  fit_biomass(biomass ~ year + sar_prev_std + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_eigen_std),
  "space_env_conn_eigen", "sensitivity"
))

comparison_sensitivity <- bind_rows(comparison_sensitivity, record(
  fit_occurrence(present ~ survey + year + sediment + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_transitivity_std),
  fit_biomass(biomass ~ year + sar_prev_std + depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std + conn_transitivity_std),
  "space_env_conn_transitivity", "sensitivity"
))

print(comparison_sensitivity)

# 08 Save comparison tables ----
saveRDS(comparison_main, here("data", "sdm", "main", "model_comparison.rds"))
saveRDS(comparison_sensitivity, here("data", "sdm", "sensitivity", "model_comparison.rds"))
