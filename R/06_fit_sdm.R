# Fit and compare SDMs for cockle presence and biomass.
#
# Six model structures crossing space (spatial random field), environment, and
# connectivity, each fitted for both responses:
#   space            survey + spatial field
#   space_env        survey + env + spatial field
#   env              survey + env                     (no field)
#   space_conn       survey + connectivity + spatial field
#   space_env_conn   survey + env + connectivity + spatial field
#   env_conn         survey + env + connectivity      (no field)
# Every model carries the survey (gear) intercept; env = depth + temp + oxy +
# sal + shear_max. The connectivity slot is rotated over three predictors
# (eigenvector centrality, in-strength, log biomass supply), so each
# connectivity model is fitted three times.
#   present ~ ...   binomial (logit)   presence / absence
#   biomass ~ ...   Tweedie (log)      biomass
# Continuous covariates are z-scored here (model-input scaling); all models share
# the same complete-case dataset so their AIC values are comparable.

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
library(rnaturalearth)

# compatibility shim: sdmTMBextra 0.0.5 calls fmesher::fm_identical_CRS, renamed
# to fm_crs_is_identical in fmesher >= 0.7; register the old name as an alias.
# required for add_barrier_mesh() below - do not remove unless sdmTMBextra is updated
local({
  exps <- .getNamespaceInfo(asNamespace("fmesher"), "exports")
  if (!exists("fm_identical_CRS", envir = exps, inherits = FALSE)) {
    assign("fm_identical_CRS", "fm_crs_is_identical", envir = exps)
  }
})

# variance inflation factor for each predictor (1 / (1 - R^2) vs the others)
vif <- function(data, vars) {
  map_dbl(set_names(vars), function(v) {
    others <- setdiff(vars, v)
    r2 <- summary(lm(reformulate(others, v), data = data))$r.squared
    1 / (1 - r2)
  })
}

# ---- Model structures to compare -------------------------------------------
# `survey` is a factor (gear) carried in every model; env is the environmental
# block; the spatial field is toggled per model via `spatial`. Connectivity
# models rotate one of three connectivity predictors into the RHS.
survey_term <- "survey"
env_terms <- "depth_std + temp_std + oxy_std + sal_std + shear_max_std"
conn_predictors <- c(
  eigen = "eigen_centrality_std",
  in_strength = "in_strength_std",
  biomass_supply_log = "biomass_supply_log_std"
)

# base structures: rhs is the non-connectivity part; conn = TRUE models get a
# connectivity predictor appended (once per conn_predictors)
model_specs <- tibble::tribble(
  ~model, ~spatial, ~rhs, ~conn,
  "space", "on", survey_term, FALSE,
  "space_env", "on", paste(survey_term, env_terms, sep = " + "), FALSE,
  "env", "off", paste(survey_term, env_terms, sep = " + "), FALSE,
  "space_conn", "on", survey_term, TRUE,
  "space_env_conn", "on", paste(survey_term, env_terms, sep = " + "), TRUE,
  "env_conn", "off", paste(survey_term, env_terms, sep = " + "), TRUE
)

# expand connectivity models over the rotating predictors -> concrete predictor
# sets, each with its own fixed spatial setting
models <- model_specs |>
  mutate(conn_name = if_else(conn, list(names(conn_predictors)), list(NA_character_))) |>
  unnest(conn_name) |>
  mutate(
    id = if_else(conn, paste(model, conn_name, sep = "_"), model),
    rhs = if_else(conn, paste(rhs, conn_predictors[conn_name], sep = " + "), rhs)
  )
predictor_sets <- set_names(models$rhs, models$id)

# 01 Load data ----
dat <- readRDS(here("data", "final", "cockles_connectivity.rds")) |>
  mutate(survey = factor(survey))

# 02 Standardise covariates ----
std_vars <- c(
  "depth", "temp", "sal", "oxy", "chla", "shear_max", "shear_mean", "phyto",
  "biomass_supply", "biomass_supply_log", "present_supply",
  "biomass_export", "biomass_export_log", "present_export",
  "in_strength", "local_retention", "eigen_centrality"
)
dat <- dat |>
  mutate(across(all_of(std_vars), ~ as.numeric(scale(.x)), .names = "{.col}_std"))

# 03 Shared model dataset ----
# keep rows complete on every covariate used by any model, so all fits share n
used_vars <- predictor_sets |>
  map(~ all.vars(as.formula(paste("~", .x)))) |>
  unlist() |>
  unique()
dat <- dat |>
  filter(
    !is.na(biomass), !is.na(present), !is.na(x_utm), !is.na(y_utm),
    if_all(all_of(used_vars), ~ !is.na(.x))
  )

# 04 Collinearity ----
# correlation across all predictors (overview), then VIF WITHIN each model -
# pooling all covariates would inflate VIF for redundant alternatives (e.g.
# biomass_supply vs biomass_supply_log vs in_strength) that never share a model
coll_vars <- intersect(used_vars, paste0(std_vars, "_std"))
cat("Predictor correlation matrix:\n")
print(round(cor(dat[, coll_vars]), 2))

# VIF within each distinct predictor block (space on / off does not change it,
# so dedupe by the covariate set)
cat("\nVIF within each predictor set (>~5 = concerning):\n")
seen <- character()
for (pname in names(predictor_sets)) {
  vars <- intersect(all.vars(as.formula(paste("~", predictor_sets[[pname]]))), coll_vars)
  key <- paste(sort(vars), collapse = "+")
  if (length(vars) >= 2 && !key %in% seen) {
    seen <- c(seen, key)
    cat("--", pname, "--\n")
    print(round(vif(dat, vars), 2))
  }
}

# 05 Coastline barrier mesh ----
land <- ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf") |>
  st_make_valid()
region <- st_bbox(st_transform(
  st_as_sf(dat |> transmute(x = x_utm * 1000, y = y_utm * 1000), coords = c("x", "y"), crs = 32632),
  4326
))
region["xmin"] <- region["xmin"] - 0.4
region["ymin"] <- region["ymin"] - 0.4
region["xmax"] <- region["xmax"] + 0.4
region["ymax"] <- region["ymax"] + 0.4
land_region <- suppressWarnings(st_crop(land, region)) |> st_transform(32632)

mesh <- make_mesh(dat, c("x_utm", "y_utm"), cutoff = 1) # 1 km min edge
barrier_mesh <- sdmTMBextra::add_barrier_mesh(
  mesh, land_region,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)

# 06 Fit all models ----
# spatial field is on / off per model (see model_specs); "space" vs matching
# no-field model shows how much residual spatial autocorrelation the field soaks up
responses <- list(
  present = binomial(link = "logit"),
  biomass = tweedie(link = "log")
)

fits <- list()
comparison <- tibble()
for (rname in names(responses)) {
  for (i in seq_len(nrow(models))) {
    pname <- models$id[i]
    sp <- models$spatial[i]
    model_id <- paste(rname, pname, sep = "_")
    form <- as.formula(paste(rname, "~", predictor_sets[[pname]]))
    fit <- tryCatch(
      sdmTMB(form, data = dat, mesh = barrier_mesh, spatial = sp, family = responses[[rname]]),
      error = function(e) {
        message("failed: ", model_id, " - ", conditionMessage(e))
        NULL
      }
    )
    converged <- !is.null(fit) && isTRUE(suppressMessages(sanity(fit))$all_ok)
    fits[[model_id]] <- fit
    comparison <- bind_rows(comparison, tibble(
      response = rname,
      model = pname,
      structure = models$model[i],
      spatial = sp,
      formula = deparse1(form),
      converged = converged,
      aic = if (converged) AIC(fit) else NA_real_
    ))
  }
}

# rank within each response by AIC
comparison <- comparison |>
  group_by(response) |>
  mutate(delta_aic = aic - min(aic, na.rm = TRUE)) |>
  arrange(response, aic) |>
  ungroup()

print(comparison, n = Inf)

# 07 Coefficients of every model ----
for (rname in names(responses)) {
  for (pname in models$id) {
    model_id <- paste(rname, pname, sep = "_")
    cat("\n===", rname, "~", pname, "===\n")
    f <- fits[[model_id]]
    if (is.null(f)) {
      cat("(failed to converge)\n")
    } else {
      print(tidy(f, effects = "fixed", conf.int = TRUE))
    }
  }
}

# 08 Save ----
saveRDS(fits, here("data", "intermediate", "sdm_fits.rds"))
saveRDS(comparison, here("data", "final", "sdm_model_comparison.rds"))
