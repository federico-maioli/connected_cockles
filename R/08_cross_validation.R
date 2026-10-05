# Spatially blocked 10-fold cross-validation of the four main models from
# 07_fit_sdm.R (space, space_env, space_conn, space_env_conn): does adding the
# environment or in-strength improve out-of-sample prediction?
#
# Blocks are the 2 x 2 km connectivity grid cells: all stations in a cell are
# held out together, and cells are randomly assigned to 10 folds, the same
# folds for every model. Held-out predictions use predictive = "mle-mvn",
# which integrates over parameter uncertainty by sampling from the
# multivariate normal of the MLE (sdmTMB >= 1.1.0.9024, development version).
# Models are compared on the expected log predictive density (ELPD) with
# loo::loo_compare(), and pairwise for each added block (ELPD difference and
# its standard error across observations), as in the sdmTMB
# cross-validation vignette.
#
# Saved in data/sdm/cv/: cv_results.rds (per-model sdmTMB_cv output, models
# not kept), elpd_compare.rds (loo_compare table) and elpd_pairwise.rds.

library(tidyverse)
library(here)
library(sdmTMB)

source(here("R", "helpers.R")) # fmesher compatibility shim

k_folds <- 10
n_col <- 65
x_min <- 450074
y_min <- 6258093
cell_size <- 2000

# 01 Data, mesh and folds ----
# the model dataset and barrier mesh of the fitted full model, so the
# cross-validation uses exactly the data and mesh of 07_fit_sdm.R
fit_full <- readRDS(here("data", "sdm", "main", "space_env_conn.rds"))
dat <- fit_full$data
mesh <- fit_full$spde

cell <- (floor((dat$x_utm * 1000 - x_min) / cell_size) + 1) + floor((dat$y_utm * 1000 - y_min) / cell_size) * n_col
set.seed(1)
cells <- unique(cell)
cell_fold <- setNames(sample(rep_len(1:k_folds, length(cells))), cells)
dat$fold <- unname(cell_fold[as.character(cell)])
count(dat, fold) |>
  print()

# 02 Cross-validate ----
env <- "depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std"
models <- tribble(
  ~model, ~rhs,
  "space", "survey + year",
  "space_env", paste("survey + year +", env),
  "space_conn", "survey + year + conn_in_strength_std",
  "space_env_conn", paste("survey + year +", env, "+ conn_in_strength_std")
)

cv_results <- models |>
  mutate(cv = map(rhs, \(rhs) {
    sdmTMB_cv(
      as.formula(paste("biomass ~", rhs)),
      data = dat, mesh = mesh, spatial = "on", family = delta_gamma(),
      fold_ids = dat$fold, k_folds = k_folds, save_models = FALSE, predictive = "mle-mvn"
    )
  }))

# 03 Compare ----
elpd <- set_names(map(cv_results$cv, loo::elpd), cv_results$model)

elpd_compare <- loo::loo_compare(elpd) |>
  as_tibble() |>
  select(model, elpd_diff, se_diff, elpd, se_elpd) |>
  left_join(
    cv_results |> transmute(model, converged = map_lgl(cv, \(x) all(x$converged))),
    by = "model"
  )
print(elpd_compare)

# each added block: ELPD of the model with it minus the model without it
pairs <- tribble(
  ~comparison, ~with, ~without,
  "Environment added to space", "space_env", "space",
  "In-strength added to space", "space_conn", "space",
  "In-strength added to space + environment", "space_env_conn", "space_env",
  "Environment added to space + in-strength", "space_env_conn", "space_conn"
)
elpd_pairwise <- pairs |>
  mutate(cmp = map2(with, without, \(a, b) {
    diff <- elpd[[a]]$pointwise[, "elpd"] - elpd[[b]]$pointwise[, "elpd"]
    tibble(elpd_diff = sum(diff), se_diff = sqrt(length(diff)) * sd(diff))
  })) |>
  unnest(cmp)
print(elpd_pairwise)

# 04 Save ----
dir.create(here("data", "sdm", "cv"), showWarnings = FALSE, recursive = TRUE)
saveRDS(cv_results, here("data", "sdm", "cv", "cv_results.rds"))
saveRDS(elpd_compare, here("data", "sdm", "cv", "elpd_compare.rds"))
saveRDS(elpd_pairwise, here("data", "sdm", "cv", "elpd_pairwise.rds"))
