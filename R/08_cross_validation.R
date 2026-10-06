# Spatially blocked 10-fold cross-validation of the four main models from
# 07_fit_sdm.R (space, space_env, space_conn, space_env_conn): does adding the
# environment or in-strength improve out-of-sample prediction of presence, of
# biomass where present, or of both?
#
# The two parts of the delta-gamma model share no parameters, so each is
# cross-validated separately with the same predictors as in 07_fit_sdm.R:
# presence (binomial, biomass > 0, all stations) and biomass where present
# (gamma, stations with cockles). The total is their sum for each station
# (the log predictive density of a delta model is the presence term plus, for
# stations with cockles, the biomass term).
#
# Blocks are the 2 x 2 km connectivity grid cells: all stations in a cell are
# held out together, and cells are randomly assigned to 10 folds, the same
# folds for every model. Presence uses all cells; biomass where present uses
# its own folds, drawn over the cells holding stations with cockles. Held-out predictions use
# predictive = "mle-mvn", which integrates over parameter uncertainty by
# sampling from the multivariate normal of the MLE (sdmTMB >= 1.1.0.9024,
# development version). Models are compared on the expected log predictive
# density (ELPD) with loo::loo_compare(), and pairwise for each added block
# (ELPD difference and its standard error across observations), as in the
# sdmTMB cross-validation vignette.
#
# Saved in data/sdm/cv/: cv_results.rds (per-model sdmTMB_cv output of each
# part, models not kept), elpd_compare.rds and elpd_pairwise.rds (one block
# per part: presence, biomass, total).

library(tidyverse)
library(here)
library(sdmTMB)

source(here("R", "helpers.R")) # fmesher compatibility shim, needed by add_barrier_mesh() below

k_folds <- 10
n_draws <- 500 # mle-mvn draws per held-out prediction (default 100 gave noisy ELPD differences)
n_col <- 65
x_min <- 450074
y_min <- 6258093
cell_size <- 2000

# 01 Data, meshes and folds ----
# the model dataset and barrier mesh of the fitted full model, so the
# cross-validation uses exactly the data and mesh of 07_fit_sdm.R
fit_full <- readRDS(here("data", "sdm", "main", "space_env_conn.rds"))
dat <- fit_full$data |>
  mutate(present_biomass = as.integer(biomass > 0))
mesh <- fit_full$spde

cell <- (floor((dat$x_utm * 1000 - x_min) / cell_size) + 1) + floor((dat$y_utm * 1000 - y_min) / cell_size) * n_col
set.seed(1)
cells <- unique(cell)
cell_fold <- setNames(sample(rep_len(1:k_folds, length(cells))), cells)
dat$fold <- unname(cell_fold[as.character(cell)])
count(dat, fold) |>
  print()

# biomass where present: the stations with cockles, on the same mesh and
# barrier, with their own folds - the 2 x 2 km cells holding at least one
# station with cockles, randomly assigned to the 10 folds, so every fold holds
# out a similar number of occupied cells (the presence folds, drawn over all
# cells, would leave some biomass folds with very few occupied cells)
dat$row <- seq_len(nrow(dat))
dat$cell <- cell
dat_pos <- filter(dat, biomass > 0)
set.seed(1)
cells_pos <- unique(dat_pos$cell)
cell_fold_pos <- setNames(sample(rep_len(1:k_folds, length(cells_pos))), cells_pos)
dat_pos$fold <- unname(cell_fold_pos[as.character(dat_pos$cell)])
dat_pos |>
  summarise(stations = n(), cells = n_distinct(cell), .by = fold) |>
  arrange(fold) |>
  print()
mesh_pos <- sdmTMBextra::add_barrier_mesh(
  make_mesh(dat_pos, c("x_utm", "y_utm"), mesh = mesh$mesh),
  readRDS(here("data", "mesh", "mesh.rds"))$land_barrier,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)

# 02 Cross-validate ----
env <- "depth_std + I(depth_std^2) + temp_std + sal_std + oxy_std + shear_max_std"
models <- tribble(
  ~model, ~rhs,
  "space", "survey + year",
  "space_env", paste("survey + year +", env),
  "space_conn", "survey + year + conn_in_strength_std",
  "space_env_conn", paste("survey + year +", env, "+ conn_in_strength_std")
)

# the mle-mvn draws are random, so each cross-validation is seeded
set.seed(1)
cv_results <- models |>
  mutate(
    cv_presence = map(rhs, \(rhs) {
      sdmTMB_cv(
        as.formula(paste("present_biomass ~", rhs)),
        data = dat, mesh = mesh, spatial = "on", family = binomial(),
        fold_ids = dat$fold, k_folds = k_folds, save_models = FALSE, predictive = "mle-mvn", nsim = n_draws
      )
    }),
    cv_biomass = map(rhs, \(rhs) {
      sdmTMB_cv(
        as.formula(paste("biomass ~", rhs)),
        data = dat_pos, mesh = mesh_pos, spatial = "on", family = Gamma(link = "log"),
        fold_ids = dat_pos$fold, k_folds = k_folds, save_models = FALSE, predictive = "mle-mvn", nsim = n_draws
      )
    })
  )

# 03 Compare ----
# pointwise ELPD of each part; the total adds the biomass term to the
# presence term of the stations with cockles
elpd_presence <- set_names(map(cv_results$cv_presence, loo::elpd), cv_results$model)
elpd_biomass <- set_names(map(cv_results$cv_biomass, loo::elpd), cv_results$model)
elpd_total <- map2(elpd_presence, elpd_biomass, \(p, b) {
  pointwise <- p$pointwise[, "elpd"]
  pointwise[dat_pos$row] <- pointwise[dat_pos$row] + b$pointwise[, "elpd"]
  loo::elpd(matrix(pointwise, nrow = 1))
})
elpd_parts <- list(presence = elpd_presence, biomass = elpd_biomass, total = elpd_total)

converged <- cv_results |>
  transmute(model, converged = map2_lgl(cv_presence, cv_biomass, \(p, b) all(p$converged) && all(b$converged)))

elpd_compare <- imap(elpd_parts, \(elpd, part) {
  loo::loo_compare(elpd) |>
    as_tibble() |>
    select(model, elpd_diff, se_diff, elpd, se_elpd) |>
    mutate(part = part)
}) |>
  list_rbind() |>
  left_join(converged, by = "model")
print(elpd_compare, n = Inf)

# each added block: ELPD of the model with it minus the model without it
pairs <- tribble(
  ~comparison, ~with, ~without,
  "Environment added to space", "space_env", "space",
  "In-strength added to space", "space_conn", "space",
  "In-strength added to space + environment", "space_env_conn", "space_env",
  "Environment added to space + in-strength", "space_env_conn", "space_conn"
)
elpd_pairwise <- imap(elpd_parts, \(elpd, part) {
  pairs |>
    mutate(cmp = map2(with, without, \(a, b) {
      diff <- elpd[[a]]$pointwise[, "elpd"] - elpd[[b]]$pointwise[, "elpd"]
      tibble(elpd_diff = sum(diff), se_diff = sqrt(length(diff)) * sd(diff))
    })) |>
    unnest(cmp) |>
    mutate(part = part)
}) |>
  list_rbind()
print(elpd_pairwise, n = Inf)

# 04 Save ----
dir.create(here("data", "sdm", "cv"), showWarnings = FALSE, recursive = TRUE)
saveRDS(cv_results, here("data", "sdm", "cv", "cv_results.rds"))
saveRDS(elpd_compare, here("data", "sdm", "cv", "elpd_compare.rds"))
saveRDS(elpd_pairwise, here("data", "sdm", "cv", "elpd_pairwise.rds"))
