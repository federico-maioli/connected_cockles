# Shared helpers for the cockle SDM pipeline. Callers must already have
# library(tidyverse), library(here), library(sf), and library(sdmTMB) loaded.

# compatibility shim: sdmTMBextra 0.0.5 calls fmesher::fm_identical_CRS, renamed
# to fm_crs_is_identical in fmesher >= 0.7; register the old name as an alias.
# required for add_barrier_mesh() below - do not remove unless sdmTMBextra is updated
local({
  exps <- .getNamespaceInfo(asNamespace("fmesher"), "exports")
  if (!exists("fm_identical_CRS", envir = exps, inherits = FALSE)) {
    assign("fm_identical_CRS", "fm_crs_is_identical", envir = exps)
  }
})

# coastline barrier mesh for the cockle survey extent: local coastline (UTM
# 32N, metres) cropped to the data extent with a margin, an fmesher mesh with a
# capped maximum edge, and a land barrier so correlation does not cross
# headlands. dat needs x_utm/y_utm in km, as in cockles_clean.rds.
build_barrier_mesh <- function(dat, margin = 30000) {
  land_utm <- st_read(
    here("data", "raw", "boundaries", "land_small_utm", "land_small_utm.shp"),
    quiet = TRUE
  ) |>
    st_make_valid()

  region <- st_bbox(
    c(
      xmin = min(dat$x_utm) * 1000 - margin,
      ymin = min(dat$y_utm) * 1000 - margin,
      xmax = max(dat$x_utm) * 1000 + margin,
      ymax = max(dat$y_utm) * 1000 + margin
    ),
    crs = st_crs(32632)
  )
  land_region <- suppressWarnings(st_crop(land_utm, region))

  # spatial range is ~5 km, so a 2 km inner edge resolves the field; 20 km
  # outer offset keeps the boundary away from the data (dat coords are in km)
  inla_mesh <- fmesher::fm_mesh_2d_inla(
    loc = cbind(dat$x_utm, dat$y_utm),
    max.edge = c(2, 10),
    offset = c(5, 20),
    cutoff = 1
  )
  mesh <- make_mesh(dat, c("x_utm", "y_utm"), mesh = inla_mesh)
  sdmTMBextra::add_barrier_mesh(
    mesh, land_region,
    range_fraction = 0.1,
    proj_scaling = 1000,
    plot = FALSE
  )
}

# Nakagawa & Schielzeth marginal/conditional R2 and a variance share per
# component, for one fitted sdmTMB model (Bernoulli, nbinom2, or Tweedie).
# See 03_variance_partitioning.R for the full derivation and a check against
# insight::get_variance(). blocks: named list of model-matrix columns to
# report separately, e.g. list(env = c("depth_std", "temp_std")); columns not
# named in any block are pooled into "other".
#
# The fixed-effects share is split among blocks HMSC-style (Tikhonov et al.
# 2020 MEE, computeVariancePartitioning): each block's own variance is taken
# from its own columns only, ignoring covariance with every other block, which
# keeps every block's share >= 0 even under collinearity (unlike splitting by
# each block's covariance with the total, which is exact but can go negative
# - see 03_variance_partitioning.R). Unlike HMSC, which normalises the blocks
# to sum to its own total fixed-effect variance, they are rescaled here to sum
# exactly to var_fixed - the same total this function reports elsewhere and
# the one checked against insight::get_variance() in test_nakagawa_sdmtmb.R -
# so the reported marginal R2 is unaffected by how it's split among blocks.
nakagawa_sdmtmb <- function(fit, blocks = NULL) {
  if (isTRUE(fit$family$delta)) stop("delta models are not supported")
  fam <- fit$family$family
  link <- fit$family$link
  rp <- tidy(fit, effects = "ran_pars")
  ran_par <- function(term) rp$estimate[rp$term == term]

  # fixed-effects linear predictor, term by term
  X <- fit$tmb_data$X_ij[[1]]
  coefs <- tidy(fit, effects = "fixed") |>
    select(term, estimate) |>
    deframe()
  lp <- sweep(X, 2, coefs[colnames(X)], "*")
  eta <- rowSums(lp)

  slopes <- setdiff(colnames(X), "(Intercept)")
  if (is.null(blocks)) blocks <- list(fixed = slopes)
  other <- setdiff(slopes, unlist(blocks))
  if (length(other) > 0) blocks$other <- other
  var_fixed <- var(eta)
  var_blocks_own <- map_dbl(blocks, \(cols) var(rowSums(lp[, cols, drop = FALSE])))
  # a model with no predictors (e.g. response ~ 1) has var_blocks_own all 0;
  # rescaling would divide 0 by 0, so leave every block at 0 instead of NaN
  var_blocks <- if (sum(var_blocks_own) > 0) var_blocks_own / sum(var_blocks_own) * var_fixed else var_blocks_own

  # random effects; sum() of a missing term is 0, so absent fields drop out
  var_spatial <- sum(ran_par("sigma_O")^2)
  var_spatiotemporal <- sum(ran_par("sigma_E")^2)
  var_intercept <- sum(rp$estimate[startsWith(rp$term, "sd__")]^2)
  var_random <- var_spatial + var_spatiotemporal + var_intercept

  mu <- mean(exp(eta)) * exp(0.5 * var_random)
  var_dist <- switch(fam,
    binomial = switch(link,
      logit = pi^2 / 3,
      probit = 1,
      cloglog = pi^2 / 6,
      stop("unsupported link: ", link)
    ),
    nbinom2 = log1p(1 / mu + 1 / ran_par("phi")),
    tweedie = log1p(ran_par("phi") * mu^(ran_par("tweedie_p") - 2)),
    stop("unsupported family: ", fam)
  )

  var_all <- c(
    var_blocks,
    spatial = var_spatial,
    spatiotemporal = var_spatiotemporal,
    random_intercept = var_intercept,
    distribution = var_dist
  )
  total <- var_fixed + var_random + var_dist
  tibble(
    component = names(var_all),
    variance = unname(var_all),
    share = variance / total,
    r2_marginal = var_fixed / total,
    r2_conditional = (var_fixed + var_random) / total
  )
}
