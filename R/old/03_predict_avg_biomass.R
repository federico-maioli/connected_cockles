# Model average cockle biomass and probability of presence from depth, then
# predict both onto the 2 km connectivity grid.
#
# Two spatial sdmTMB models on an fmesher mesh with a land barrier (spatial
# correlation does not cross land) and survey gear as a factor:
#   biomass ~ depth + survey   Tweedie (log link)  -> average biomass per cell
#   present ~ depth + survey   binomial (logit)    -> probability of presence per cell
# Predictions are made for the "Stock" survey and attached back onto the grid in
# its original id order.

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
library(patchwork)

# compatibility shim: sdmTMBextra 0.0.5 calls fmesher::fm_identical_CRS, renamed
# to fm_crs_is_identical in fmesher >= 0.7; register the old name as an alias.
# required for add_barrier_mesh() below - do not remove unless sdmTMBextra is updated
local({
  exps <- .getNamespaceInfo(asNamespace("fmesher"), "exports")
  if (!exists("fm_identical_CRS", envir = exps, inherits = FALSE)) {
    assign("fm_identical_CRS", "fm_crs_is_identical", envir = exps)
  }
})

# 01 Load data ----
# cleaned survey points already carry raster depth and coordinates in km
dat <- readRDS(here("data", "derived", "cockles_clean.rds")) |>
  filter(!is.na(biomass), !is.na(present), !is.na(depth), !is.na(x_utm), !is.na(y_utm)) |>
  mutate(survey = factor(survey))

grid <- readRDS(here("data", "derived", "spatial_grid.rds"))

# local coastline for the land mask and plots, already UTM 32N in metres
land_utm <- st_read(
  here("data", "boundaries", "land_small_utm", "land_small_utm.shp"),
  quiet = TRUE
) |>
  st_make_valid()

# 02 Coastline and spatial mesh ----
# crop the coastline to the survey and grid extent with a 30 km margin, used for
# the land mask and the plots (dat is in km, grid in metres)
margin <- 30000
region <- st_bbox(
  c(
    xmin = min(min(dat$x_utm) * 1000, min(grid$x_utm)) - margin,
    ymin = min(min(dat$y_utm) * 1000, min(grid$y_utm)) - margin,
    xmax = max(max(dat$x_utm) * 1000, max(grid$x_utm)) + margin,
    ymax = max(max(dat$y_utm) * 1000, max(grid$y_utm)) + margin
  ),
  crs = st_crs(32632)
)
land_region <- suppressWarnings(st_crop(land_utm, region))

# fmesher mesh with an explicit maximum edge length: the spatial range is ~5 km,
# so a 2 km inner edge resolves the field, a 10 km outer edge keeps the extension
# cheap, and the 20 km outer offset holds the boundary away from the data
inla_mesh <- fmesher::fm_mesh_2d_inla(
  loc = cbind(dat$x_utm, dat$y_utm), # coordinates in km
  max.edge = c(2, 10), # max triangle edge; inner and outer meshes
  offset = c(5, 20), # inner and outer border widths
  cutoff = 1 # minimum triangle edge length
)
mesh <- make_mesh(dat, c("x_utm", "y_utm"), mesh = inla_mesh)

# barrier mesh: correlation is downweighted across land, so it does not cross
# the fjord's headlands (proj_scaling 1000 because the mesh is built in km while
# the land polygon is in metres). fmesher supplies the mesh geometry; the barrier
# FEM itself comes from INLAspacetime via sdmTMBextra
barrier_mesh <- sdmTMBextra::add_barrier_mesh(
  mesh, land_region,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)
cat(
  "mesh vertices:", mesh$mesh$n,
  "| water triangles:", length(barrier_mesh$normal_triangles),
  "| land triangles:", length(barrier_mesh$barrier_triangles), "\n"
)

# 03 Biomass model ----
fit_biomass <- sdmTMB(
  biomass ~ depth + survey,
  data = dat,
  mesh = barrier_mesh,
  spatial = "on",
  family = tweedie(link = "log")
)
print(sanity(fit_biomass))

# check Tweedie fit: randomized-quantile residuals should be ~ N(0, 1)
res_bio <- residuals(fit_biomass)
qqnorm(res_bio)
qqline(res_bio)

# 04 Presence model ----
fit_present <- sdmTMB(
  present ~ depth + survey,
  data = dat,
  mesh = barrier_mesh,
  spatial = "on",
  family = binomial(link = "logit")
)
print(sanity(fit_present))

# 05 Predict onto the grid (survey = Stock) ----
# grid coordinates are in metres; the model is in km, so convert for prediction.
# only wet cells (non-missing depth) can be predicted
wet <- !is.na(grid$depth)
pred_grid <- grid |>
  filter(wet) |>
  mutate(
    x_utm = x_utm / 1000,
    y_utm = y_utm / 1000,
    survey = factor("Stock", levels = levels(dat$survey))
  )

# the Tweedie link is log, so predict() returns biomass on the log scale and
# exp() puts it back on the response scale; both are kept
log_avg_biomass <- predict(fit_biomass, newdata = pred_grid)$est # log scale
avg_biomass <- exp(log_avg_biomass) # response scale
prob_present <- plogis(predict(fit_present, newdata = pred_grid)$est) # probability

# uncertainty on the biomass prediction only: draws from the joint precision
# matrix, so the fixed effects and the spatial field are sampled together and
# their covariance is respected. nsim returns log-scale draws (cells x nsim).
# The two SDs are summarised on their own scale and are NOT interconvertible:
# exp() is non-linear, and with sigma_O ~ 4.6 the exponentiated draws are
# lognormal with a heavy upper tail, so avg_biomass_sd is large and skewed while
# log_avg_biomass_sd is a roughly multiplicative error
n_sim <- 500
set.seed(1)
sim_link <- predict(fit_biomass, newdata = pred_grid, nsim = n_sim)
avg_biomass_sd <- apply(exp(sim_link), 1, sd)
log_avg_biomass_sd <- apply(sim_link, 1, sd)
rm(sim_link)

# 06 Mask cells far from any survey point ----
# spatial field extrapolates; cells > max_dist_km from data are set to 0
max_dist_km <- 4
obs <- cbind(dat$x_utm, dat$y_utm)
nn_dist <- apply(cbind(pred_grid$x_utm, pred_grid$y_utm), 1, function(p) {
  sqrt(min((obs[, 1] - p[1])^2 + (obs[, 2] - p[2])^2))
})
far <- nn_dist > max_dist_km
avg_biomass[far] <- 0
prob_present[far] <- 0
# the response-scale estimates are forced to 0 so the grid stays aligned with the
# flow matrix, but a forced 0 is an assumption rather than a prediction, so it
# carries no uncertainty. The log-scale columns are NA rather than 0 there too:
# log(0) is undefined, and a literal 0 would wrongly read as a biomass of 1
avg_biomass_sd[far] <- NA_real_
log_avg_biomass[far] <- NA_real_
log_avg_biomass_sd[far] <- NA_real_

# 07 Attach to grid ----
grid$avg_biomass <- 0
grid$prob_present <- 0
grid$avg_biomass[wet] <- avg_biomass
grid$prob_present[wet] <- prob_present
# NA everywhere by default, so dry cells (never predicted) are NA too
for (nm in c("avg_biomass_sd", "log_avg_biomass", "log_avg_biomass_sd")) {
  grid[[nm]] <- NA_real_
  grid[[nm]][wet] <- get(nm)
}

# 08 Zero out land cells ----
# use the local land_small_utm coastline to set any on-land cell to 0
on_land <- lengths(st_intersects(
  st_as_sf(grid, coords = c("x_utm", "y_utm"), crs = 32632),
  land_utm
)) > 0
grid$avg_biomass[on_land] <- 0
grid$prob_present[on_land] <- 0
grid[on_land, c("avg_biomass_sd", "log_avg_biomass", "log_avg_biomass_sd")] <- NA_real_

kept <- grid$avg_biomass > 0
cat(sprintf(
  paste0(
    "cells: %d | predicted: %d | masked (log columns and SDs = NA): %d\n",
    "  median SD, response scale: %.0f | median SD, log scale: %.2f\n"
  ),
  nrow(grid), sum(kept), sum(!kept),
  median(grid$avg_biomass_sd[kept]), median(grid$log_avg_biomass_sd[kept])
))

# 09 Save ----
# land / masked cells are 0 (not dropped), so the grid keeps all 2340 cells and
# stays aligned cell-for-cell with the flow matrix
flow <- readRDS(here("data", "derived", "flow_matrix_all.rds"))
# the response-scale estimates must be complete; the log columns and the SDs are
# NA on masked cells by design, and present on every predicted cell
stopifnot(
  nrow(grid) == nrow(flow),
  !anyNA(grid$avg_biomass), !anyNA(grid$prob_present),
  !anyNA(grid[kept, c("avg_biomass_sd", "log_avg_biomass", "log_avg_biomass_sd")])
)

saveRDS(grid, here("data", "derived", "avg_biomass_grid.rds"))

# 10 Visual check of the predictions ----
# diagnostic only, not saved: a quick look at whether the predicted surfaces are
# sensible before they are fed to 04. Publication figures are built by the
# fig*/figS* scripts from the saved outputs, never here.
# tiles for the wet cells, coastline drawn on top; two panels share the extent
plot_dat <- grid |> filter(!is.na(depth))
xlim <- range(grid$x_utm)
ylim <- range(grid$y_utm)

p_biomass <- ggplot() +
  geom_tile(data = plot_dat, aes(x_utm, y_utm, fill = avg_biomass)) +
  geom_sf(data = land_region, fill = "grey85", colour = "grey60", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  scale_fill_viridis_c(trans='sqrt') +
  labs(title = "Average biomass", x = NULL, y = NULL, fill = "Biomass") +
  theme_light()

p_present <- ggplot() +
  geom_tile(data = plot_dat, aes(x_utm, y_utm, fill = prob_present)) +
  geom_sf(data = land_region, fill = "grey85", colour = "grey60", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  scale_fill_viridis_c(option = "mako") +
  labs(title = "Probability of presence", x = NULL, y = NULL, fill = "P(presence)") +
  theme_light()

# the log-scale panels: masked cells are NA and drop out of the tile layer, so
# these show exactly the cells that carry a prediction
p_log_biomass <- ggplot() +
  geom_tile(
    data = filter(plot_dat, !is.na(log_avg_biomass)),
    aes(x_utm, y_utm, fill = log_avg_biomass)
  ) +
  geom_sf(data = land_region, fill = "grey85", colour = "grey60", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  scale_fill_viridis_c() +
  labs(title = "Average biomass (log)", x = NULL, y = NULL, fill = "log biomass") +
  theme_light()

p_log_biomass_sd <- ggplot() +
  geom_tile(
    data = filter(plot_dat, !is.na(log_avg_biomass_sd)),
    aes(x_utm, y_utm, fill = log_avg_biomass_sd)
  ) +
  geom_sf(data = land_region, fill = "grey85", colour = "grey60", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  scale_fill_viridis_c(option = "rocket") +
  labs(title = "Biomass SD (log)", x = NULL, y = NULL, fill = "SD") +
  theme_light()

(p_biomass + p_present) / (p_log_biomass + p_log_biomass_sd)






