# Model average cockle biomass and probability of presence from depth, then
# predict both onto the 2 km connectivity grid.
#
# Two spatial sdmTMB models with a coastline barrier mesh (spatial correlation
# does not cross land) and survey gear as a factor:
#   biomass ~ depth + survey   Tweedie (log link)  -> average biomass per cell
#   present ~ depth + survey   binomial (logit)    -> probability of presence per cell
# Predictions are made for the "Stock" survey and attached back onto the grid in
# its original id order.

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
library(rnaturalearth)
library(patchwork)

# 01 Load data ----
# cleaned survey points already carry raster depth and coordinates in km
dat <- readRDS(here("data", "final", "cockles_clean.rds")) |>
  filter(!is.na(biomass), !is.na(present), !is.na(depth), !is.na(x_utm), !is.na(y_utm)) |>
  mutate(survey = factor(survey))

grid <- readRDS(here("data", "final", "spatial_grid.rds"))

# 02 Coastline barrier mesh ----
# Natural Earth land, cropped to the grid area and projected to UTM metres
land <- ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf") |>
  st_make_valid()

region <- st_bbox(st_transform(st_as_sf(grid, coords = c("x_utm", "y_utm"), crs = 32632), 4326))
region["xmin"] <- region["xmin"] - 0.4
region["ymin"] <- region["ymin"] - 0.4
region["xmax"] <- region["xmax"] + 0.4
region["ymax"] <- region["ymax"] + 0.4
land_region <- suppressWarnings(st_crop(land, region)) |> st_transform(32632)

# barrier mesh: correlation is downweighted across land (proj_scaling 1000
# because the mesh is built in km while the land polygon is in metres)
mesh <- make_mesh(dat, c("x_utm", "y_utm"), cutoff = 1) # 1 km min edge
barrier_mesh <- sdmTMBextra::add_barrier_mesh(
  mesh, land_region,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
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

avg_biomass <- exp(predict(fit_biomass, newdata = pred_grid)$est) # response scale
prob_present <- plogis(predict(fit_present, newdata = pred_grid)$est) # probability

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

# 07 Attach to grid ----
grid$avg_biomass <- 0
grid$prob_present <- 0
grid$avg_biomass[wet] <- avg_biomass
grid$prob_present[wet] <- prob_present

# 08 Zero out land cells ----
# use the same Natural Earth land to set any on-land cell to 0
on_land <- lengths(st_intersects(
  st_as_sf(grid, coords = c("x_utm", "y_utm"), crs = 32632),
  land_region
)) > 0
grid$avg_biomass[on_land] <- 0
grid$prob_present[on_land] <- 0

# 09 Save ----
# land / masked cells are 0 (not dropped), so the grid keeps all 2340 cells and
# stays aligned cell-for-cell with the flow matrix
flow <- readRDS(here("data", "intermediate", "flow_matrix.rds"))
stopifnot(nrow(grid) == nrow(flow), !anyNA(grid$avg_biomass), !anyNA(grid$prob_present))

saveRDS(grid, here("data", "final", "avg_biomass_grid.rds"))

# 10 Plot predictions ----
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

p_biomass + p_present

