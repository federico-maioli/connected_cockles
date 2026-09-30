# Fit a "suitability" surface for biomass and presence: intercept + spatial
# field only, no environmental covariates. Predict it onto the connectivity
# grid; the predicted values become the weights used to collapse the raw MIKE
# connectivity matrix into weighted connectivity metrics in
# 04_weight_connectivity.R.
#
# No covariates: checked directly, environmental terms were unstable once the
# spatial field is on (temp/sal/oxy lose significance - classic spatial
# confounding) and most of what an environment-only fit was capturing was
# actually spatial structure. A pure spatial field is a more honest
# description of what the data supports.

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)

source(here("R", "helpers.R")) # fmesher compatibility shim, needed by add_barrier_mesh() below

# 01 Load data ----
dat <- readRDS(here("data", "cockles", "derived", "cockles_env.rds"))
grid <- readRDS(here("data", "grid", "grid_env.rds"))

dat <- dat |>
  filter(!is.na(biomass), !is.na(present), !is.na(x_utm), !is.na(y_utm))

# 02 Coastline barrier mesh ----
land <- st_read(
  here("data", "boundaries", "land_small_utm", "land_small_utm.shp"),
  quiet = TRUE
) |>
  st_make_valid()

region <- st_bbox(
  c(
    xmin = min(dat$x_utm) * 1000 - 30000,
    ymin = min(dat$y_utm) * 1000 - 30000,
    xmax = max(dat$x_utm) * 1000 + 30000,
    ymax = max(dat$y_utm) * 1000 + 30000
  ),
  crs = st_crs(32632)
)
land_region <- suppressWarnings(st_crop(land, region))

inla_mesh <- fmesher::fm_mesh_2d_inla(
  loc = cbind(dat$x_utm, dat$y_utm),
  max.edge = c(2, 10),
  offset = c(5, 20),
  cutoff = 1
)
barrier_mesh <- sdmTMBextra::add_barrier_mesh(
  make_mesh(dat, c("x_utm", "y_utm"), mesh = inla_mesh),
  land_region,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)

# 03 Fit suitability models ----
fit_biomass <- sdmTMB(
  biomass ~ 1,
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)

fit_presence <- sdmTMB(
  present ~ 1,
  data = dat, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
)

# 04 Predict onto grid cells with at least one survey observation ----
# predict only where the grid has direct support from the data, not
# everywhere the mesh covers
n_col <- 65
x_min <- 450074
y_min <- 6258093
cell_size <- 2000

sampled_col <- floor((dat$x_utm * 1000 - x_min) / cell_size) + 1
sampled_row <- floor((dat$y_utm * 1000 - y_min) / cell_size) + 1
sampled_ids <- unique(sampled_col + (sampled_row - 1) * n_col)

grid_complete <- grid |>
  filter(id %in% sampled_ids) |>
  mutate(x_utm = x / 1000, y_utm = y / 1000)

pred_biomass <- predict(fit_biomass, newdata = grid_complete, type = "response")
pred_presence <- predict(fit_presence, newdata = grid_complete, type = "response")

# every grid cell is kept (NA where there was no survey observation); only
# the columns later scripts need, so this stays a small lookup table rather
# than another copy of the whole grid
grid <- grid |>
  select(id, x, y) |>
  left_join(select(pred_biomass, id, suit_biomass = est), by = "id") |>
  left_join(select(pred_presence, id, suit_presence = est), by = "id")

# 05 Save ----
dir.create(here("data", "sdm", "suitability"), showWarnings = FALSE, recursive = TRUE)
saveRDS(fit_biomass, here("data", "sdm", "suitability", "biomass_space.rds"))
saveRDS(fit_presence, here("data", "sdm", "suitability", "present_space.rds"))
saveRDS(grid, here("data", "grid", "grid_suitability.rds"))

# 06 Plot ----
grid_sampled <- grid |> filter(!is.na(suit_biomass))

ggplot() +
  geom_sf(data = land, fill = "grey80", colour = NA) +
  geom_tile(data = grid_sampled, aes(x, y, fill = suit_biomass), width = 2000, height = 2000) +
  coord_sf(crs = 32632, xlim = range(grid$x), ylim = range(grid$y), expand = FALSE) +
  scale_fill_viridis_c(na.value = NA, trans = "pseudo_log", name = "Suitability\n(biomass)") +
  theme_light()

ggplot() +
  geom_sf(data = land, fill = "grey80", colour = NA) +
  geom_tile(data = grid_sampled, aes(x, y, fill = suit_presence), width = 2000, height = 2000) +
  coord_sf(crs = 32632, xlim = range(grid$x), ylim = range(grid$y), expand = FALSE) +
  scale_fill_viridis_c(na.value = NA, name = "Suitability\n(presence)") +
  theme_light()
