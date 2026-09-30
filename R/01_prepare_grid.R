# Build the regular 2 x 2 km grid used by the MIKE/ABM larval transport model
# (65 columns x 36 rows, UTM 32N metres) - the grid the raw MIKE connectivity
# matrix is defined on, and the one every later script predicts or matches onto.
#
# Ported from Flemming Thorbjørn Hansen's original connectivity-analysis
# script (2024-06-01) - only the grid definition. Cell id runs left-to-right
# within a row, then bottom-to-top (row 1 = southernmost row), matching the
# row/column order of the connectivity matrix: id = col + (row - 1) * n_col.
# Verified against ConMetrics_PresAbs.csv: this formula reproduces its x/y for
# every one of the 2,340 cells exactly (0 m mismatch).
#
# Also extracts the same environmental covariates as 02_clean_cockles.R, with
# the same column names, onto this grid, from the DTU Aqua raster set
# (data/env/Asc4dtuaqua).

library(tidyverse)
library(here)
library(terra)
library(sf)

# 01 Grid definition ----
n_col <- 65
n_row <- 36
cell_size <- 2000
x_min <- 450074
x_max <- 580074
y_min <- 6258093
y_max <- 6330093

stopifnot(
  (x_max - x_min) / n_col == cell_size,
  (y_max - y_min) / n_row == cell_size
)

# 02 Cell centres ----
grid <- expand_grid(row = 1:n_row, col = 1:n_col) |>
  mutate(
    id = col + (row - 1) * n_col,
    x = x_min + (col - 0.5) * cell_size,
    y = y_min + (row - 0.5) * cell_size
  ) |>
  select(id, row, col, x, y) |>
  arrange(id)

# 03 Environmental covariates ----
# raster CRS (EPSG:3034) differs from the grid (UTM 32N), so the points are
# projected to each raster before extracting - same rasters and column names
# as 02_clean_cockles.R
pts <- vect(grid, geom = c("x", "y"), crs = "EPSG:32632")

env_dir <- here("data", "env", "Asc4dtuaqua")
env_files <- c(
  depth = "depth.asc",
  temp = "temp_bot_mean.asc",
  sal = "salt_bot_mean.asc",
  oxy = "do4mgl_bot_mean.asc",
  chla = "chla_bot_mean.asc",
  shear_max = "tw_bot_mean_of_max.asc",
  shear_mean = "tw_bot_mean.asc",
  phyto = "pc_bot_stdev.asc",
  sediment = "sediment.asc"
)

for (nm in names(env_files)) {
  r <- rast(file.path(env_dir, env_files[nm]))
  pts_proj <- project(pts, r)
  grid[[nm]] <- extract(r, pts_proj, ID = FALSE)[, 1]
}

# depth.asc is negative below sea level; flip so larger values mean deeper
grid$depth <- -grid$depth

# 04 Checks ----
glimpse(grid)

# 05 Save ----
dir.create(here("data", "grid"), showWarnings = FALSE)
saveRDS(grid, here("data", "grid", "grid_env.rds"))

# 06 Plot ----
land <- st_read(
  here("data", "boundaries", "land_small_utm", "land_small_utm.shp"),
  quiet = TRUE
) |>
  st_make_valid()

# standardized (z-score) so all nine covariates, on very different units,
# can share one fill scale; every facet keeps the same map extent, so
# facet_wrap() doesn't need (and can't combine with coord_sf() and) free scales
grid_long <- grid |>
  mutate(across(all_of(names(env_files)), ~ as.numeric(scale(.x)))) |>
  pivot_longer(all_of(names(env_files)), names_to = "covariate", values_to = "value")

ggplot() +
  geom_sf(data = land, fill = "grey80", colour = NA) +
  geom_tile(data = grid_long, aes(x, y, fill = value), width = cell_size, height = cell_size) +
  coord_sf(crs = 32632, xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE) +
  facet_wrap(~covariate) +
  scale_fill_viridis_c(na.value = NA, name = "Value (standardized)") +
  theme_light() +
  theme(strip.background = element_blank())
