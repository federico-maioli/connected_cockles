# Attach the grid connectivity metrics to the cockle survey observations.
#
# The weighted connectivity lives on the 2 km grid (one value per cell); each
# survey point falls inside one grid cell. The grid is a regular UTM (EPSG:32632)
# raster, so we rasterise the connectivity columns and sample them at the points.

library(tidyverse)
library(here)
library(terra)

# 01 Load inputs ----
dat <- readRDS(here("data", "final", "cockles_clean.rds"))
grid <- readRDS(here("data", "final", "connectivity_weighted.rds"))

# connectivity metrics to attach to each observation
conn_vars <- c(
  "biomass_supply", "biomass_export",
  "biomass_supply_log", "biomass_export_log",
  "present_supply", "present_export",
  "in_strength", "out_strength", "local_retention", "eigen_centrality"
)

# 02 Rasterise the connectivity grid ----
# regular 2 km grid of cell centres -> multi-layer raster (one layer per metric)
conn_rast <- rast(
  as.data.frame(grid[, c("x_utm", "y_utm", conn_vars)]),
  type = "xyz",
  crs = "EPSG:32632"
)

# 03 Extract at survey points ----
# cockles_clean coordinates are in km -> back to metres to match the grid
pts <- dat |>
  mutate(x_m = x_utm * 1000, y_m = y_utm * 1000) |>
  vect(geom = c("x_m", "y_m"), crs = "EPSG:32632")

conn_vals <- terra::extract(conn_rast, pts, ID = FALSE)

dat <- bind_cols(dat, conn_vals)

# 04 Save ----
# cockles_clean plus the connectivity covariates (kept as a separate file so
# re-running 00 does not clobber it)
saveRDS(dat, here("data", "final", "cockles_connectivity.rds"))
