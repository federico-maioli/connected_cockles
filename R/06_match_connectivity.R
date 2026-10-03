# Match the connectivity metrics (05_weight_connectivity.R) onto the cleaned
# survey data and the connectivity grid, by extracting each metric's raster
# at each point. The rasters share the grid's CRS and cell size, so
# extracting at a grid cell's own centre just recovers that cell's own value.

library(tidyverse)
library(here)
library(terra)

# 01 Load data ----
dat <- readRDS(here("data", "cockles", "derived", "cockles_env.rds"))
grid <- readRDS(here("data", "grid", "grid_env.rds"))

conn_dir <- here("data", "connectivity", "derived")
conn_files <- list.files(conn_dir, pattern = "\\.asc$", full.names = TRUE)
conn_names <- paste0("conn_", tools::file_path_sans_ext(basename(conn_files)))

# 02 Extract onto the survey data ----
# dat's x_utm/y_utm are in km (02_clean_cockles.R); convert to m to match the
# rasters, same CRS (UTM 32N) as the grid
pts_dat <- vect(cbind(dat$x_utm * 1000, dat$y_utm * 1000), crs = "EPSG:32632")

for (i in seq_along(conn_files)) {
  r <- rast(conn_files[i])
  dat[[conn_names[i]]] <- extract(r, pts_dat, ID = FALSE)[, 1]
}

# 03 Extract onto the grid ----
pts_grid <- vect(grid, geom = c("x", "y"), crs = "EPSG:32632")

for (i in seq_along(conn_files)) {
  r <- rast(conn_files[i])
  grid[[conn_names[i]]] <- extract(r, pts_grid, ID = FALSE)[, 1]
}

# 04 Checks ----
glimpse(dat)
glimpse(grid)

# 05 Save ----
saveRDS(dat, here("data", "cockles", "derived", "cockles_env_conn.rds"))
saveRDS(grid, here("data", "grid", "grid_env_conn.rds"))
