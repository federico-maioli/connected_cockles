# Build the 2 x 2 km grid that matches the ABM connectivity matrix.
#
# 65 x 36 = 2340 cells, UTM zone 32N. Node id (1..2340) is the SAME order as the
# rows/cols of the flow matrix, so anything predicted on this grid (e.g. cockle
# abundance/biomass) lines up cell-for-cell with connectivity.
#
# id is row-major: id 1..65 = bottom row (y = ymin side), increasing left->right,
# then upward. Verified against sites.csv (id 11 -> centre 471074, 6259093).

library(tidyverse)
library(here)
library(terra)

# 01 Grid definition ----
nx <- 65 # columns (mgresx)
ny <- 36 # rows (mgresy)
cell <- 2000 # cell size, metres
xmin <- 450074
ymin <- 6258093
crs_utm <- 32632 # UTM zone 32N (Limfjorden / Denmark)

# 02 Cell index to centre coordinates ----
ids <- seq_len(nx * ny)
col <- ((ids - 1) %% nx) + 1
row <- ((ids - 1) %/% nx) + 1 # row 1 = bottom (ymin side)
xc <- xmin + cell * (col - 1) + cell / 2
yc <- ymin + cell * (row - 1) + cell / 2

grid <- tibble(
  id = ids,
  col = col,
  row = row,
  x_utm = xc,
  y_utm = yc
)

# 03 Sample depth from the depth model ----
# same 50 m depth model used for the survey points (EPSG:3034, positive metres);
# reproject the grid centres to the raster and sample depth at each cell
depth_rast <- rast(here("data", "raw", "environment", "ddm_50m.dybde.tiff"))
pts <- project(vect(grid, geom = c("x_utm", "y_utm"), crs = paste0("EPSG:", crs_utm)), depth_rast)
grid$depth <- extract(depth_rast, pts, ID = FALSE)[, 1]

# 04 Check grid matches flow matrix ----
# one grid cell per flow-matrix row/column so predictions map 1:1 onto connectivity
flow <- readRDS(here("data", "derived", "flow_matrix_all.rds"))
stopifnot(nrow(grid) == nrow(flow), nrow(grid) == ncol(flow))

# 05 Save ----
saveRDS(grid, here("data", "derived", "spatial_grid.rds"))
