# Attach the grid connectivity metrics to the cockle survey observations.
#
# The weighted connectivity lives on the 2 km grid (one value per cell); each
# survey point falls inside one grid cell. The grid is a regular UTM (EPSG:32632)
# raster, so we rasterise the connectivity columns and sample them at the points.
#
# Two outputs:
#   cockles_connectivity.rds              pooled ("all") run - the main dataset
#   cockles_connectivity_sensitivity.rds  every run stacked long, one block per
#                                         year plus "all", for the dispersal
#                                         sensitivity check
# The biomass weights are the same in every run (a time-averaged prediction), so
# the year-specific blocks differ only in the flow matrix used.

library(tidyverse)
library(here)
library(terra)

# rasterise one run's connectivity grid and sample it at the survey points
extract_conn <- function(run, pts, conn_vars) {
  grid <- readRDS(here("data", "derived", paste0("connectivity_weighted_", run, ".rds")))
  conn_rast <- rast(
    as.data.frame(grid[, c("x_utm", "y_utm", conn_vars)]),
    type = "xyz",
    crs = "EPSG:32632"
  )
  terra::extract(conn_rast, pts, ID = FALSE) |>
    rename(grid_id = id)
}

# 01 Load inputs ----
dat <- readRDS(here("data", "derived", "cockles_clean.rds"))

# connectivity metrics to attach to each observation
# the cell id comes along too, so each observation can be traced back to its
# grid cell (and points sharing a cell can be identified)
conn_vars <- c(
  "id",
  "biomass_in_strength", "log_biomass_in_strength", "presence_in_strength",
  "deg_in", "in_strength", "local_retention",
  "eigen_centrality", "closeness_centrality"
)
runs <- c(as.character(2010:2016), "all")

# 02 Survey points ----
# cockles_clean coordinates are in km -> back to metres to match the grid
pts <- dat |>
  mutate(x_m = x_utm * 1000, y_m = y_utm * 1000) |>
  vect(geom = c("x_m", "y_m"), crs = "EPSG:32632")

# 03 Main dataset from the pooled run ----
dat_all <- bind_cols(dat, extract_conn("all", pts, conn_vars))

# points falling outside the grid get no values
cat("observations:", nrow(dat_all), "| outside the grid:", sum(is.na(dat_all$grid_id)), "\n")
cat("grid cells hit:", n_distinct(dat_all$grid_id, na.rm = TRUE), "\n")

saveRDS(dat_all, here("data", "derived", "cockles_connectivity.rds"))

# 04 Sensitivity dataset: every run stacked ----
# one block of observations per run, tagged by `run`, so a model can be refitted
# per dispersal year with filter(run == "2012")
sens <- runs |>
  set_names() |>
  map(\(r) bind_cols(dat, extract_conn(r, pts, conn_vars))) |>
  list_rbind(names_to = "run")

cat("\nsensitivity rows:", nrow(sens), "=", n_distinct(sens$run), "runs x", nrow(dat), "observations\n")

saveRDS(sens, here("data", "derived", "cockles_connectivity_sensitivity.rds"))

# 05 Spread across runs ----
# how much does each metric move at the survey points when the dispersal year changes?
sens |>
  filter(run != "all") |>
  pivot_longer(all_of(setdiff(conn_vars, "id")), names_to = "metric") |>
  summarise(mean = mean(value), sd = sd(value), .by = c(metric, grid_id)) |>
  summarise(median_cv = median(sd / mean, na.rm = TRUE), .by = metric) |>
  arrange(median_cv) |>
  print()
