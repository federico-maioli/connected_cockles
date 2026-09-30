library(tidyverse)
library(here)
library(sf)
library(terra)
library(janitor)

# 01 Load data ----
dat_raw <- st_read(
  here("data", "cockles", "raw", "DataFrame4spatstats_UTM_v03.shp"),
  quiet = TRUE
)

# 02 Extract environmental covariates ----
# same DTU Aqua rasters and column names as 01_prepare_grid.R, so survey
# points and the connectivity grid use identical covariate values. Extracted
# columns are suffixed "_env" so they don't collide with the raw shapefile's
# own columns of the same name (e.g. "Depth", which is entirely empty) -
# renamed to their final names in the select() below
pts <- vect(dat_raw)

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
  dat_raw[[paste0(nm, "_env")]] <- extract(r, pts_proj, ID = FALSE)[, 1]
}

# depth.asc is negative below sea level; flip so larger values mean deeper
dat_raw$depth_env <- -dat_raw$depth_env

# 03 Clean and recode ----
dat <- dat_raw |>
  st_drop_geometry() |>
  clean_names() |>
  # Commercial survey is not usable (density not recorded), drop it
  filter(survey != "Commercial") |>
  mutate(
    # utm coordinates from m to km
    x_utm = x_utm / 1000,
    y_utm = y_utm / 1000,
    # "HighRes" and "Stock" surveys share the same protocol -> one stratum
    survey = if_else(survey == "HighRes", "Stock", survey),
    # presence / absence from density
    present = if_else(density > 0, 1, 0),
    year = as.factor(year)
  )

# 04 Rename and select ----
# connectivity metrics are matched on separately, in 05_match_connectivity.R
dat <- dat |>
  select(
    survey, year, x_utm, y_utm,
    density, biomass, present,
    depth = depth_env,
    temp = temp_env,
    sal = sal_env,
    oxy = oxy_env,
    chla = chla_env,
    shear_max = shear_max_env,
    shear_mean = shear_mean_env,
    phyto = phyto_env,
    sediment = sediment_env,
  )

# 05 Checks ----
glimpse(dat)

# 06 Save ----
dir.create(here("data", "cockles", "derived"), showWarnings = FALSE)
saveRDS(dat, here("data", "cockles", "derived", "cockles_env.rds"))
  