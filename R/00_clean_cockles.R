library(tidyverse)
library(here)
library(sf)
library(terra)
library(janitor)

# 01 Load data ----
dat_raw <- st_read(
  here("data", "raw", "cockles", "DataFrame4spatstats_UTM_v03.shp"),
  quiet = TRUE
)

depth_rast <- rast(here("data", "raw", "environment", "ddm_50m.dybde.tiff"))

# 02 Extract depth from the depth model ----
# depth raster (EPSG:3034, positive metres) is in a different CRS than the
# survey points; reproject the points to the raster and sample depth at each one
pts <- project(vect(dat_raw), depth_rast)
dat_raw$depth_env <- extract(depth_rast, pts, ID = FALSE)[, 1]

# 03 Clean and recode ----
dat <- dat_raw |>
  st_drop_geometry() |>
  clean_names() |>
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
# rename inside select so the raw columns we want are kept and relabelled
dat <- dat |>
  select(
    survey, year, x_utm, y_utm,
    density, biomass, present,
    depth = depth_env, # sampled from the 50 m depth model

    temp = temp_mean,
    sal = salt_mean,
    oxy = do4mglavg, # days/year with dissolved oxygen < 4 mg/l
    chla = chla_bot_m,
    shear_max = twcavgmax, # max shear stress from waves and currents
    shear_mean = twcavg_all, # mean shear stress from waves and currents
    phyto = pcavgstd, # phytoplankton carbon, standard deviation
    sediment, # sediment class (Emodnet)
    ba_strin, bp_strin, pa_strin # larval in-strength (connectivity)
  )

# 05 Checks ----
glimpse(dat)

# 06 Save ----
out_path <- here("data", "derived", "cockles_clean.rds")
saveRDS(dat, out_path)