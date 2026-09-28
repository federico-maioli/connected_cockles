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
conn_raw <- read_csv(here("data", "connectivity", "ConMetrics_PresAbs.csv"), show_col_types = FALSE)


# 02 Extract depth from the depth model ----
# depth raster (EPSG:3034, positive metres) is in a different CRS than the
# survey points; reproject the points to the raster and sample depth at each one
pts <- project(vect(dat_raw), depth_rast)
dat_raw$depth_env <- extract(depth_rast, pts, ID = FALSE)[, 1]

# 03 Extract connectivity metrics from the grid ----
# connectivity grid (2 km cells, EPSG:32632) shares the survey CRS; rasterise
# the metrics of interest and sample them at each survey point
conn_vars <- c(
  in_degree = "PresAbs_degin",
  in_strength = "PresAbs_strin",
  in_closeness = "PresAbs_clsin",
  eigen = "PresAbs_eig",
  transitivity = "PresAbs_tr"
)
conn_rast <- rast(as.data.frame(conn_raw)[, c("x", "y", conn_vars)], type = "xyz", crs = "EPSG:32632")
conn_pts <- project(vect(dat_raw), conn_rast)
conn_ext <- extract(conn_rast, conn_pts, ID = FALSE) |>
  rename_with(~ paste0("conn_", names(conn_vars)))
dat_raw <- bind_cols(dat_raw, conn_ext)

# 04 Clean and recode ----
conn_vars_clean <- make_clean_names(names(conn_ext))

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

# 05 Rename and select ----
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
    all_of(conn_vars_clean), # connectivity metrics (ConMetrics_PresAbs)
  )

# 06 Checks ----
glimpse(dat)

# 07 Save ----
out_path <- here("data", "derived", "cockles_clean.rds")
saveRDS(dat, out_path)
