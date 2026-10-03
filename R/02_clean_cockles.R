library(tidyverse)
library(here)
library(sf)
library(terra)
library(janitor)

# KSKV years added to the Stock / HighRes survey data, and the maximum
# distance (m) from a Stock / HighRes station for a KSKV station to be kept,
# so the spatial domain stays that of the cockle surveys (drops the
# south-eastern basins, which only KSKV covers)
kskv_years <- 2018:2023
kskv_max_dist <- 4000

# 01 Load data ----
# cockle surveys (Stock, HighRes; UTM 32N, metres)
cockles <- st_read(
  here("data", "cockles", "raw", "DataFrame4spatstats_UTM_v03.shp"),
  quiet = TRUE
) |>
  st_set_crs(32632) |>
  clean_names() |>
  # Commercial survey is not usable (density not recorded), drop it
  filter(survey != "Commercial") |>
  transmute(survey, year, density, biomass)

# KSKV blue mussel and oyster dredge surveys: one row per station and species
# caught; a row with another species code has zero cockle weight. Collapsed to
# one row per station: cockle (HMS) biomass and density summed (zero if none),
# biomass from kg/m2 to g/m2, position = midpoint of the dredge track
kskv <- read_csv(here("data", "cockles", "raw", "KSKV_HMS_2010-2026.csv"), show_col_types = FALSE) |>
  filter(year %in% kskv_years) |>
  summarise(
    lat = mean((latPosStartDec + latPosEndDec) / 2),
    lon = mean((lonPosStartDec + lonPosEndDec) / 2),
    biomass = sum(`Biomass kgm2`[speciesCode == "HMS"]) * 1000,
    density = sum(`Density #m2`[speciesCode == "HMS"]),
    .by = c(year, cruise, station)
  ) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(32632) |>
  transmute(survey = "KSKV", year, density, biomass)

near_cockle_station <- apply(st_distance(kskv, cockles), 1, min) <= kskv_max_dist
kskv <- kskv[near_cockle_station, ]

dat_raw <- bind_rows(cockles, kskv)

# 02 Extract environmental covariates ----
# same DTU Aqua rasters and column names as 01_prepare_grid.R, so survey
# points and the connectivity grid use identical covariate values
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
  dat_raw[[nm]] <- extract(r, pts_proj, ID = FALSE)[, 1]
}

# depth.asc is negative below sea level; flip so larger values mean deeper
dat_raw$depth <- -dat_raw$depth

# 03 Clean and recode ----
dat <- dat_raw |>
  mutate(
    # utm coordinates in km
    x_utm = st_coordinates(geometry)[, 1] / 1000,
    y_utm = st_coordinates(geometry)[, 2] / 1000
  ) |>
  st_drop_geometry() |>
  mutate(
    # surveys with their own gear and protocol. Stock (grab, 2021-2025) is the
    # reference; HighRes (2021-2022) uses the same grab protocol, so it joins
    # it. The 2018 Stock survey used a suction dredge over a large area, so it
    # is its own level
    survey = case_when(
      survey == "Stock" & year == 2018 ~ "Stock2018",
      survey == "HighRes" ~ "Stock",
      .default = survey
    ),
    survey = factor(survey, levels = c("Stock", "Stock2018", "KSKV")),
    # presence / absence from density
    present = if_else(density > 0, 1, 0),
    year = as.factor(year)
  ) |>
  # connectivity metrics are matched on separately, in 06_match_connectivity.R
  select(
    survey, year, x_utm, y_utm,
    density, biomass, present,
    depth, temp, sal, oxy, chla, shear_max, shear_mean, phyto, sediment
  )

# 04 Checks ----
glimpse(dat)
count(dat, survey, year)

# 05 Save ----
dir.create(here("data", "cockles", "derived"), showWarnings = FALSE)
saveRDS(dat, here("data", "cockles", "derived", "cockles_env.rds"))
