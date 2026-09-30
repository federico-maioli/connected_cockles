# Is the negative in-strength effect driven by a few high-inflow cells?
# The full model (environment + spatial field + presence-weighted in-strength,
# as in 06_fit_sdm.R) is refitted without
#   a) the survey samples in the hub cells of the raw connectivity (cells that
#      are the main destination of many release cells, see
#      plot_connectivity.R) - the two strongest hubs (1919, 1854) have no
#      samples, so this drops the sampled ones among the top five (1259, 865,
#      996)
#   b) the survey samples in the top 10% of survey cells by in-strength
# and the in-strength coefficient is compared with the full-data fit.
# Standardisation is kept from the full data, so coefficients are comparable.

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)

source(here("R", "helpers.R")) # fmesher compatibility shim, needed by add_barrier_mesh() below

hub_cells <- c(1259, 865, 996)

# 01 Data ----
# the model rows of the full fit: complete cases, standardised covariates
dat <- readRDS(here("data", "sdm", "main", "biomass_space_env_conn.rds"))$data |>
  mutate(id = (floor((x_utm * 1000 - 450074) / 2000) + 1) + floor((y_utm * 1000 - 6258093) / 2000) * 65)

cell_instrength <- dat |> distinct(id, conn_in_strength)
top_cells <- cell_instrength |>
  filter(conn_in_strength >= quantile(conn_in_strength, 0.9)) |>
  pull(id)

subsets <- list(
  drop_hubs = filter(dat, !id %in% hub_cells),
  drop_top10 = filter(dat, !id %in% top_cells)
)

land_utm <- st_read(here("data", "boundaries", "land_small_utm", "land_small_utm.shp"), quiet = TRUE) |>
  st_make_valid()

# 02 Refit the full model on each subset ----
# same barrier mesh recipe as 06_fit_sdm.R, rebuilt on each subset's locations
coefs <- map(names(subsets), function(s) {
  d <- subsets[[s]]
  region <- st_bbox(
    c(
      xmin = min(d$x_utm) * 1000 - 30000, ymin = min(d$y_utm) * 1000 - 30000,
      xmax = max(d$x_utm) * 1000 + 30000, ymax = max(d$y_utm) * 1000 + 30000
    ),
    crs = st_crs(32632)
  )
  inla_mesh <- fmesher::fm_mesh_2d_inla(loc = cbind(d$x_utm, d$y_utm), max.edge = c(2, 10), offset = c(5, 20), cutoff = 1)
  barrier_mesh <- sdmTMBextra::add_barrier_mesh(
    make_mesh(d, c("x_utm", "y_utm"), mesh = inla_mesh), suppressWarnings(st_crop(land_utm, region)),
    range_fraction = 0.1, proj_scaling = 1000, plot = FALSE
  )

  fit_biomass <- sdmTMB(
    biomass ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
    data = d, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
  )
  fit_present <- sdmTMB(
    present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std + conn_in_strength_std,
    data = d, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit")
  )

  bind_rows(
    tidy(fit_biomass, conf.int = TRUE) |> mutate(response = "biomass"),
    tidy(fit_present, conf.int = TRUE) |> mutate(response = "present")
  ) |>
    filter(term == "conn_in_strength_std") |>
    mutate(subset = s, n = nrow(d), converged = c(sanity(fit_biomass, silent = TRUE)$all_ok, sanity(fit_present, silent = TRUE)$all_ok))
}) |>
  list_rbind()

# 03 Compare with the full-data fits ----
full <- map(c("biomass", "present"), \(r) {
  tidy(readRDS(here("data", "sdm", "main", paste0(r, "_space_env_conn.rds"))), conf.int = TRUE) |>
    filter(term == "conn_in_strength_std") |>
    mutate(response = r, subset = "full", n = nrow(dat), converged = TRUE)
}) |>
  list_rbind()

bind_rows(full, coefs) |>
  select(response, subset, n, converged, estimate, conf.low, conf.high) |>
  mutate(across(c(estimate, conf.low, conf.high), \(x) round(x, 3))) |>
  arrange(response, factor(subset, levels = c("full", "drop_hubs", "drop_top10"))) |>
  print()
