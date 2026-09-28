# Is there an environmental variable, currently unused, that predicts cockles
# better than what's already in 02_fit_sdm.R (depth, temp, sal, oxy, shear_max)?
#
# Two pools of candidates:
#   - already extracted into cockles_clean.rds but not in the model: chla,
#     phyto, shear_mean, sediment
#   - in the raw shapefile but never extracted at all: carbon flux to the
#     seabed (cflux, a food-supply proxy for a filter feeder), current speed
#     (hcsp - inferred from the name, not documented in the shapefile
#     metadata), GEBCO bathymetry (a coarser alternative to the depth already
#     used), and water-column ("vertical") chlorophyll, vs. the bottom
#     chlorophyll already extracted
# Several raw columns turned out to be entirely empty (DOM, Blue_musse,
# Fishing_ef, Sediment_t, Env__*) - dead columns in the source file, not usable.
#
# This is a quick univariate screen (correlation with biomass, single-variable
# AIC vs. an intercept-only model), no spatial field - just to see which
# candidates are worth fitting properly in 02_fit_sdm.R. A variable can look
# good here and still turn out to be standing in for space once the field is
# in the model (see figS9_correlation.R, figS10_coeff_stability.R for how
# that played out for connectivity).

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)

# 01 Rebuild the Stock-only dataset with the extra candidate columns ----
# same cleaning as 01_clean_cockles.R (drop Commercial, no survey term needed)
dat_raw <- st_read(
  here("data", "raw", "cockles", "DataFrame4spatstats_UTM_v03.shp"),
  quiet = TRUE
) |>
  st_drop_geometry() |>
  janitor::clean_names() |>
  filter(survey != "Commercial") |>
  mutate(present = if_else(density > 0, 1, 0))

candidates <- c(
  chla = "chla_bot_m", phyto = "pcavgstd", shear_mean = "twcavg_all", sediment = "sediment",
  cflux_avg = "cfluxavgkg", cflux_std = "cflux_stdk",
  current_avg = "hcspavgavg", current_std = "hcspavgstd",
  depth_gebco = "gebco_2020", chla_vert = "chla_ver_m", chla_vert_std = "chla_ver_st"
)

dat <- dat_raw |>
  select(present, biomass, all_of(candidates))

# 02 Univariate screen ----
# correlation with biomass (Spearman - biomass is zero-inflated and skewed),
# and how much a single-predictor GLM beats an intercept-only model (delta
# AIC), for both responses
# screens a column of `dat` by its friendly name (select(all_of(candidates))
# already renamed the raw columns to these names)
screen <- function(dat, var_name) {
  d <- tibble(present = dat$present, biomass = dat$biomass, x = dat[[var_name]])
  present_fit <- glm(present ~ x, data = d, family = binomial())
  present_null <- glm(present ~ 1, data = d, family = binomial())
  biomass_fit <- glm(log1p(biomass) ~ x, data = d, family = gaussian())
  biomass_null <- glm(log1p(biomass) ~ 1, data = d, family = gaussian())
  tibble(
    variable = var_name,
    cor_spearman_biomass = cor(d$x, d$biomass, method = "spearman"),
    delta_aic_present = AIC(present_null) - AIC(present_fit),
    delta_aic_biomass = AIC(biomass_null) - AIC(biomass_fit)
  )
}

results <- map(names(candidates), \(nm) screen(dat, nm)) |>
  list_rbind() |>
  arrange(desc(delta_aic_biomass))

print(results, n = Inf)

# 03 Same screen for the variables already in the model, for comparison ----
# from cockles_clean.rds, not raw: the shapefile's own Depth/depthdk fields are
# dead (all NA, like DOM/Blue_musse/etc.) - the real depth in the pipeline
# comes from the DEM raster extraction in 01_clean_cockles.R
used <- c("depth", "temp", "sal", "oxy", "shear_max")
dat_used <- readRDS(here("data", "derived", "cockles_clean.rds")) |> select(present, biomass, all_of(used))
used_results <- map(used, \(nm) screen(dat_used, nm)) |>
  list_rbind() |>
  arrange(desc(delta_aic_biomass))

cat("\nfor comparison, the variables already in the model:\n")
print(used_results, n = Inf)

# 04 Space + env + connectivity (in-strength), with sediment and chla added ----
# same predictors as 02_fit_sdm.R's main "space_env_conn" model, plus the two
# candidates that beat the univariate screen above: sediment (factor, mostly
# class 1 vs 6) and chla (bottom chlorophyll, standardized like the other
# continuous covariates)
source(here("R", "helpers.R"))

dat_full <- readRDS(here("data", "derived", "cockles_clean.rds")) |>
  mutate(
    across(c(depth, temp, sal, oxy, shear_max, chla, conn_in_strength), ~ as.numeric(scale(.x)), .names = "{.col}_std"),
    sediment = factor(sediment)
  )

barrier_mesh <- build_barrier_mesh(dat_full)

form_present <- present ~ depth_std + temp_std + sal_std + oxy_std + shear_max_std +
  conn_in_strength_std + chla_std + sediment
form_biomass <- update(form_present, biomass ~ .)

fit_present <- sdmTMB(form_present, data = dat_full, mesh = barrier_mesh, spatial = "on", family = binomial(link = "logit"))
fit_biomass <- sdmTMB(form_biomass, data = dat_full, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log"))

cat("\npresent ~ ... + chla + sediment, spatial field on:\n")
print(sanity(fit_present)$all_ok)
print(tidy(fit_present, effects = "fixed", conf.int = TRUE), n = Inf)

cat("\nbiomass ~ ... + chla + sediment, spatial field on:\n")
print(sanity(fit_biomass)$all_ok)
print(tidy(fit_biomass, effects = "fixed", conf.int = TRUE), n = Inf)
