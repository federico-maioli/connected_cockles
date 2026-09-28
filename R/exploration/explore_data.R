# Explore the cleaned cockle data and check the spatiotemporal structure
# needed for the SDM in 02_fit_sdm.R.
#
# Three questions:
#   - is a shared spatial field enough, or does presence/biomass shift
#     enough from year to year to need a spatiotemporal field (iid vs ar1)?
#   - Commercial and Stock are two different survey protocols covering
#     different (mostly non-overlapping) years - does the model need a
#     survey term?
#   - once we settle on a single shared spatial field, does adding year as a
#     fixed factor help (a real year signal) or complicate (confounded with
#     survey coverage, since 2019/2020 are Commercial-only)?
# Mesh and barrier setup follow 02_fit_sdm.R so the two scripts stay comparable.

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)

# compatibility shim: sdmTMBextra 0.0.5 calls fmesher::fm_identical_CRS, renamed
# to fm_crs_is_identical in fmesher >= 0.7; register the old name as an alias.
# required for add_barrier_mesh() below - do not remove unless sdmTMBextra is updated
local({
  exps <- .getNamespaceInfo(asNamespace("fmesher"), "exports")
  if (!exists("fm_identical_CRS", envir = exps, inherits = FALSE)) {
    assign("fm_identical_CRS", "fm_crs_is_identical", envir = exps)
  }
})

# 01 Load data ----
dat <- readRDS(here("data", "derived", "cockles_clean.rds")) |>
  mutate(
    survey = factor(survey),
    year_i = as.integer(as.character(year)) # sdmTMB needs a numeric time column for AR(1)
  )

count(dat, survey, year)

# presence rate by year/survey - flags years with zero presences (complete
# separation risk for `year` as a factor, see section 06)
dat |>
  count(survey, year, present) |>
  pivot_wider(names_from = present, values_from = n, values_fill = 0, names_prefix = "present_") |>
  print(n = Inf)

# 02 Plot coverage and presence by year ----
ggplot(dat, aes(x_utm, y_utm, colour = survey)) +
  geom_point(size = 0.6) +
  coord_equal() +
  facet_wrap(~year) +
  theme_light() +
  theme(strip.background = element_blank())
ggsave(here("output", "figs", "explore_survey_coverage.png"), width = 9, height = 7)

ggplot(dat, aes(x_utm, y_utm, colour = factor(present))) +
  geom_point(size = 0.6) +
  coord_equal() +
  facet_wrap(~year) +
  scale_colour_manual(values = c("0" = "grey80", "1" = "firebrick"), name = "present") +
  theme_light() +
  theme(strip.background = element_blank())
ggsave(here("output", "figs", "explore_presence_by_year.png"), width = 9, height = 7)

# 03 Coastline barrier mesh ----
# same as 02_fit_sdm.R: local coastline (UTM 32N, metres), cropped to the
# survey extent with a 30 km margin, fmesher mesh with a capped max edge,
# and a land barrier so correlation does not leak across headlands
land_utm <- st_read(
  here("data", "raw", "boundaries", "land_small_utm", "land_small_utm.shp"),
  quiet = TRUE
) |>
  st_make_valid()

margin <- 30000
region <- st_bbox(
  c(
    xmin = min(dat$x_utm) * 1000 - margin,
    ymin = min(dat$y_utm) * 1000 - margin,
    xmax = max(dat$x_utm) * 1000 + margin,
    ymax = max(dat$y_utm) * 1000 + margin
  ),
  crs = st_crs(32632)
)
land_region <- suppressWarnings(st_crop(land_utm, region))

inla_mesh <- fmesher::fm_mesh_2d_inla(
  loc = cbind(dat$x_utm, dat$y_utm),
  max.edge = c(2, 10),
  offset = c(5, 20),
  cutoff = 1
)
mesh <- make_mesh(dat, c("x_utm", "y_utm"), mesh = inla_mesh)
barrier_mesh <- sdmTMBextra::add_barrier_mesh(
  mesh, land_region,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)

# 04 Spatial vs spatiotemporal field ----
# same fixed effects (survey + year) throughout, so AIC differences come only
# from the field: a shared spatial field, or one that also varies by year
st_specs <- tribble(
  ~id, ~spatiotemporal,
  "spatial_only", "off",
  "spatiotemporal_iid", "iid",
  "spatiotemporal_ar1", "ar1"
)

st_fits <- list()
st_comparison <- tibble()
for (i in seq_len(nrow(st_specs))) {
  id <- st_specs$id[i]
  fit <- sdmTMB(
    present ~ survey + year,
    data = dat,
    mesh = barrier_mesh,
    time = "year_i",
    spatial = "on",
    spatiotemporal = st_specs$spatiotemporal[i],
    family = binomial(link = "logit")
  )
  st_fits[[id]] <- fit
  cat("\n--", id, "--\n")
  ok <- sanity(fit) # prints its own checklist, e.g. sigma_O collapsing to ~0
  st_comparison <- bind_rows(st_comparison, tibble(
    id = id,
    converged = isTRUE(ok$all_ok),
    aic = AIC(fit)
  ))
}

cat("\nSpatial vs spatiotemporal field:\n")
print(st_comparison |> arrange(aic))

# from here on: one shared spatial field only (spatial = "on", spatiotemporal
# = "off") - no persistent field vs. year-varying field tradeoff to worry about

# 05 Does the model need a survey term? ----
fit_survey_year <- st_fits[["spatial_only"]] # already fit as survey + year above
fit_year_only <- sdmTMB(
  present ~ year,
  data = dat,
  mesh = barrier_mesh,
  spatial = "on",
  family = binomial(link = "logit")
)

cat("\nSurvey term check (shared spatial field):\n")
tibble(
  model = c("survey + year", "year only"),
  aic = c(AIC(fit_survey_year), AIC(fit_year_only))
) |> print()

cat("\nsurvey coefficient:\n")
print(tidy(fit_survey_year, effects = "fixed", conf.int = TRUE) |> filter(str_starts(term, "survey")))

# 06 Does year as a factor help or complicate? ----
# AIC-wise year helps a lot, but 2019, 2020 and 2021-Commercial have ZERO
# presences (see the table in section 01) - with survey already in the model,
# year2019/year2020 are trying to estimate a within-Commercial log-odds from
# all-zero data, which is unidentifiable: complete separation, not just a
# large SE (estimate ~-20, SE ~2900 below). Fix is to drop those two
# all-zero years, not to keep year as a plain factor.
fit_survey_only <- sdmTMB(
  present ~ survey,
  data = dat,
  mesh = barrier_mesh,
  spatial = "on",
  family = binomial(link = "logit")
)

cat("\nYear term check (shared spatial field):\n")
tibble(
  model = c("survey + year", "survey only"),
  aic = c(AIC(fit_survey_year), AIC(fit_survey_only))
) |> print()

cat("\nyear coefficients (2019/2020 ~ -20 with SE in the thousands = complete separation, not a real estimate):\n")
print(tidy(fit_survey_year, effects = "fixed", conf.int = TRUE) |> filter(str_starts(term, "year")))

# refit dropping the two all-zero years, to see the year effect without the
# separation problem
dat_nosep <- dat |> filter(!year %in% c("2019", "2020"))
mesh_nosep <- make_mesh(dat_nosep, c("x_utm", "y_utm"), mesh = inla_mesh)
barrier_mesh_nosep <- sdmTMBextra::add_barrier_mesh(
  mesh_nosep, land_region,
  range_fraction = 0.1,
  proj_scaling = 1000,
  plot = FALSE
)
fit_year_nosep <- sdmTMB(
  present ~ survey + year,
  data = dat_nosep,
  mesh = barrier_mesh_nosep,
  spatial = "on",
  family = binomial(link = "logit")
)

cat("\nyear coefficients after dropping 2019/2020 (all-zero years):\n")
print(tidy(fit_year_nosep, effects = "fixed", conf.int = TRUE) |> filter(str_starts(term, "year")))
