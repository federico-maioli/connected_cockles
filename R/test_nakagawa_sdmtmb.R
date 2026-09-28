# Validation for nakagawa_sdmtmb() (R/helpers.R): does it give the right
# variance components? Checked two ways:
#   01 fit it on sdmTMB's own example data (pcod, yelloweye), covering a
#      spatiotemporal Tweedie field, a Bernoulli model with a year intercept,
#      and a negative binomial (nbinom2) model - so every component the
#      function computes gets exercised at least once.
#   02 refit the same models without the spatial field, so they can also be
#      fitted in lme4 / glmmTMB, and compare variances directly against
#      insight::get_variance() (what performance::r2_nakagawa() itself calls).
# Not part of the pipeline - run this after changing nakagawa_sdmtmb() to
# confirm it still agrees with the reference implementation.

library(tidyverse)
library(here)
library(sdmTMB)

source(here("R", "helpers.R"))

# our variances next to insight's for the same model fitted in lme4 / glmmTMB
compare_insight <- function(fit_sdm, fit_ref) {
  mine <- nakagawa_sdmtmb(fit_sdm)
  ref <- insight::get_variance(fit_ref)
  tibble(
    component = c("fixed", "random_intercept", "distribution"),
    sdmTMB = mine$variance[match(component, mine$component)],
    insight = c(ref$var.fixed, ref$var.random, ref$var.distribution),
    rel_diff = abs(sdmTMB - insight) / insight
  )
}

# 01 Test on the sdmTMB example data ----
# pcod: Tweedie for density (with a spatiotemporal field) and Bernoulli for
# presence (with a year intercept); yelloweye: negative binomial for catch counts
pcod <- pcod |> mutate(fyear = factor(year))
yelloweye <- yelloweye |>
  mutate(depth_c = as.numeric(scale(log(depth))), fyear = factor(year))

mesh_pcod <- make_mesh(pcod, c("X", "Y"), cutoff = 10)
mesh_ye <- make_mesh(yelloweye, c("X", "Y"), cutoff = 10)

fit_tweedie <- sdmTMB(
  density ~ depth_scaled + depth_scaled2,
  data = pcod, mesh = mesh_pcod, family = tweedie(),
  spatial = "on", spatiotemporal = "iid", time = "year"
)
fit_bernoulli <- sdmTMB(
  present ~ depth_scaled + depth_scaled2 + (1 | fyear),
  data = pcod, mesh = mesh_pcod, family = binomial(),
  spatial = "on", spatiotemporal = "off"
)
fit_nbinom <- sdmTMB(
  catch_count ~ depth_c + I(depth_c^2),
  data = yelloweye, mesh = mesh_ye, family = nbinom2(),
  spatial = "on", spatiotemporal = "off"
)

test_fits <- list(tweedie = fit_tweedie, bernoulli = fit_bernoulli, nbinom2 = fit_nbinom)
stopifnot(all(map_lgl(test_fits, \(f) sanity(f, silent = TRUE)$all_ok)))

test_fits |>
  map(nakagawa_sdmtmb) |>
  list_rbind(names_to = "model") |>
  filter(variance > 0) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

# 02 Check against insight ----
# same models without the spatial fields, plus a year intercept, so lme4 and
# glmmTMB can fit them too; insight::get_variance() is what r2_nakagawa() calls
checks <- list(
  bernoulli = compare_insight(
    sdmTMB(present ~ depth_scaled + depth_scaled2 + (1 | fyear),
      data = pcod, family = binomial(), spatial = "off"
    ),
    lme4::glmer(present ~ depth_scaled + depth_scaled2 + (1 | fyear),
      data = pcod, family = binomial()
    )
  ),
  nbinom2 = compare_insight(
    sdmTMB(catch_count ~ depth_c + I(depth_c^2) + (1 | fyear),
      data = yelloweye, family = nbinom2(), spatial = "off"
    ),
    suppressWarnings(lme4::glmer.nb(catch_count ~ depth_c + I(depth_c^2) + (1 | fyear),
      data = yelloweye
    ))
  ),
  tweedie = compare_insight(
    sdmTMB(density ~ depth_scaled + depth_scaled2 + (1 | fyear),
      data = pcod, family = tweedie(), spatial = "off"
    ),
    glmmTMB::glmmTMB(density ~ depth_scaled + depth_scaled2 + (1 | fyear),
      data = pcod, family = glmmTMB::tweedie()
    )
  )
)

checks |>
  list_rbind(names_to = "model") |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  print(n = Inf)
stopifnot(all(map_lgl(checks, \(x) all(x$rel_diff < 0.02))))

cat("\nnakagawa_sdmtmb() matches insight::get_variance() within 2% on every component.\n")
