# Supplementary figure: stability of the in-strength connectivity coefficient
# across model specifications, as a heatmap.
#
# Rows (environment block):
#   full model         all 5 environment terms + connectivity
#   - <one variable>    full model, minus one environment term (x5)
#   connectivity only  connectivity, no environment at all
# Columns: spatial field on vs off (as in figS5_space_confounding.R, but here
# crossed with the leave-one-out environment specs too). Mesh and shim come
# from R/helpers.R (build_barrier_mesh()).
# Cell fill = estimate; a dagger marks a cell whose 95% CI crosses zero. If
# colour barely changes across the grid, the coefficient is not being driven
# by any single covariate, or by whether the field soaks up spatial structure.

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)

source(here("R", "helpers.R"))

# 01 Load data ----
dat <- readRDS(here("data", "derived", "cockles_clean.rds"))

env_vars <- c("depth", "temp", "sal", "oxy", "shear_max")
env_labs <- c(depth = "Depth", temp = "Temperature", sal = "Salinity", oxy = "Oxygen", shear_max = "Shear stress")

dat <- dat |>
  mutate(across(all_of(c(env_vars, "conn_in_strength")), ~ as.numeric(scale(.x)), .names = "{.col}_std"))

env_terms <- paste0(env_vars, "_std")
conn_term <- "conn_in_strength_std"

# same complete-case set for every spec, so they are all fitted on the same data
dat <- dat |>
  filter(
    !is.na(biomass), !is.na(present), !is.na(x_utm), !is.na(y_utm),
    if_all(all_of(c(env_terms, conn_term)), ~ !is.na(.x))
  )

# 02 Coastline barrier mesh ----
barrier_mesh <- build_barrier_mesh(dat)

# 03 Specifications ----
# one row per spec: an id, a plot label, and the predictors on the RHS
# (connectivity is in every spec); displayed full model -> drop one at a time
# -> connectivity alone, so the environment block shrinks down the plot
specs <- bind_rows(
  tibble(id = "full", label = "Full model", predictors = list(c(env_terms, conn_term))),
  map(env_vars, function(v) {
    tibble(
      id = paste0("drop_", v),
      label = paste("–", env_labs[[v]]),
      predictors = list(c(setdiff(env_terms, paste0(v, "_std")), conn_term))
    )
  }) |> list_rbind(),
  tibble(id = "conn_only", label = "Connectivity only (no environment)", predictors = list(conn_term))
)

# 04 Fit every spec, for both responses and with/without the spatial field ----
families <- list(present = binomial(link = "logit"), biomass = tweedie(link = "log"))
spatial_settings <- c("on", "off")

fits <- list()
for (response_name in names(families)) {
  for (spatial in spatial_settings) {
    for (i in seq_len(nrow(specs))) {
      form <- reformulate(specs$predictors[[i]], response = response_name)
      fit_id <- paste(response_name, spatial, specs$id[i], sep = "_")
      fits[[fit_id]] <- tryCatch(
        sdmTMB(form, data = dat, mesh = barrier_mesh, spatial = spatial, family = families[[response_name]]),
        error = function(e) {
          message("failed: ", fit_id, " - ", conditionMessage(e))
          NULL
        }
      )
    }
  }
}

# 05 Pull the connectivity coefficient out of every spec ----
coefs <- tibble()
for (response_name in names(families)) {
  for (spatial in spatial_settings) {
    for (i in seq_len(nrow(specs))) {
      fit <- fits[[paste(response_name, spatial, specs$id[i], sep = "_")]]
      converged <- !is.null(fit) && isTRUE(suppressMessages(sanity(fit))$all_ok)
      est <- if (converged) {
        tidy(fit, effects = "fixed", conf.int = TRUE) |> filter(term == conn_term)
      } else {
        tibble(estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_)
      }
      coefs <- bind_rows(coefs, tibble(
        response = response_name,
        spatial = spatial,
        id = specs$id[i],
        label = specs$label[i],
        converged = converged,
        estimate = est$estimate,
        conf.low = est$conf.low,
        conf.high = est$conf.high
      ))
    }
  }
}

print(coefs, n = Inf)

# 06 Plot ----
spec_order <- rev(specs$label) # full model on top, connectivity-only at the bottom
plot_data <- coefs |>
  filter(converged) |>
  mutate(
    label = factor(label, levels = spec_order),
    spatial_lab = factor(recode(spatial, "on" = "Spatial field: on", "off" = "Spatial field: off"),
      levels = c("Spatial field: on", "Spatial field: off")
    ),
    crosses_zero = conf.low <= 0 & conf.high >= 0,
    response_lab = recode(response, biomass = "Biomass (log link)", present = "Presence (logit link)"),
    cell_label = paste0(sprintf("%.2f", estimate), if_else(crosses_zero, "†", ""))
  )

fill_limit <- max(abs(plot_data$estimate))

p <- ggplot(plot_data, aes(spatial_lab, label, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = cell_label), size = 3.1) +
  facet_wrap(~response_lab) +
  scale_fill_distiller(palette = "RdBu", limits = c(-fill_limit, fill_limit), name = "Estimate") +
  labs(x = NULL, y = NULL, caption = "† 95% CI crosses zero") +
  theme_light(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black"),
    panel.grid = element_blank(),
    plot.caption = element_text(colour = "grey30", hjust = 0)
  )

# 07 Save ----
dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "figS10_coeff_stability.png"), p,
  width = 9, height = 4.4, dpi = 600, bg = "white"
)
