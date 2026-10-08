# All supplementary figures except S7 (predictions, figS7_predictions.R), one
# section each, numbered in the order they are cited in main.tex:
#   S1   survey stations by year and survey
#   S2   correlation among all covariates (environment and connectivity)
#   S3   connectivity metrics on the 2 km grid cells with stations
#   S4   the barrier mesh
#   S5   cross-validation folds, occurrence and positive biomass
#   S6   simulation-based residuals of the full model
#   S8   connectivity coefficient of each metric in the full model
# (The covariate maps are in the main text, fig2_predictors.R.)
#
# Connectivity metrics are the presence-weighted ones used in 07_fit_sdm.R;
# every model is a delta-gamma model read from its saved fit (data/sdm/main
# and data/sdm/sensitivity), nothing is refitted here. Model results are shown
# per component: 1 presence (logit link), 2 biomass where present (log link).

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
library(patchwork)
library(ggstats)

conn_labs <- c(
  conn_in_degree = "In-degree",
  conn_in_strength = "In-strength",
  conn_in_closeness = "In-closeness",
  conn_eigen = "Eigenvector centrality",
  conn_transitivity = "Transitivity"
)

# fitted model saved by 07_fit_sdm.R: the in-strength models live in main/,
# the full model with the other metrics in sensitivity/
read_fit <- function(folder, model) {
  readRDS(here("data", "sdm", folder, paste0(model, ".rds")))
}

# delta-gamma components, in model order (component 1, component 2)
part_labs <- c(presence = "Presence (logit link)", biomass = "Biomass where present (log link)")
facet_theme <- theme(
  strip.background = element_blank(),
  strip.text = element_text(colour = "black", size = 11)
)
tag_theme <- theme(plot.tag = element_text(size = 12, face = "bold", colour = "grey20"))

# map theme for the grid maps
map_theme <- theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 9.5, hjust = 0.5, colour = "grey20"),
    legend.position = "right",
    legend.key.width = unit(3, "mm"), legend.key.height = unit(8, "mm"),
    legend.text = element_text(size = 7, colour = "grey30"),
    plot.margin = margin(4, 4, 4, 4)
  )

# centre of the 2 x 2 km connectivity grid cell holding a station (km in, m out)
cell_centre_x <- function(x_utm) 450074 + floor((x_utm * 1000 - 450074) / 2000) * 2000 + 1000
cell_centre_y <- function(y_utm) 6258093 + floor((y_utm * 1000 - 6258093) / 2000) * 2000 + 1000

# 01 Shared data ----
water <- st_read(here("data", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

# the model dataset (complete cases, standardised covariates, presence-weighted
# connectivity), taken from a saved fit so it matches the models exactly
fit_full <- read_fit("main", "space_env_conn")
dat_model <- fit_full$data

# grid cells holding at least one survey station, and a common map extent
sampled_cells <- dat_model |>
  transmute(x = cell_centre_x(x_utm), y = cell_centre_y(y_utm)) |>
  distinct()
xlim_map <- range(dat_model$x_utm * 1000) + c(-6000, 6000)
ylim_map <- range(dat_model$y_utm * 1000) + c(-6000, 6000)

# white mask = the plot frame minus the fjord, drawn over the tiles so the 2 km
# cells are clipped to the coastline instead of spilling onto land
frame <- st_as_sfc(st_bbox(
  c(xmin = xlim_map[1], ymin = ylim_map[1], xmax = xlim_map[2], ymax = ylim_map[2]),
  crs = st_crs(32632)
))
land_mask <- st_difference(frame, st_union(water))

dir.create(here("output", "figs", "supp"), showWarnings = FALSE, recursive = TRUE)

# 02 Figure S1: survey stations by year ----
# every station used in the models, one panel per year (three per row),
# coloured by survey
survey_labs <- c(Stock = "Grab (2021-2025)", Stock2018 = "Suction dredge (2018)", KSKV = "KSKV dredge")
stations_s1 <- dat_model |>
  mutate(survey = factor(survey_labs[as.character(survey)], levels = survey_labs))

fig_s1 <- ggplot() +
  geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.25) +
  geom_point(data = stations_s1, aes(x_utm * 1000, y_utm * 1000, colour = survey), size = 0.5, alpha = 0.6) +
  facet_wrap(~year, ncol = 3) +
  scale_colour_manual(values = c("#1F4E79", "#E69F00", "#CC79A7"), name = NULL) +
  coord_sf(crs = 32632, xlim = xlim_map, ylim = ylim_map, expand = FALSE) +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(colour = "grey15", margin = margin(3, 0, 3, 0)),
    legend.position = "bottom"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
ggsave(here("output", "figs", "supp", "figS1_survey_stations.png"), fig_s1, width = 11, height = 11, dpi = 600, bg = "white")

# 03 Figure S2: correlation among all covariates ----
# the five environmental predictors and the five connectivity metrics, at the
# survey stations; Spearman, since the connectivity metrics are skewed
s2_vars <- c(
  Depth = "depth", Temperature = "temp", Salinity = "sal", Oxygen = "oxy", `Shear stress` = "shear_max",
  `In-degree` = "conn_in_degree", `In-strength` = "conn_in_strength", `In-closeness` = "conn_in_closeness",
  `Eigenvector centrality` = "conn_eigen", Transitivity = "conn_transitivity"
)

cor_s2 <- dat_model |>
  select(all_of(s2_vars)) |>
  cor(method = "spearman", use = "complete.obs") |>
  as_tibble(rownames = "var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "rho") |>
  mutate(
    var1 = factor(var1, levels = names(s2_vars)),
    var2 = factor(var2, levels = rev(names(s2_vars)))
  )

fig_s2 <- ggplot(cor_s2, aes(var1, var2, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.8) +
  # separates environment from connectivity
  geom_vline(xintercept = 5.5, colour = "grey30", linewidth = 0.4) +
  geom_hline(yintercept = 5.5, colour = "grey30", linewidth = 0.4) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), direction = 1, name = "Spearman\ncorrelation") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_light(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(here("output", "figs", "supp", "figS2_correlation.png"), fig_s2, width = 8, height = 6.8, dpi = 600, bg = "white")

# 04 Figure S3: connectivity metrics on the grid ----
# the grid cells holding survey stations only; each metric keeps its own
# single-hue ramp (light = low, dark = high); skewed ones use a square-root
# scale, which keeps the zeros that a log scale would drop
# transitivity is undefined in cells receiving no larvae (in-degree 0, stored
# as 0); shown in grey so they don't stretch the scale
grid_s3 <- readRDS(here("data", "grid", "grid_env_conn.rds")) |>
  semi_join(sampled_cells, by = c("x", "y")) |>
  mutate(conn_transitivity_presence = if_else(conn_in_degree_presence == 0, NA, conn_transitivity_presence))

metric_titles <- c(
  conn_in_degree_presence = "In-degree",
  conn_in_strength_presence = "In-strength",
  conn_in_closeness_presence = "In-closeness",
  conn_eigen_presence = "Eigenvector centrality",
  conn_transitivity_presence = "Transitivity"
)
metric_scales <- list(
  conn_in_degree_presence = scale_fill_gradient(low = "#d0e1f2", high = "#08306b", name = NULL),
  conn_in_strength_presence = scale_fill_gradient(low = "#efedf5", high = "#3f007d", name = NULL, transform = "sqrt"),
  conn_in_closeness_presence = scale_fill_gradient(low = "#dbf1f0", high = "#01665e", name = NULL),
  conn_eigen_presence = scale_fill_gradient(low = "#fee6ce", high = "#7f2704", name = NULL, transform = "sqrt"),
  conn_transitivity_presence = scale_fill_gradient(low = "#e5f5e0", high = "#005a32", name = NULL, na.value = "grey80")
)

metric_panels <- map(names(metric_titles), function(col) {
  ggplot() +
    geom_tile(data = grid_s3, aes(x, y, fill = .data[[col]]), width = 2000, height = 2000) +
    geom_sf(data = land_mask, fill = "white", colour = NA) +
    geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.22) +
    metric_scales[[col]] +
    coord_sf(xlim = xlim_map, ylim = ylim_map, crs = 32632, expand = FALSE) +
    labs(title = metric_titles[[col]]) +
    map_theme
})

fig_s3 <- wrap_plots(metric_panels, ncol = 3) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  tag_theme
ggsave(here("output", "figs", "supp", "figS3_connectivity_metrics.png"), fig_s3, width = 13, height = 7, dpi = 600, bg = "white")

# 05 Figure S4: barrier mesh ----
# the mesh (03_build_mesh.R) as used in the fitted models: water triangles
# carry the spatial field, land triangles are the barrier
barrier_mesh <- fit_full$spde
triangles <- barrier_mesh$mesh$graph$tv
mesh_sf <- map(seq_len(nrow(triangles)), \(i) st_polygon(list(barrier_mesh$mesh$loc[triangles[i, c(1:3, 1)], 1:2] * 1000))) |>
  st_sfc(crs = 32632) |>
  st_sf(geometry = _) |>
  mutate(barrier = seq_len(n()) %in% barrier_mesh$barrier_triangles)

fig_s4 <- ggplot() +
  geom_sf(data = mesh_sf, aes(fill = barrier), colour = "grey45", linewidth = 0.1) +
  geom_sf(data = water, fill = NA, colour = "black", linewidth = 0.35) +
  geom_point(data = dat_model, aes(x_utm * 1000, y_utm * 1000, colour = "Survey station"), size = 0.12) +
  scale_fill_manual(values = c(`FALSE` = "#d6e6f2", `TRUE` = "grey85"), labels = c("Water triangle", "Land triangle (barrier)"), name = NULL) +
  scale_colour_manual(values = c(`Survey station` = "#B2182B"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  coord_sf(crs = 32632, xlim = range(dat_model$x_utm * 1000) + c(-12000, 12000), ylim = range(dat_model$y_utm * 1000) + c(-8000, 8000), expand = FALSE) +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "bottom",
    plot.margin = margin(8, 8, 8, 8)
  )
ggsave(here("output", "figs", "supp", "figS4_mesh.png"), fig_s4, width = 10, height = 8, dpi = 600, bg = "white")

# 06 Figure S5: cross-validation folds ----
# the 2 x 2 km grid cells used as cross-validation blocks, coloured by fold, as
# saved by 08_cross_validation.R: (a) occurrence, all cells; (b) positive
# biomass, its own folds over the cells holding stations with cockles
cv_results <- readRDS(here("data", "sdm", "cv", "cv_results.rds"))
fold_cells <- function(cv_data) {
  cv_data |>
    transmute(x = cell_centre_x(x_utm), y = cell_centre_y(y_utm), fold = factor(fold, levels = 1:10)) |>
    distinct()
}
fold_panel <- function(cells, title) {
  ggplot() +
    geom_tile(data = cells, aes(x, y, fill = fold), width = 2000, height = 2000) +
    geom_sf(data = water, fill = NA, colour = "grey40", linewidth = 0.25) +
    scale_fill_brewer(palette = "Paired", name = "Fold", drop = FALSE) +
    coord_sf(crs = 32632, xlim = xlim_map, ylim = ylim_map, expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = 11, hjust = 0.5, colour = "grey20"),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 1))
}

fig_s5 <- (fold_panel(fold_cells(cv_results$cv_presence[[1]]$data), "Occurrence") +
  fold_panel(fold_cells(cv_results$cv_biomass[[1]]$data), "Biomass where present")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(legend.position = "bottom") &
  tag_theme
ggsave(here("output", "figs", "supp", "figS5_cv_folds.png"), fig_s5, width = 13, height = 6.5, dpi = 600, bg = "white")

# 07 Figure S6: residuals of the full model ----
# simulation-based (DHARMa) residuals of the delta-gamma model as a whole:
# 500 biomass draws per station from the fitted model with parameter
# uncertainty (type = "mle-mvn"), and each observation's quantile among its
# draws. integerResponse = TRUE randomises the quantile within ties, needed
# for the point mass at zero (it changes nothing for positive biomass, which
# has no ties). Uniform on 0-1 if the model is right: Q-Q plot against the
# uniform, and each station's residual on the map
set.seed(1)
sim_s6 <- simulate(fit_full, nsim = 500, type = "mle-mvn")
cat(
  "S6 zeros: observed", round(mean(dat_model$biomass == 0), 3),
  "simulated", round(mean(sim_s6 == 0), 3), "\n"
)
dharma_s6 <- dharma_residuals(sim_s6, fit_full, plot = FALSE, return_DHARMa = TRUE, integerResponse = TRUE)
resid_s6 <- dat_model |>
  select(x_utm, y_utm) |>
  mutate(resid = dharma_s6$scaledResiduals)

qq_s6 <- tibble(observed = sort(resid_s6$resid), expected = seq_along(observed) / (length(observed) + 1)) |>
  ggplot(aes(expected, observed)) +
  geom_point(size = 0.4, alpha = 0.4, colour = "grey25") +
  geom_abline(colour = "#B2182B") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Expected (uniform)", y = "Observed residual") +
  theme_light(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# largest departures drawn last, on top
map_s6 <- ggplot() +
  geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.22) +
  geom_point(
    data = arrange(resid_s6, abs(resid - 0.5)),
    aes(x_utm * 1000, y_utm * 1000, colour = resid), size = 0.2, alpha = 0.5
  ) +
  scale_colour_distiller(palette = "RdBu", limits = c(0, 1), direction = 1, name = "Residual") +
  coord_sf(xlim = xlim_map, ylim = ylim_map, crs = 32632, expand = FALSE) +
  map_theme +
  theme(legend.title = element_text(size = 9, colour = "grey20"))

fig_s6 <- qq_s6 + map_s6 +
  plot_layout(widths = c(1, 1.6)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  tag_theme
ggsave(here("output", "figs", "supp", "figS6_residuals.png"), fig_s6, width = 11, height = 5, dpi = 600, bg = "white")

# 08 Figure S8: connectivity coefficient of each metric ----
# the full model (space + environment + connectivity) refitted with each of
# the five connectivity metrics in turn, per component
coefs_s8 <- expand_grid(metric = names(conn_labs), component = 1:2) |>
  mutate(coef = map2(metric, component, \(metric, component) {
    fit <- if (metric == "conn_in_strength") {
      fit_full
    } else {
      read_fit("sensitivity", paste0("space_env_conn_", str_remove(metric, "^conn_")))
    }
    tidy(fit, effects = "fixed", model = component, conf.int = TRUE) |>
      filter(term == paste0(metric, "_std"))
  })) |>
  unnest(coef) |>
  mutate(
    metric = factor(conn_labs[metric], levels = rev(unname(conn_labs))),
    part = factor(unname(part_labs)[component], levels = unname(part_labs))
  )

fig_s8 <- ggplot(coefs_s8, aes(estimate, metric)) +
  geom_stripped_rows(colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0, linewidth = 0.9, colour = "#E69F00") +
  geom_point(size = 2.4, colour = "#E69F00") +
  facet_wrap(~part) +
  labs(x = "Standardized connectivity coefficient", y = NULL) +
  theme_light(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
  facet_theme
ggsave(here("output", "figs", "supp", "figS8_connectivity_coefficients.png"), fig_s8, width = 9, height = 4, dpi = 600, bg = "white")
