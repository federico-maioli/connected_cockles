# All supplementary figures except S6 (predictions, figS6_predictions.R), one
# section each, numbered in the order they are cited in the manuscript:
#   S1   survey stations by year and survey
#   S2   the barrier mesh
#   S3   cross-validation folds (2 km grid cells)
#   S4   model residuals (not yet made)
#   S5   connectivity coefficient of each metric in the full model
#   S7   connectivity vs environment and the residual spatial field
#   S8   connectivity metrics on the 2 km grid
#   S9   correlation among all model predictors
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

# map theme for S1
map_theme <- theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 9.5, hjust = 0.5, colour = "grey20"),
    legend.position = "right",
    legend.key.width = unit(3, "mm"), legend.key.height = unit(8, "mm"),
    legend.text = element_text(size = 7, colour = "grey30"),
    plot.tag = element_text(size = 12, face = "bold", colour = "grey20"),
    plot.margin = margin(4, 4, 4, 4)
  )

# 01 Shared data ----
water <- st_read(here("data", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

# the model dataset (complete cases, standardised covariates, presence-weighted
# connectivity), taken from a saved fit so it matches the models exactly
dat_model <- read_fit("main", "space_env")$data

dir.create(here("output", "figs", "supp"), showWarnings = FALSE, recursive = TRUE)

# the 2 km grid with connectivity (S1): only wet cells
# inside Limfjorden - the grid extends past the fjord into the North Sea,
# where cells have no larval links
grid <- readRDS(here("data", "grid", "grid_env_conn.rds")) |>
  filter(!is.na(depth))
inside <- lengths(st_intersects(st_as_sf(grid, coords = c("x", "y"), crs = 32632), water)) > 0
grid <- grid[inside, ]

xlim_grid <- range(grid$x) + c(-4000, 4000)
ylim_grid <- range(grid$y) + c(-4000, 4000)

# white mask = the plot frame minus the fjord, drawn over the tiles so the 2 km
# cells are clipped to the coastline instead of spilling onto land
frame <- st_as_sfc(st_bbox(
  c(xmin = xlim_grid[1], ymin = ylim_grid[1], xmax = xlim_grid[2], ymax = ylim_grid[2]),
  crs = st_crs(32632)
))
land_mask <- st_difference(frame, st_union(water))

# 02 Figure S1: connectivity metrics on the grid ----
# each metric keeps its own single-hue ramp (light = low, dark = high); skewed
# ones use a square-root scale, which keeps the zeros that a log scale would drop
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
  conn_transitivity_presence = scale_fill_gradient(low = "#e5f5e0", high = "#005a32", name = NULL)
)

metric_panels <- map(names(metric_titles), function(col) {
  ggplot() +
    geom_tile(data = grid, aes(x, y, fill = .data[[col]]), width = 2000, height = 2000) +
    geom_sf(data = land_mask, fill = "white", colour = NA) +
    geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.22) +
    metric_scales[[col]] +
    coord_sf(xlim = xlim_grid, ylim = ylim_grid, crs = 32632, expand = FALSE) +
    labs(title = metric_titles[[col]]) +
    map_theme
})

fig_s1 <- wrap_plots(metric_panels, ncol = 3) + plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave(here("output", "figs", "supp", "figS8_connectivity_metrics.png"), fig_s1, width = 13, height = 5.6, dpi = 600, bg = "white")

# 03 Figure S2: correlation among all model predictors ----
# environment and connectivity together, so cross-block collinearity is visible
s2_vars <- c(
  Depth = "depth", Temperature = "temp", Salinity = "sal", Oxygen = "oxy", `Shear stress` = "shear_max",
  `In-degree` = "conn_in_degree", `In-strength` = "conn_in_strength", `In-closeness` = "conn_in_closeness",
  Eigenvector = "conn_eigen", Transitivity = "conn_transitivity"
)

cor_s2 <- dat_model |>
  select(all_of(s2_vars)) |>
  cor(use = "complete.obs") |>
  as_tibble(rownames = "var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "r") |>
  mutate(
    var1 = factor(var1, levels = names(s2_vars)),
    var2 = factor(var2, levels = rev(names(s2_vars)))
  )

fig_s2 <- ggplot(cor_s2, aes(var1, var2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.8) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), direction = 1, name = "Pearson r") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_light(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(here("output", "figs", "supp", "figS9_correlation.png"), fig_s2, width = 7, height = 6, dpi = 600, bg = "white")

# 04 Figure S3: connectivity coefficient of each metric ----
# the full model (space + environment + connectivity) refitted with each of
# the five connectivity metrics in turn, per component
coefs_s3 <- expand_grid(metric = names(conn_labs), component = 1:2) |>
  mutate(coef = map2(metric, component, \(metric, component) {
    fit <- if (metric == "conn_in_strength") {
      read_fit("main", "space_env_conn")
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

fig_s3 <- ggplot(coefs_s3, aes(estimate, metric)) +
  geom_stripped_rows(colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0, linewidth = 0.9, colour = "#E69F00") +
  geom_point(size = 2.4, colour = "#E69F00") +
  facet_wrap(~part) +
  labs(x = "Standardized connectivity coefficient", y = NULL) +
  theme_light(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
  facet_theme
ggsave(here("output", "figs", "supp", "figS5_connectivity_coefficients.png"), fig_s3, width = 9, height = 4, dpi = 600, bg = "white")

# 05 Figure S4: connectivity vs environment and the spatial field ----
# are the connectivity metrics collinear with the environment, or standing in
# for residual spatial structure? The last two columns are the spatial random
# fields of the Space + Environment model, one per component - not
# environmental covariates, hence the separator. Computed at the survey stations (the model rows), so it
# is weighted towards heavily sampled cells. Spearman: the metrics are skewed
s4_conn <- conn_labs
s4_env <- c(
  depth = "Depth", temp = "Temperature", oxy = "Oxygen", sal = "Salinity", shear_max = "Shear stress",
  field_presence = "Spatial field (presence)", field_biomass = "Spatial field (biomass)"
)

pred_s4 <- predict(read_fit("main", "space_env"))
dat_s4 <- dat_model
dat_s4$field_presence <- pred_s4$est_rf1
dat_s4$field_biomass <- pred_s4$est_rf2

cors_s4 <- cor(dat_s4[, names(s4_conn)], dat_s4[, names(s4_env)], method = "spearman", use = "pairwise.complete.obs") |>
  as.data.frame() |>
  rownames_to_column("metric") |>
  pivot_longer(-metric, names_to = "covariate", values_to = "rho") |>
  mutate(
    metric = factor(s4_conn[metric], levels = rev(unname(s4_conn))),
    covariate = factor(s4_env[covariate], levels = unname(s4_env))
  )

fig_s4 <- ggplot(cors_s4, aes(covariate, metric, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.1) +
  geom_vline(xintercept = 5.5, colour = "grey30", linewidth = 0.5) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), direction = 1, name = "Spearman correlation", guide = guide_colourbar(title.position = "top")) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_light(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    legend.key.height = unit(3, "mm"), legend.key.width = unit(14, "mm")
  )
ggsave(here("output", "figs", "supp", "figS7_connectivity_environment.png"), fig_s4, width = 6.4, height = 5.6, dpi = 600, bg = "white")


# 06 Figure S5: barrier mesh ----
# the mesh (03_build_mesh.R) as used in the fitted models: water triangles
# carry the spatial field, land triangles are the barrier
barrier_mesh <- read_fit("main", "space_env_conn")$spde
triangles <- barrier_mesh$mesh$graph$tv
mesh_sf <- map(seq_len(nrow(triangles)), \(i) st_polygon(list(barrier_mesh$mesh$loc[triangles[i, c(1:3, 1)], 1:2] * 1000))) |>
  st_sfc(crs = 32632) |>
  st_sf(geometry = _) |>
  mutate(barrier = seq_len(n()) %in% barrier_mesh$barrier_triangles)

fig_s5 <- ggplot() +
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
ggsave(here("output", "figs", "supp", "figS2_mesh.png"), fig_s5, width = 10, height = 8, dpi = 600, bg = "white")

# 07 Figure S6: survey stations by year ----
# every station used in the models, one panel per year, coloured by survey
survey_labs <- c(Stock = "Grab (2021-2025)", Stock2018 = "Suction dredge (2018)", KSKV = "KSKV dredge")
stations_s6 <- dat_model |>
  mutate(survey = factor(survey_labs[as.character(survey)], levels = survey_labs))

fig_s6 <- ggplot() +
  geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.25) +
  geom_point(data = stations_s6, aes(x_utm * 1000, y_utm * 1000, colour = survey), size = 0.5, alpha = 0.6) +
  facet_wrap(~year, ncol = 4) +
  scale_colour_manual(values = c("#1F4E79", "#E69F00", "#CC79A7"), name = NULL) +
  coord_sf(crs = 32632, xlim = range(dat_model$x_utm * 1000) + c(-6000, 6000), ylim = range(dat_model$y_utm * 1000) + c(-6000, 6000), expand = FALSE) +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(colour = "grey15", margin = margin(3, 0, 3, 0)),
    legend.position = "bottom"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
ggsave(here("output", "figs", "supp", "figS1_survey_stations.png"), fig_s6, width = 11, height = 6.5, dpi = 600, bg = "white")

# 08 Figure S7: cross-validation folds ----
# the 2 x 2 km grid cells used as cross-validation blocks, coloured by the fold
# they were assigned to, as saved by 08_cross_validation.R
folds_s7 <- readRDS(here("data", "sdm", "cv", "cv_results.rds"))$cv_presence[[1]]$data |>
  mutate(
    x = 450074 + floor((x_utm * 1000 - 450074) / 2000) * 2000 + 1000,
    y = 6258093 + floor((y_utm * 1000 - 6258093) / 2000) * 2000 + 1000
  ) |>
  distinct(x, y, fold)

fig_s7 <- ggplot() +
  geom_tile(data = folds_s7, aes(x, y, fill = factor(fold)), width = 2000, height = 2000) +
  geom_sf(data = water, fill = NA, colour = "grey40", linewidth = 0.25) +
  scale_fill_brewer(palette = "Paired", name = "Fold") +
  coord_sf(crs = 32632, xlim = range(folds_s7$x) + c(-6000, 6000), ylim = range(folds_s7$y) + c(-6000, 6000), expand = FALSE) +
  theme_void(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))
ggsave(here("output", "figs", "supp", "figS3_cv_folds.png"), fig_s7, width = 9, height = 7, dpi = 600, bg = "white")
