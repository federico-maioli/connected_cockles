# All supplementary figures except S4 (predictions, figS4_predictions.R) and S8 (accessibility, figS8_accessibility.R), one
# section each, numbered in the order they are cited in main.tex:
#   S1   survey stations by year and survey, presence and absence
#   S2   connectivity metrics on the 2 km grid cells with stations
#   S3   correlation among all covariates (environment and connectivity)
#   S5   the barrier mesh
#   S6   cross-validation folds, occurrence and positive biomass
#   S7   simulation-based residuals of the full model
#   S9   connectivity coefficient of each metric in the full model
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

# the two parts of each model, in order (component 1, component 2)
part_labs <- c(presence = "Occurrence (logit link)", biomass = "Biomass where present (log link)")
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
# (the occurrence part holds every station; list(occurrence, biomass))
fit_full <- read_fit("main", "space_env_conn")
dat_model <- fit_full$occurrence$data

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
# every station used in the models, one panel per year, coloured by survey
# (same colours as Figure 1); absences as small, semi-transparent open
# circles, presences (biomass > 0) as larger filled, white-outlined circles
# drawn on top. Full-page figure
# (7 x 7.8 in) at 600 dpi with an 8 pt base font
survey_labs <- c(Stock = "Grab (2021-2025)", Stock2018 = "Suction dredge (2018)", KSKV = "KSKV dredge (2018-2025)")
survey_cols <- set_names(c("#1F4E79", "#E69F00", "#CC79A7"), survey_labs)
stations_s1 <- dat_model |>
  mutate(
    survey = factor(survey_labs[as.character(survey)], levels = survey_labs),
    cockles = factor(if_else(present == 1, "Present", "Absent"), levels = c("Absent", "Present"))
  ) |>
  arrange(cockles)

fig_s1 <- ggplot() +
  geom_sf(data = water, fill = "grey97", colour = "grey55", linewidth = 0.2) +
  # absences: small, semi-transparent open circles, so densely sampled grab
  # beds read as a light area
  geom_point(
    data = filter(stations_s1, cockles == "Absent"),
    aes(x_utm * 1000, y_utm * 1000, colour = survey, shape = cockles),
    size = 0.4, stroke = 0.25, alpha = 0.3
  ) +
  # presences: filled circles with a thin white outline, so they stand out
  # even on top of dense absences of the same colour
  geom_point(
    data = filter(stations_s1, cockles == "Present"),
    aes(x_utm * 1000, y_utm * 1000, fill = survey, shape = cockles),
    size = 1.1, stroke = 0.15, colour = "white"
  ) +
  facet_wrap(~year, ncol = 3) +
  scale_colour_manual(values = survey_cols, name = NULL, aesthetics = c("colour", "fill")) +
  scale_shape_manual(values = c(Absent = 1, Present = 21), name = NULL) +
  coord_sf(crs = 32632, xlim = xlim_map, ylim = ylim_map, expand = FALSE) +
  theme_void(base_size = 8) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(size = 9, face = "bold", colour = "grey15", margin = margin(3, 0, 2, 0)),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "mm"),
    legend.text = element_text(size = 8),
    panel.spacing = unit(3, "mm"),
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(
    # one survey legend (colour and fill share the scale): filled keys
    fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 2.8, colour = "white", alpha = 1, stroke = 0.2)),
    colour = "none",
    shape = guide_legend(order = 2, override.aes = list(size = 2.5, colour = "grey30", fill = "grey30", stroke = 0.6, alpha = 1))
  )
ggsave(here("output", "figs", "supp", "figS1_survey_stations.png"), fig_s1, width = 7, height = 7.8, dpi = 600, bg = "white")

# 03 Figure S3: correlation among all covariates ----
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
ggsave(here("output", "figs", "supp", "figS3_correlation.png"), fig_s2, width = 8, height = 6.8, dpi = 600, bg = "white")

# 04 Figure S2: connectivity metrics on the grid ----
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
ggsave(here("output", "figs", "supp", "figS2_connectivity_metrics.png"), fig_s3, width = 13, height = 7, dpi = 600, bg = "white")

# 05 Figure S5: barrier mesh ----
# the mesh (03_build_mesh.R) as used in the fitted models: water triangles
# carry the spatial field, land triangles are the barrier
barrier_mesh <- fit_full$occurrence$spde
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
ggsave(here("output", "figs", "supp", "figS5_mesh.png"), fig_s4, width = 10, height = 8, dpi = 600, bg = "white")

# 06 Figure S6: cross-validation folds ----
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
  fold_panel(fold_cells(cv_results$cv_biomass[[1]]$data), "Biomass where present (grab)")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(legend.position = "bottom") &
  tag_theme
ggsave(here("output", "figs", "supp", "figS6_cv_folds.png"), fig_s5, width = 13, height = 6.5, dpi = 600, bg = "white")

# 07 Figure S7: residuals of the full model ----
# simulation-based (DHARMa) residuals of each part: 500 draws per station from
# the fitted part with parameter uncertainty (type = "mle-mvn"), and each
# observation's quantile among its draws. integerResponse = TRUE randomises
# the quantile within ties, needed for the 0/1 occurrence data (it changes
# nothing for positive biomass, which has no ties). Uniform on 0-1 if the
# model is right: Q-Q plot against the uniform, and each station's residual
# on the map. Occurrence: every station; biomass where present: grab stations
# with cockles
set.seed(1)
resid_s6 <- imap(c(presence = "occurrence", biomass = "biomass"), \(part, nm) {
  fit <- fit_full[[part]]
  sim <- simulate(fit, nsim = 500, type = "mle-mvn")
  res <- dharma_residuals(sim, fit, plot = FALSE, return_DHARMa = TRUE, integerResponse = TRUE)
  fit$data |>
    select(x_utm, y_utm) |>
    mutate(resid = res$scaledResiduals, part = part_labs[[nm]])
}) |>
  list_rbind() |>
  mutate(part = factor(part, levels = unname(part_labs)))

# shared strip style for both columns, with room above and below the title
# so it is not clipped; 7 in wide (double column) at 600 dpi, 8 pt text
strip_s6 <- theme(
  strip.background = element_blank(),
  strip.text = element_text(colour = "black", size = 8, margin = margin(t = 2, b = 4)),
  strip.clip = "off"
)

qq_s6 <- resid_s6 |>
  group_by(part) |>
  arrange(resid, .by_group = TRUE) |>
  mutate(expected = row_number() / (n() + 1)) |>
  ungroup() |>
  ggplot(aes(expected, resid)) +
  geom_point(size = 0.3, alpha = 0.4, colour = "grey25") +
  geom_abline(colour = "#B2182B", linewidth = 0.4) +
  facet_wrap(~part, ncol = 1) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Expected (uniform)", y = "Observed residual") +
  theme_light(base_size = 8) +
  theme(panel.grid.minor = element_blank()) +
  strip_s6

# largest departures drawn last, on top; extent: the sampled area with a 3 km
# margin, as in Figure 2
map_s6 <- ggplot() +
  geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.15) +
  geom_point(
    data = arrange(resid_s6, abs(resid - 0.5)),
    aes(x_utm * 1000, y_utm * 1000, colour = resid), size = 0.15, alpha = 0.5
  ) +
  facet_wrap(~part, ncol = 1) +
  scale_colour_distiller(palette = "RdBu", limits = c(0, 1), direction = 1, name = "Residual", breaks = c(0, 0.5, 1)) +
  coord_sf(
    xlim = range(dat_model$x_utm * 1000) + c(-3000, 3000), ylim = range(dat_model$y_utm * 1000) + c(-3000, 3000),
    crs = 32632, expand = FALSE
  ) +
  guides(colour = guide_colourbar(title.position = "left", title.vjust = 1, theme = theme(
    legend.key.width = unit(16, "mm"), legend.key.height = unit(1.4, "mm")
  ))) +
  theme_void(base_size = 8) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "bottom",
    legend.title = element_text(size = 7, colour = "grey20"),
    legend.text = element_text(size = 6.5, colour = "grey30")
  ) +
  strip_s6

fig_s6 <- qq_s6 + map_s6 +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 9, face = "bold", colour = "grey20"))
ggsave(here("output", "figs", "supp", "figS7_residuals.png"), fig_s6, width = 7, height = 6.4, dpi = 600, bg = "white")

# 08 Figure S9: connectivity coefficient of each metric ----
# the full model (space + environment + connectivity) refitted with each of
# the five connectivity metrics in turn, per component
coefs_s8 <- expand_grid(metric = names(conn_labs), component = 1:2) |>
  mutate(coef = map2(metric, component, \(metric, component) {
    fit <- if (metric == "conn_in_strength") {
      fit_full
    } else {
      read_fit("sensitivity", paste0("space_env_conn_", str_remove(metric, "^conn_")))
    }
    tidy(fit[[component]], effects = "fixed", conf.int = TRUE) |>
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
ggsave(here("output", "figs", "supp", "figS9_connectivity_coefficients.png"), fig_s8, width = 9, height = 4, dpi = 600, bg = "white")
