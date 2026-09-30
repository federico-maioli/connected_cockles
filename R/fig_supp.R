# All supplementary figures, one section each:
#   S1   connectivity metrics on the 2 km grid
#   S3   correlation among all model predictors
#   S4   connectivity coefficient with vs without the spatial field
#   S5   variance partitioning across the 5 connectivity metrics
#   S6   stability of the in-strength coefficient across specifications
#   S7   connectivity vs environment and the residual spatial field
# (S2 is the presence version of fig3, made in fig3_coeff.R; the covariate
# maps are in the main text, fig2_predictors.R.)
#
# Connectivity metrics are the presence-weighted ones used in 06_fit_sdm.R;
# every model shown is read from its saved fit (data/sdm/main and
# data/sdm/sensitivity), nothing is refitted here.

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

# fitted model saved by 06_fit_sdm.R: in-strength models live in main/, the
# other metrics and the leave-one-out specs in sensitivity/
read_fit <- function(folder, resp, model) {
  readRDS(here("data", "sdm", folder, paste0(resp, "_", model, ".rds")))
}

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
dat_model <- read_fit("main", "biomass", "space_env")$data

comparison <- bind_rows(
  readRDS(here("data", "sdm", "main", "model_comparison.rds")),
  readRDS(here("data", "sdm", "sensitivity", "model_comparison.rds"))
)

dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)

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
  conn_transitivity_presence = "Transitivity",
  conn_in_strength_flipped = "In-strength, flipped"
)
metric_scales <- list(
  conn_in_degree_presence = scale_fill_gradient(low = "#d0e1f2", high = "#08306b", name = NULL),
  conn_in_strength_presence = scale_fill_gradient(low = "#efedf5", high = "#3f007d", name = NULL, transform = "sqrt"),
  conn_in_closeness_presence = scale_fill_gradient(low = "#dbf1f0", high = "#01665e", name = NULL),
  conn_eigen_presence = scale_fill_gradient(low = "#fee6ce", high = "#7f2704", name = NULL, transform = "sqrt"),
  conn_transitivity_presence = scale_fill_gradient(low = "#e5f5e0", high = "#005a32", name = NULL),
  conn_in_strength_flipped = scale_fill_gradient(low = "#fde0dd", high = "#ae017e", name = NULL, transform = "sqrt")
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
ggsave(here("output", "figs", "figS1_connectivity_metrics.png"), fig_s1, width = 13, height = 5.6, dpi = 600, bg = "white")

# 03 Figure S3: correlation among all model predictors ----
# environment and connectivity together, so cross-block collinearity is visible
s3_vars <- c(
  Depth = "depth", Temperature = "temp", Salinity = "sal", Oxygen = "oxy", `Shear stress` = "shear_max",
  `In-degree` = "conn_in_degree", `In-strength` = "conn_in_strength", `In-closeness` = "conn_in_closeness",
  Eigenvector = "conn_eigen", Transitivity = "conn_transitivity"
)

cor_s3 <- dat_model |>
  select(all_of(s3_vars)) |>
  cor(use = "complete.obs") |>
  as_tibble(rownames = "var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "r") |>
  mutate(
    var1 = factor(var1, levels = names(s3_vars)),
    var2 = factor(var2, levels = rev(names(s3_vars)))
  )

fig_s3 <- ggplot(cor_s3, aes(var1, var2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.8) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), direction = 1, name = "Pearson r") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_light(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(here("output", "figs", "figS3_correlation.png"), fig_s3, width = 7, height = 6, dpi = 600, bg = "white")

# 04 Figure S4: connectivity coefficient with vs without the spatial field ----
# the same Environment + Connectivity model with and without the field
# (space_env_conn vs env_conn). If connectivity were standing in for broad
# spatial structure, dropping the field would inflate its coefficient
structure_labs <- c(space_env_conn = "With spatial field", env_conn = "Without spatial field")

s4_models <- expand_grid(
  response = c("biomass", "present"),
  structure = names(structure_labs),
  metric = names(conn_labs)
) |>
  mutate(
    folder = if_else(metric == "conn_in_strength", "main", "sensitivity"),
    model = if_else(metric == "conn_in_strength", structure, paste0(structure, "_", str_remove(metric, "^conn_")))
  ) |>
  semi_join(filter(comparison, converged), by = c("response", "model"))

coefs_s4 <- s4_models |>
  mutate(coef = pmap(list(folder, response, model, metric), \(folder, response, model, metric) {
    tidy(read_fit(folder, response, model), effects = "fixed", conf.int = TRUE) |>
      filter(term == paste0(metric, "_std"))
  })) |>
  unnest(coef) |>
  mutate(
    metric = factor(conn_labs[metric], levels = rev(unname(conn_labs))),
    structure = factor(structure_labs[structure], levels = unname(structure_labs)),
    response = factor(response, levels = c("present", "biomass"), labels = c("Presence (logit)", "Biomass (log)"))
  )

fig_s4 <- ggplot(coefs_s4, aes(estimate, metric, colour = structure)) +
  geom_stripped_rows(colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0, position = position_dodge(width = 0.55), linewidth = 0.9) +
  geom_point(position = position_dodge(width = 0.55), size = 2.4) +
  facet_wrap(~response, scales = "free_x") +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
  labs(x = "Standardized connectivity coefficient", y = NULL) +
  theme_light(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "grey20", face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(here("output", "figs", "figS4_space_confounding.png"), fig_s4, width = 9, height = 4.5, dpi = 600, bg = "white")

# how much each connectivity coefficient moves when the field is dropped
coefs_s4 |>
  select(response, metric, structure, estimate) |>
  pivot_wider(names_from = structure, values_from = estimate) |>
  rename(with_field = `With spatial field`, without_field = `Without spatial field`) |>
  mutate(across(c(with_field, without_field), \(x) round(x, 3)), ratio = round(without_field / with_field, 2)) |>
  arrange(response, metric) |>
  print(n = Inf)

# 05 Figure S5: variance partitioning across connectivity metrics ----
# spatial field on, both responses. Unexplained (distribution) variance is not
# drawn, as in fig2 - it is the empty remainder of each bar. Segments >= 15%
# are labelled in place; thinner ones are listed to the right of the bar
part_cols <- c("Spatial field" = "#0072B2", "Environment" = "#009E73", "Connectivity" = "#E69F00")
label_cols <- c("Spatial field" = "white", "Environment" = "white", "Connectivity" = "grey15")

s5 <- expand_grid(response = c("biomass", "present"), metric = names(conn_labs)) |>
  mutate(variance = map2(response, metric, \(resp, m) {
    if (m == "conn_in_strength") {
      readRDS(here("data", "sdm", "main", paste0(resp, "_space_env_conn_variance.rds")))
    } else {
      readRDS(here("data", "sdm", "sensitivity", paste0(resp, "_space_env_conn_", str_remove(m, "^conn_"), "_variance.rds")))
    }
  })) |>
  unnest(variance) |>
  filter(share > 0, component != "distribution") |>
  mutate(
    component = factor(recode(component, spatial = "Spatial field"), levels = names(part_cols)),
    group = factor(conn_labs[metric], levels = rev(unname(conn_labs))),
    response_lab = recode(response, biomass = "Biomass (log link)", present = "Presence (logit link)")
  )

s5_labs <- s5 |>
  arrange(response_lab, group, component) |>
  mutate(x_mid = cumsum(share) - share / 2, .by = c(response_lab, group))
s5_narrow <- s5_labs |>
  filter(share < 0.15) |>
  summarise(label = paste(component, scales::percent(share, accuracy = 0.1), collapse = "\n"), .by = c(response_lab, group))

fig_s5 <- ggplot(s5, aes(x = share, y = group, fill = component)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
  geom_text(
    data = filter(s5_labs, share >= 0.15),
    aes(x = x_mid, y = group, label = scales::percent(share, accuracy = 0.1), colour = component),
    inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  geom_text(
    data = s5_narrow, aes(x = 1.02, y = group, label = label),
    inherit.aes = FALSE, size = 2.6, colour = "grey25", hjust = 0, vjust = 0.5, lineheight = 0.85
  ) +
  scale_fill_manual(values = part_cols, name = NULL) +
  scale_colour_manual(values = label_cols, guide = "none") +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = scales::percent, limits = c(0, 1.4), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  facet_wrap(~response_lab) +
  labs(x = expression("Share of total variance (Nakagawa " * R^2 * ")"), y = NULL) +
  theme_light(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(colour = "black")
  ) +
  guides(fill = guide_legend(nrow = 1))
ggsave(here("output", "figs", "figS5_variance_metrics.png"), fig_s5, width = 9.5, height = 4.4, dpi = 600, bg = "white")

# 06 Figure S6: stability of the in-strength coefficient ----
# rows: full model, full minus one environment term (x5), connectivity only;
# columns: spatial field on / off. A dagger marks a 95% CI crossing zero.
# Connectivity only is fitted in 06 with the spatial field only, so its
# field-off cell is empty
env_labs <- c(depth = "Depth", temp = "Temperature", sal = "Salinity", oxy = "Oxygen", shear_max = "Shear stress")

s6_specs <- bind_rows(
  tibble(label = "Full model", spatial = c("on", "off"), folder = "main", model = c("space_env_conn", "env_conn")),
  map(names(env_labs), \(v) tibble(
    label = paste("–", env_labs[[v]]), spatial = c("on", "off"), folder = "sensitivity",
    model = paste0(c("space_env_conn_drop_", "env_conn_drop_"), v)
  )) |> list_rbind(),
  tibble(label = "Connectivity only (no environment)", spatial = "on", folder = "main", model = "space_conn")
)

coefs_s6 <- expand_grid(response = c("biomass", "present"), s6_specs) |>
  semi_join(filter(comparison, converged), by = c("response", "model")) |>
  mutate(coef = pmap(list(folder, response, model), \(folder, response, model) {
    tidy(read_fit(folder, response, model), effects = "fixed", conf.int = TRUE) |>
      filter(term == "conn_in_strength_std")
  })) |>
  unnest(coef) |>
  mutate(
    label = factor(label, levels = rev(unique(s6_specs$label))),
    spatial_lab = factor(paste("Spatial field:", spatial), levels = c("Spatial field: on", "Spatial field: off")),
    response_lab = recode(response, biomass = "Biomass (log link)", present = "Presence (logit link)"),
    cell_label = paste0(sprintf("%.2f", estimate), if_else(conf.low <= 0 & conf.high >= 0, "†", ""))
  )

fill_limit <- max(abs(coefs_s6$estimate))

fig_s6 <- ggplot(coefs_s6, aes(spatial_lab, label, fill = estimate)) +
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
ggsave(here("output", "figs", "figS6_coeff_stability.png"), fig_s6, width = 9, height = 4.4, dpi = 600, bg = "white")

# 07 Figure S7: connectivity vs environment and the spatial field ----
# are the connectivity metrics collinear with the environment, or standing in
# for residual spatial structure? The last column is the spatial random field
# of the Space + Environment biomass model - not an environmental covariate,
# hence the separator. Computed at the survey stations (the model rows), so it
# is weighted towards heavily sampled cells. Spearman: the metrics are skewed
s7_conn <- c(conn_labs, conn_in_strength_flipped = "In-strength, flipped")
s7_env <- c(depth = "Depth", temp = "Temperature", oxy = "Oxygen", sal = "Salinity", shear_max = "Shear stress", spatial_field = "Spatial field")

dat_s7 <- dat_model
dat_s7$spatial_field <- predict(read_fit("main", "biomass", "space_env"))$est_rf

cors_s7 <- cor(dat_s7[, names(s7_conn)], dat_s7[, names(s7_env)], method = "spearman", use = "pairwise.complete.obs") |>
  as.data.frame() |>
  rownames_to_column("metric") |>
  pivot_longer(-metric, names_to = "covariate", values_to = "rho") |>
  mutate(
    metric = factor(s7_conn[metric], levels = rev(unname(s7_conn))),
    covariate = factor(s7_env[covariate], levels = unname(s7_env))
  )

fig_s7 <- ggplot(cors_s7, aes(covariate, metric, fill = rho)) +
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
ggsave(here("output", "figs", "figS7_connectivity_environment.png"), fig_s7, width = 6.4, height = 5.6, dpi = 600, bg = "white")

