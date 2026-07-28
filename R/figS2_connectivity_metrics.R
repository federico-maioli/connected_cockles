# Supplementary figure: connectivity metrics mapped on the 2 km grid.
#
# Five metrics from the pooled ("all") flow matrix, one panel each, styled like
# Figure 1 and the covariate maps (white background, water outline only). Each
# metric keeps its own single-hue sequential scale (light = low, dark = high)
# because the metrics are on very different scales.
#
# Only wet cells (non-missing depth) are drawn; cells with no larval links sit
# at 0, which is the lightest end of each ramp.

library(tidyverse)
library(here)
library(sf)
library(patchwork)

# single-hue sequential ramp per metric, light (low) -> dark (high)
metric_scale <- function(col, tr) {
  switch(col,
    log_biomass_in_strength = scale_fill_gradient(low = "#e5f5e0", high = "#005a32", name = NULL, transform = tr),
    deg_in = scale_fill_gradient(low = "#d0e1f2", high = "#08306b", name = NULL, transform = tr),
    in_strength = scale_fill_gradient(low = "#efedf5", high = "#3f007d", name = NULL, transform = tr),
    eigen_centrality = scale_fill_gradient(low = "#fee6ce", high = "#7f2704", name = NULL, transform = tr),
    closeness_centrality = scale_fill_gradient(low = "#dbf1f0", high = "#01665e", name = NULL, transform = tr)
  )
}

metric_panel <- function(col) {
  ggplot() +
    geom_tile(data = grid, aes(x_utm, y_utm, fill = .data[[col]])) +
    geom_sf(data = mask, fill = "white", colour = NA) +
    geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.22) +
    metric_scale(col, transforms[[col]]) +
    coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
    labs(title = titles[[col]]) +
    theme_void(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = 9.5, hjust = 0.5, colour = "grey20"),
      legend.position = "right",
      legend.key.width = unit(3, "mm"), legend.key.height = unit(8, "mm"),
      legend.text = element_text(size = 7, colour = "grey30"),
      plot.tag = element_text(size = 12, face = "bold", colour = "grey20"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# 01 Data ----
water <- st_read(here("data", "raw", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

# pooled flow matrix run; only wet cells carry a prediction
grid <- readRDS(here("data", "derived", "connectivity_weighted_all.rds")) |>
  filter(!is.na(depth))

# keep only cells inside Limfjorden: the ABM domain extends past the fjord into
# the North Sea, and those outer cells have no larval links (129 of the 130
# cells outside the polygon have deg_in = 0), so drawing them implies data where
# there is none
inside <- lengths(st_intersects(
  st_as_sf(grid, coords = c("x_utm", "y_utm"), crs = 32632),
  water
)) > 0
grid <- grid[inside, ]

# frame on the mapped cells (two cells of margin) rather than the full polygon,
# which reaches further east than any grid cell
xlim <- range(grid$x_utm) + c(-4000, 4000)
ylim <- range(grid$y_utm) + c(-4000, 4000)

# white mask = the plot frame minus the fjord, drawn over the tiles so the 2 km
# cells are clipped to the coastline instead of spilling onto land
frame <- st_as_sfc(st_bbox(
  c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
  crs = st_crs(32632)
))
mask <- st_difference(frame, st_union(water))

# 02 Metrics to map ----
# skewed metrics get a square-root scale so the low end stays readable;
# sqrt handles the zeros that log / pseudo-log would compress
metrics <- c(
  "log_biomass_in_strength", "deg_in", "in_strength",
  "eigen_centrality", "closeness_centrality"
)
titles <- list(
  log_biomass_in_strength = "Biomass-weighted in-strength (log)",
  deg_in = "In-degree",
  in_strength = "In-strength",
  eigen_centrality = "Eigenvector centrality",
  closeness_centrality = "Closeness centrality"
)
transforms <- c(
  log_biomass_in_strength = "sqrt",
  deg_in = "identity",
  in_strength = "identity",
  eigen_centrality = "sqrt",
  closeness_centrality = "identity"
)

# 03 Build panels ----
panels <- map(metrics, metric_panel)

# 04 Combine and save ----
# 3 columns (3 on top, 2 below); tag panels a-e
fig <- wrap_plots(panels, ncol = 3) +
  plot_annotation(tag_levels = "a")

dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(
  here("output", "figs", "figS2_connectivity_metrics.png"), fig,
  width = 13, height = 5.6, dpi = 600, bg = "white"
)
