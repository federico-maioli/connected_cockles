# Figure 1: (a) study area / survey biomass and (b) larval settlement footprints.
# Panel a: survey biomass over the coastline, with a location inset top-left.
# Panel b: for four release cells spread over the basin, the full settlement
#   probability field from the pooled flow matrix. Showing where larvae actually
#   land makes the reach of dispersal visible directly, rather than asserting it
#   through a community-detection label.

library(tidyverse)
library(here)
library(sf)
library(rcartocolor)
library(ggspatial)
library(cowplot)
library(patchwork)

biomass_lab <- expression(Biomass ~ (g / m^2))
biomass_breaks <- c(0, 100, 1000, 10000)
biomass_cols <- rev(carto_pal(7, "SunsetDark"))

land_fill <- "#e6e6e6"
land_line <- "#bcbcbc"

# 01 Shared data ----
dat <- readRDS(here("data", "derived", "cockles_connectivity.rds")) |>
  filter(!is.na(biomass), !is.na(x_utm), !is.na(y_utm))
coastline <- st_read(here("data", "raw", "boundaries", "land_small_utm", "land_small_utm.shp"), quiet = TRUE) |>
  st_make_valid()
grid <- readRDS(here("data", "derived", "connectivity_weighted_all.rds"))
flow <- readRDS(here("data", "derived", "flow_matrix_all.rds"))
# Limfjord water outline for the fjord shape (the main-map coastline)
water <- st_read(here("data", "raw", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

pts <- dat |>
  mutate(x_m = x_utm * 1000, y_m = y_utm * 1000) |>
  arrange(biomass) |>
  st_as_sf(coords = c("x_m", "y_m"), crs = 32632)
pts_abs <- filter(pts, biomass == 0)
pts_pres <- filter(pts, biomass > 0)

xlim <- range(dat$x_utm * 1000) + c(-16000, 16000)
ylim <- range(dat$y_utm * 1000) + c(-6000, 6000)

# 02 Panel A: survey biomass ----
p_main <- ggplot() +
  geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.3) +
  geom_sf(data = pts_abs, colour = "grey78", size = 0.3, alpha = 0.6) +
  geom_sf(data = pts_pres, aes(colour = biomass), size = 1.1, alpha = 0.85) +
  scale_colour_gradientn(
    colours = biomass_cols, trans = "pseudo_log", name = biomass_lab,
    breaks = biomass_breaks, labels = scales::comma(biomass_breaks)
  ) +
  annotation_scale(location = "bl", width_hint = 0.2, height = unit(0.15, "cm"), text_cex = 0.7, line_col = "grey40", text_col = "grey40") +
  annotation_north_arrow(location = "br", which_north = "true", height = unit(0.9, "cm"), width = unit(0.7, "cm"), style = north_arrow_minimal(line_col = "grey40", text_col = "grey40", fill = "grey40")) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "inside", legend.position.inside = c(0.9, 0.58),
    legend.title = element_text(size = 9, colour = "grey20"),
    legend.text = element_text(size = 8, colour = "grey30"),
    legend.key.height = unit(9, "mm"), legend.key.width = unit(4, "mm")
  ) +
  guides(colour = guide_colourbar(title.position = "top"))

# location inset, top-left: where the basin sits in the North Sea transition
inset_land <- st_transform(coastline, 4326)
# label positions are resolved here rather than by geom_sf_text, which would
# recompute them at render time and warn about lon/lat input on every run
country_labs <- inset_land |>
  filter(NAME %in% c("Denmark", "Germany", "Sweden", "Norway")) |>
  group_by(NAME) |>
  summarise(.groups = "drop") |>
  st_crop(xmin = 4.6, ymin = 54, xmax = 12.4, ymax = 57.7) |>
  st_point_on_surface() |>
  suppressWarnings()
country_labs <- bind_cols(
  st_drop_geometry(country_labs),
  as_tibble(st_coordinates(country_labs))
)
# the red box marks the extent actually drawn in the panels (xlim/ylim), not just
# the survey-point bounding box, so it matches what the other figures show
study_box <- st_bbox(
  c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
  crs = st_crs(32632)
) |>
  st_as_sfc() |>
  st_transform(4326)

p_inset <- ggplot() +
  geom_sf(data = inset_land, fill = land_fill, colour = land_line, linewidth = 0.2) +
  geom_text(data = country_labs, aes(X, Y, label = NAME), size = 2.3, colour = "grey35") +
  geom_sf(data = study_box, fill = NA, colour = "#c0392b", linewidth = 0.8) +
  coord_sf(xlim = c(4, 13), ylim = c(53.5, 58), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey40", linewidth = 0.6)
  )

panel_a <- ggdraw() +
  draw_plot(p_main) +
  draw_plot(p_inset, x = 0.02, y = 0.70, width = 0.28, height = 0.28)

# 03 Panel B: settlement footprints ----
# Release cells are chosen systematically, not hand-picked: among cells whose
# total export exceeds the median, three are taken at the 12th, 50th and 88th
# percentile of along-fjord position, plus the westernmost cell of the northern
# arm (top 20% by northing) - the main axis alone leaves the whole northern lobe
# unrepresented. Panels are then ordered west to east.
out_strength <- rowSums(flow)
src_pool <- which(out_strength > quantile(out_strength[out_strength > 0], 0.5))

src_axis <- quantile(grid$x_utm[src_pool], c(0.12, 0.5, 0.88)) |>
  sapply(\(q) src_pool[which.min(abs(grid$x_utm[src_pool] - q))])
north_arm <- src_pool[grid$y_utm[src_pool] > quantile(grid$y_utm[src_pool], 0.80)]
src_north <- north_arm[which.min(grid$x_utm[north_arm])]

srcs <- c(src_axis, src_north)
srcs <- srcs[order(grid$x_utm[srcs])] # west -> east

# panel titles, in west-to-east order; swap in the local basin names here
src_labs <- c("Release 1", "Release 2", "Release 3", "Release 4")

fields <- map2(srcs, src_labs, \(i, l) {
  tibble(x = grid$x_utm, y = grid$y_utm, p = flow[i, ], basin = l)
}) |>
  bind_rows() |>
  filter(p > 0) |>
  mutate(basin = factor(basin, levels = src_labs))

src_pts <- st_as_sf(
  tibble(basin = factor(src_labs, levels = src_labs), x = grid$x_utm[srcs], y = grid$y_utm[srcs]),
  coords = c("x", "y"), crs = 32632
)

# area holding 90% of each release's settlement, quoted in the caption
walk2(srcs, src_labs, \(i, l) {
  w <- flow[i, ]
  n <- which(cumsum(sort(w, decreasing = TRUE)) / sum(w) >= 0.9)[1]
  cat(sprintf("%-15s 90%% of larvae settle within %3d cells (%4.0f km2)\n", l, n, n * 4))
})

# white mask = the plot frame minus the fjord, drawn over the tiles so the 2 km
# cells are clipped to the coastline instead of spilling onto land
frame <- st_as_sfc(st_bbox(
  c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
  crs = st_crs(32632)
))
land_mask <- st_difference(frame, st_union(water))

# tiles first, then the mask, then the outline: the raster sits under the
# coastline. A light-to-blue ramp keeps the top of the scale clear of the black
# release marker, which is drawn white-filled so it reads over any cell value
panel_b <- ggplot() +
  geom_tile(data = fields, aes(x, y, fill = p), width = 2000, height = 2000) +
  geom_sf(data = land_mask, fill = "white", colour = NA) +
  geom_sf(data = water, fill = NA, colour = "grey70", linewidth = 0.2) +
  geom_sf(data = src_pts, fill = "white", colour = "grey10", size = 2.2, stroke = 0.7, shape = 23) +
  scale_fill_distiller(palette = "Blues", direction = 1, trans = "sqrt", name = "Settlement probability") +
  facet_wrap(~basin, nrow = 1) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(size = 9.5, colour = "grey20"),
    legend.position = "bottom",
    legend.title = element_text(size = 9, colour = "grey20"),
    legend.text = element_text(size = 8, colour = "grey30"),
    legend.key.height = unit(3, "mm"), legend.key.width = unit(12, "mm")
  )

# 04 Combine and save ----
# heights match the map aspect (97 x 72 km, 1.33:1) so neither panel carries dead
# space; stacking keeps the maps large once the figure is scaled to page width.
# the tag theme is set per panel rather than through plot_annotation(), which
# warns under patchwork 1.3.1 + ggplot2 4.0
tag_theme <- theme(plot.tag = element_text(size = 12, face = "bold", colour = "grey20"))

fig1 <- (wrap_elements(panel_a) + tag_theme) / (panel_b + tag_theme) +
  plot_layout(heights = c(2.6, 1)) +
  plot_annotation(tag_levels = "a")

dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "fig1_map.png"), fig1, width = 9, height = 9.4, dpi = 600, bg = "white")
