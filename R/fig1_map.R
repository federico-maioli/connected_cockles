# Figure 1: (a) study area / survey biomass and (b) larval settlement footprints.
# Panel a: biomass predicted by the full model (Stock survey) on the 2 x 2 km
#   grid cells with survey stations, with a location inset top-left.
# Panel b: for four release cells spread over the basin, the full settlement
#   probability field from the pooled flow matrix. Showing where larvae actually
#   land makes the reach of dispersal visible directly, rather than asserting it
#   through a community-detection label.

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
library(rcartocolor)
library(ggspatial)
library(patchwork)

biomass_lab <- expression(Predicted ~ biomass ~ (g / m^2))
biomass_breaks <- c(1, 10, 100, 1000)
# predictions below this (g/m2) are shown as 0
biomass_zero <- 1
# CARTO companion ramps (warm biomass, cool settlement), matched in lightness;
# settlement starts from white so thin, low-probability spread fades out
biomass_cols <- carto_pal(7, "BurgYl")
settle_cols <- c("white", carto_pal(7, "TealGrn"))

land_fill <- "#e6e6e6"
land_line <- "#bcbcbc"

# 01 Shared data ----
dat <- readRDS(here("data", "cockles", "derived", "cockles_env.rds")) |>
  filter(!is.na(biomass), !is.na(x_utm), !is.na(y_utm))
coastline <- st_read(here("data", "boundaries", "land_small_utm", "land_small_utm.shp"), quiet = TRUE) |>
  st_make_valid()
grid <- readRDS(here("data", "grid", "grid_env.rds"))
# Limfjord water outline for the fjord shape (the main-map coastline)
water <- st_read(here("data", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

# pooled flow matrix: agents settling at j / agents released at i, as in
# 05_weight_connectivity.R
n_cell <- nrow(grid)
cmn <- read_csv(here("data", "connectivity", "raw", "ABM_cmn_all.csv"), col_names = FALSE, show_col_types = FALSE) |>
  as.matrix()
dimnames(cmn) <- NULL
release <- read_csv(here("data", "connectivity", "raw", "ABM_rel_all.csv"), col_names = FALSE, show_col_types = FALSE)[[1]]
flow <- ifelse(matrix(release, n_cell, n_cell) > 0, cmn / matrix(release, n_cell, n_cell), 0)

xlim <- range(dat$x_utm * 1000) + c(-16000, 16000)
ylim <- range(dat$y_utm * 1000) + c(-6000, 6000)

# white mask = the plot frame minus the fjord, drawn over the tiles so the 2 km
# cells are clipped to the coastline instead of spilling onto land
frame <- st_as_sfc(st_bbox(
  c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
  crs = st_crs(32632)
))
land_mask <- st_difference(frame, st_union(water))

# 02 Panel A: predicted biomass ----
# expected biomass (presence probability x biomass where present) from the
# full delta-gamma model (Space + Environment + Connectivity), predicted at the
# centre of every 2 km grid cell holding at least one survey station, for the
# Stock survey in the most recent year (2025). Grid covariates are standardised with the model data's own
# mean and SD, recovered from the raw and standardised columns of the fit
fit <- readRDS(here("data", "sdm", "main", "space_env_conn.rds"))

cells <- dat |>
  mutate(
    col = floor((x_utm * 1000 - 450074) / 2000) + 1,
    row = floor((y_utm * 1000 - 6258093) / 2000) + 1,
    id = col + (row - 1) * 65
  ) |>
  summarise(biomass = mean(biomass), .by = id)

std_vars <- c("depth", "temp", "sal", "oxy", "shear_max", "conn_in_strength")
newdata <- readRDS(here("data", "grid", "grid_env_conn.rds")) |>
  filter(id %in% cells$id) |>
  mutate(
    conn_in_strength = conn_in_strength_presence,
    x_utm = x / 1000, y_utm = y / 1000,
    survey = factor("Stock", levels = levels(fit$data$survey)),
    year = factor(last(levels(fit$data$year)), levels = levels(fit$data$year))
  )
for (v in std_vars) {
  s <- sd(fit$data[[v]]) / sd(fit$data[[paste0(v, "_std")]])
  m <- mean(fit$data[[v]]) - s * mean(fit$data[[paste0(v, "_std")]])
  newdata[[paste0(v, "_std")]] <- (newdata[[v]] - m) / s
}
newdata <- filter(newdata, if_all(all_of(paste0(std_vars, "_std")), \(x) !is.na(x)))

pred <- predict(fit, newdata = newdata, type = "response") |>
  mutate(est = if_else(est < biomass_zero, 0, est))
summary(pred$est)

p_main <- ggplot() +
  # zero cells get their own legend key through a dummy linetype mapping (the
  # tile outline itself is not drawn, linewidth = 0)
  geom_tile(data = filter(pred, est == 0), aes(x, y, linetype = "0"), fill = "grey85", linewidth = 0, width = 2000, height = 2000) +
  geom_tile(data = filter(pred, est > 0), aes(x, y, fill = est), width = 2000, height = 2000) +
  geom_sf(data = land_mask, fill = "white", colour = NA) +
  geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.3) +
  scale_fill_gradientn(
    colours = biomass_cols, trans = "pseudo_log", name = biomass_lab,
    breaks = biomass_breaks, labels = scales::comma(biomass_breaks),
    limits = c(biomass_zero, NA) # bar starts at the zero cut-off, so its first label is 1
  ) +
  scale_linetype_manual(values = "solid", name = NULL) +
  annotation_scale(location = "bl", width_hint = 0.2, height = unit(0.15, "cm"), text_cex = 0.7, line_col = "grey40", text_col = "grey40") +
  annotation_north_arrow(location = "br", which_north = "true", height = unit(0.9, "cm"), width = unit(0.7, "cm"), style = north_arrow_minimal(line_col = "grey40", text_col = "grey40", fill = "grey40")) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  theme_void(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  # the "0" key sits just left of the colour bar: same height as the bar, label
  # underneath, so it reads as the bar's first step
  guides(
    linetype = guide_legend(
      order = 1, label.position = "bottom",
      override.aes = list(fill = "grey85", linewidth = 0),
      theme = theme(legend.key.height = unit(3.5, "mm"), legend.key.width = unit(6, "mm"))
    ),
    fill = guide_colourbar(
      title.position = "top", title.hjust = 0.5, order = 2,
      theme = theme(legend.key.height = unit(3.5, "mm"), legend.key.width = unit(55, "mm"))
    )
  )

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

# inset added with patchwork (not cowplot) so panel a stays a regular plot
# that patchwork can align with panel b and collect legends from
panel_a <- p_main +
  inset_element(p_inset, left = 0.01, bottom = 0.665, right = 0.37, top = 0.995, ignore_tag = TRUE)

# 03 Panel B: settlement footprints ----
# Release cells are chosen systematically, not hand-picked: among cells where
# cockles were found (mean survey biomass > 0) that also release
# larvae in the model, three are taken at the 12th, 50th and 88th percentile of
# along-fjord position, plus the westernmost cell of the northern arm (top 20%
# by northing) - the main axis alone leaves the whole northern lobe
# unrepresented. Panels are then ordered west to east.
out_strength <- rowSums(flow)
src_pool <- intersect(cells$id[cells$biomass > 0], which(out_strength > 0))

src_axis <- quantile(grid$x[src_pool], c(0.12, 0.5, 0.88)) |>
  sapply(\(q) src_pool[which.min(abs(grid$x[src_pool] - q))])
north_arm <- src_pool[grid$y[src_pool] > quantile(grid$y[src_pool], 0.80)]
src_north <- north_arm[which.min(grid$x[north_arm])]

srcs <- c(src_axis, src_north)
srcs <- srcs[order(grid$x[srcs])] # west -> east

# panel titles, in west-to-east order; swap in the local basin names here
src_labs <- c("Release 1", "Release 2", "Release 3", "Release 4")

fields <- map2(srcs, src_labs, \(i, l) {
  tibble(x = grid$x, y = grid$y, p = flow[i, ], basin = l)
}) |>
  bind_rows() |>
  filter(p > 0) |>
  mutate(basin = factor(basin, levels = src_labs))

src_pts <- st_as_sf(
  tibble(basin = factor(src_labs, levels = src_labs), x = grid$x[srcs], y = grid$y[srcs]),
  coords = c("x", "y"), crs = 32632
)

# area holding 90% of each release's settlement, quoted in the caption
walk2(srcs, src_labs, \(i, l) {
  w <- flow[i, ]
  n <- which(cumsum(sort(w, decreasing = TRUE)) / sum(w) >= 0.9)[1]
  cat(sprintf("%-15s 90%% of larvae settle within %3d cells (%4.0f km2)\n", l, n, n * 4))
})

# tiles first, then the mask, then the outline: the raster sits under the
# coastline. A white-to-teal ramp keeps the top of the scale clear of the black
# release marker, which is drawn white-filled so it reads over any cell value
panel_b <- ggplot() +
  geom_tile(data = fields, aes(x, y, fill = p), width = 2000, height = 2000) +
  geom_sf(data = land_mask, fill = "white", colour = NA) +
  geom_sf(data = water, fill = NA, colour = "grey70", linewidth = 0.2) +
  geom_sf(data = src_pts, fill = "white", colour = "grey10", size = 2.2, stroke = 0.7, shape = 23) +
  scale_fill_gradientn(colours = settle_cols, trans = "sqrt", breaks = c(0.01, 0.03, 0.06), name = "Settlement probability") +
  facet_wrap(~basin, ncol = 1) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(size = 10, colour = "grey20", margin = margin(2, 0, 2, 0))
  ) +
  guides(fill = guide_colourbar(
    title.position = "top", title.hjust = 0.5, order = 3,
    theme = theme(legend.key.height = unit(3.5, "mm"), legend.key.width = unit(40, "mm"))
  ))

# 04 Combine and save ----
# panel a on the left, the four footprints stacked in a column on the right,
# aligned top and bottom; each legend sits under its own panel, on one row
# (patchwork aligns the two panels). Widths are set so the column of four
# small maps (1.33:1 each) is about as tall as the main map
fig1 <- panel_a + panel_b +
  plot_layout(widths = c(3.4, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(
    plot.tag = element_text(size = 12, face = "bold", colour = "grey20"),
    legend.position = "bottom",
    legend.box.just = "bottom",
    legend.spacing.x = unit(2, "mm"),
    legend.title = element_text(size = 11, colour = "grey15"),
    legend.text = element_text(size = 10, colour = "grey30")
  )

dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "fig1_map.png"), fig1, width = 11, height = 7.6, dpi = 600, bg = "white")
