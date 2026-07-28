# Figure 1: (A) study area / survey biomass and (B) larval dispersal (flow matrix).
# Panel A: survey biomass over the coastline with a location inset.
# Panel B: strong source -> sink dispersal links drawn as arrows on the same fjord.

library(tidyverse)
library(here)
library(sf)
library(igraph)
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
    legend.position = "inside", legend.position.inside = c(0.9, 0.78),
    legend.title = element_text(size = 9, colour = "grey20"),
    legend.text = element_text(size = 8, colour = "grey30"),
    legend.key.height = unit(9, "mm"), legend.key.width = unit(4, "mm")
  ) +
  guides(colour = guide_colourbar(title.position = "top"))

# location inset
inset_land <- st_transform(coastline, 4326)
country_labs <- inset_land |>
  filter(NAME %in% c("Denmark", "Germany", "Sweden", "Norway")) |>
  group_by(NAME) |>
  summarise(.groups = "drop") |>
  st_crop(xmin = 4.6, ymin = 54, xmax = 12.4, ymax = 57.7) |>
  suppressWarnings()
study_box <- st_as_sfc(st_bbox(st_transform(pts, 4326)))

p_inset <- ggplot() +
  geom_sf(data = inset_land, fill = land_fill, colour = land_line, linewidth = 0.2) +
  geom_sf_text(data = country_labs, aes(label = NAME), size = 2.3, colour = "grey35") +
  geom_sf(data = study_box, fill = NA, colour = "#c0392b", linewidth = 0.8) +
  coord_sf(xlim = c(4, 13), ylim = c(53.5, 58), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey40", linewidth = 0.6)
  )

panel_a <- ggdraw() +
  draw_plot(p_main) +
  draw_plot(p_inset, x = 0.02, y = 0.66, width = 0.31, height = 0.31)

# 03 Panel B: connectivity clusters and exchange ----
# Louvain community detection on the (symmetrised) connectivity matrix delineates
# dispersal-based sub-populations; arrows show inter-cluster larval exchange
# (% of a cluster's export) and circled values the within-cluster self-recruitment (%)
active <- which(rowSums(flow > 0) + colSums(flow > 0) > 0)
fa <- flow[active, active]
gu <- graph_from_adjacency_matrix((fa + t(fa)) / 2, mode = "undirected", weighted = TRUE)
set.seed(1)
cl <- rep(NA_integer_, nrow(flow))
cl[active] <- as.integer(membership(cluster_louvain(gu, weights = E(gu)$weight)))
keep <- as.integer(names(sort(table(cl), decreasing = TRUE)))[1:5]
cl <- ifelse(cl %in% keep, match(cl, keep), NA_integer_)
grid$clus <- factor(cl, levels = 1:5)

# cluster centroids, inter-cluster exchange and self-recruitment
cent <- grid |>
  filter(!is.na(clus)) |>
  group_by(clus) |>
  summarise(x = mean(x_utm), y = mean(y_utm), .groups = "drop") |>
  mutate(ci = as.integer(clus))
idx <- which(flow > 0, arr.ind = TRUE)
lk <- tibble(a = cl[idx[, 1]], b = cl[idx[, 2]], w = flow[idx]) |>
  filter(!is.na(a), !is.na(b)) |>
  group_by(a, b) |>
  summarise(w = sum(w), .groups = "drop")
tot <- lk |> group_by(a) |> summarise(tot = sum(w), .groups = "drop")
sr <- lk |> filter(a == b) |> left_join(tot, by = "a") |> transmute(ci = a, sr = round(100 * w / tot))
exch <- lk |> filter(a != b) |> left_join(tot, by = "a") |> mutate(pct = 100 * w / tot) |> filter(pct > 3) |>
  left_join(select(cent, ci, x, y), by = c("a" = "ci")) |> rename(x0 = x, y0 = y) |>
  left_join(select(cent, ci, x, y), by = c("b" = "ci")) |> rename(x1 = x, y1 = y)
cent <- cent |> left_join(sr, by = "ci")

# shorten exchange arrows so the heads sit in open space, clear of the circles
gap <- 3200
exch <- exch |> mutate(
  ux = x1 - x0, uy = y1 - y0, len = sqrt(ux^2 + uy^2), ux = ux / len, uy = uy / len,
  x0 = x0 + ux * gap, y0 = y0 + uy * gap, x1 = x1 - ux * gap, y1 = y1 - uy * gap
)
# self-recruitment as a small loop above each cluster centroid
loops <- cent |> mutate(lx0 = x - 800, ly0 = y + 2400, lx1 = x + 800, ly1 = y + 2400)

# clip cluster cells to the fjord water (drop any that fall outside the outline)
grid$in_water <- lengths(st_intersects(
  st_as_sf(grid, coords = c("x_utm", "y_utm"), crs = 32632),
  st_buffer(water, 1000)
)) > 0

panel_b <- ggplot() +
  geom_tile(data = filter(grid, !is.na(clus), in_water), aes(x_utm, y_utm, fill = clus), width = 2000, height = 2000) +
  geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.3) +
  geom_curve(
    data = exch, aes(x = x0, y = y0, xend = x1, yend = y1, linewidth = pct),
    curvature = 0.16, colour = "grey15", lineend = "round",
    arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  geom_curve(
    data = loops, aes(x = lx0, y = ly0, xend = lx1, yend = ly1, linewidth = sr),
    curvature = -2.4, colour = "grey15", lineend = "round",
    arrow = arrow(length = unit(2.4, "mm"), type = "closed")
  ) +
  geom_point(data = cent, aes(x, y), size = 7, shape = 21, fill = "white", colour = "grey30", stroke = 0.6) +
  geom_text(data = cent, aes(x, y, label = sr), size = 2.8, fontface = "bold") +
  scale_fill_carto_d(palette = "Safe", guide = "none") +
  scale_linewidth(range = c(0.3, 3.2), name = "Connectivity (%)") +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "inside", legend.position.inside = c(0.92, 0.3),
    legend.title = element_text(size = 9, colour = "grey20"),
    legend.text = element_text(size = 8, colour = "grey30"),
    legend.key.size = unit(4, "mm")
  )

# 04 Combine and save ----
fig1 <- wrap_elements(panel_a) + panel_b +
  plot_annotation(tag_levels = "A")

dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "fig1_map.png"), fig1, width = 16, height = 6.2, dpi = 600, bg = "white")
