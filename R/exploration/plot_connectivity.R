# Where does each release cell send most of its larvae? One arrow per release
# cell to its single strongest destination in the raw MIKE connectivity
# (settlement probability = agents settling / agents released), self-retention
# left out. Each arrow follows the shortest path through water cells (2 km
# cells touching Limfjorden, 8-neighbour moves) and is smoothed, so it bends
# around the coast instead of crossing land - it is still a summary of the
# link, not a simulated larval trajectory. Circles mark "hub" cells that are
# the main destination of several release cells.

library(tidyverse)
library(here)
library(sf)
library(igraph)

n_col <- 65
n_row <- 36
n_cell <- n_col * n_row

# line-of-sight straightening: from each kept cell, jump to the farthest later
# cell on the path that a straight line reaches without leaving the water, so
# the path only bends where the coast forces it to (not at every 45 degree
# grid step)
straighten <- function(x, y, water_zone) {
  keep <- 1
  i <- 1
  k <- length(x)
  while (i < k) {
    j <- (i + 1):k
    segs <- st_sfc(map(j, \(jj) st_linestring(rbind(c(x[i], y[i]), c(x[jj], y[jj])))), crs = 32632)
    ok <- as.vector(st_within(segs, water_zone, sparse = FALSE))
    i <- max(j[ok | j == i + 1])
    keep <- c(keep, i)
  }
  keep
}

# Chaikin corner cutting: rounds the remaining corners into a smooth curve,
# keeping both end points fixed
chaikin <- function(x, y, n_iter = 3) {
  for (i in seq_len(n_iter)) {
    k <- length(x)
    if (k < 3) break
    x <- c(x[1], as.vector(rbind(0.75 * x[-k] + 0.25 * x[-1], 0.25 * x[-k] + 0.75 * x[-1])), x[k])
    y <- c(y[1], as.vector(rbind(0.75 * y[-k] + 0.25 * y[-1], 0.25 * y[-k] + 0.75 * y[-1])), y[k])
  }
  tibble(x = x, y = y)
}

# 01 Load data ----
cmn <- read_csv(here("data", "connectivity", "raw", "ABM_cmn_all.csv"), col_names = FALSE, show_col_types = FALSE) |>
  as.matrix()
dimnames(cmn) <- NULL
release <- read_csv(here("data", "connectivity", "raw", "ABM_rel_all.csv"), col_names = FALSE, show_col_types = FALSE)[[1]]
grid <- readRDS(here("data", "grid", "grid_env.rds"))

water <- st_read(here("data", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

# 02 Flow matrix ----
flow <- ifelse(matrix(release, n_cell, n_cell) > 0, cmn / matrix(release, n_cell, n_cell), 0)
diag(flow) <- 0

# 03 Strongest destination per release cell ----
idx <- which(flow > 0, arr.ind = TRUE)
links <- tibble(from = idx[, 1], to = idx[, 2], p = flow[idx]) |>
  slice_max(p, n = 1, by = from, with_ties = FALSE) |>
  mutate(link = row_number())

hubs <- links |>
  count(to, name = "n_sources") |>
  filter(n_sources >= 5) |>
  mutate(x = grid$x[to], y = grid$y[to])

# 04 Water network ----
# cells whose 2 km square touches water, plus every link end point
cell_squares <- st_as_sf(grid, coords = c("x", "y"), crs = 32632) |>
  st_buffer(1000, endCapStyle = "SQUARE")
in_water <- lengths(st_intersects(cell_squares, water)) > 0
in_water[c(links$from, links$to)] <- TRUE
wet <- which(in_water)

neighbours <- expand_grid(dr = -1:1, dc = -1:1) |>
  filter(!(dr == 0 & dc == 0))

edges <- map2(neighbours$dr, neighbours$dc, \(dr, dc) {
  tibble(from = wet, row = grid$row[wet] + dr, col = grid$col[wet] + dc)
}) |>
  list_rbind() |>
  filter(between(row, 1, n_row), between(col, 1, n_col)) |>
  mutate(to = col + (row - 1) * n_col) |>
  filter(in_water[to]) |>
  mutate(weight = sqrt((grid$x[from] - grid$x[to])^2 + (grid$y[from] - grid$y[to])^2))

g <- graph_from_data_frame(
  transmute(edges, from = as.character(from), to = as.character(to), weight),
  directed = FALSE, vertices = tibble(name = as.character(wet))
) |>
  simplify(edge.attr.comb = "min")

# 05 Route each link through water ----
# straight lines may use any water within half a cell of the coast, since
# coastal cell centres themselves often sit just on land
water_zone <- water |>
  st_union() |>
  st_simplify(dTolerance = 200) |>
  st_buffer(1000)

paths <- pmap(links, \(from, to, p, link) {
  cells <- as.integer(names(shortest_paths(g, as.character(from), as.character(to))$vpath[[1]]))
  if (length(cells) < 2) cells <- c(from, to) # unreachable: fall back to a straight line
  cells <- cells[straighten(grid$x[cells], grid$y[cells], water_zone)]
  chaikin(grid$x[cells], grid$y[cells]) |>
    mutate(link = link, p = p)
}) |>
  list_rbind() |>
  arrange(p, link) # strongest drawn last, on top

# 06 Plot ----
xlim <- range(paths$x) + c(-5000, 5000)
ylim <- range(paths$y) + c(-5000, 5000)

p <- ggplot() +
  geom_sf(data = water, fill = "#eef3f7", colour = "grey60", linewidth = 0.3) +
  geom_path(
    data = paths, aes(x, y, group = link, colour = p),
    linewidth = 0.35, alpha = 0.5, lineend = "round",
    arrow = arrow(length = unit(1.3, "mm"), type = "closed")
  ) +
  geom_point(data = hubs, aes(x, y, size = n_sources), shape = 21, colour = "black", fill = NA, stroke = 0.7) +
  scale_colour_viridis_c(
    option = "magma", direction = -1, end = 0.9, trans = "log10",
    breaks = c(0.001, 0.01, 0.1), labels = c("0.001", "0.01", "0.1"), name = "Settlement\nprobability"
  ) +
  scale_size_area(max_size = 9, name = "Release cells\nsending most\nlarvae here") +
  coord_sf(crs = 32632, xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    title = "Where each release cell sends most of its larvae",
    subtitle = "One arrow per 2 km release cell to its strongest destination, routed through water (self-retention excluded)"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(colour = "grey15"),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    plot.margin = margin(8, 8, 8, 8)
  )
p

# 07 Save ----
ggsave(here("R", "exploration", "plot_connectivity.png"), p, width = 9, height = 6.5, dpi = 300, bg = "white")
