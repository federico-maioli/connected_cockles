# Weight the flow matrix by source-cell biomass and probability of presence,
# then derive per-cell connectivity metrics.
#
#   flow[i, j]   = export probability, source i -> sink j   (rows = source)
#   weight[i]    = predicted value at source i (avg biomass or P(presence))
#
#   wf[i, j] = weight[i] * flow[i, j]   (each source's outflow scaled by its own
#                                        value = realized larval export)
#
#   supply[j] = colSums(wf) = sum_i weight[i] * flow[i, j]   (incoming, "import")
#               -> value received at sink j, summed over all incoming links
#   export[i] = rowSums(wf) = weight[i] * total export prob of i   (outgoing)
#
# Flow rows/cols and the grid are both in id order (1..2340), so they align.

library(tidyverse)
library(here)
library(sf)
library(igraph)
library(rnaturalearth)
library(patchwork)

# 01 Load inputs ----
flow <- readRDS(here("data", "intermediate", "flow_matrix.rds"))
grid <- readRDS(here("data", "final", "avg_biomass_grid.rds")) # id order, has predictions

stopifnot(
  nrow(flow) == ncol(flow),
  nrow(grid) == nrow(flow),
  identical(grid$id, seq_len(nrow(flow)))
)

# 02 Weight the flow by source-cell value ----
# scale each source row i by its predicted value; colSums = incoming supply at
# each sink, rowSums = outgoing export from each source.
# log1p(biomass) is an alternative weight that compresses the biomass hotspots
# (assumes saturating larval output) while staying non-negative and keeping empty
# cells at 0 (raw log would give -Inf for zeros and negatives below 1).
wf_biomass <- sweep(flow, 1, grid$avg_biomass, "*")
wf_biomass_log <- sweep(flow, 1, log1p(grid$avg_biomass), "*")
wf_present <- sweep(flow, 1, grid$prob_present, "*")

grid$biomass_supply <- colSums(wf_biomass) # incoming biomass-weighted larval supply
grid$biomass_export <- rowSums(wf_biomass) # outgoing biomass-weighted larval export
grid$biomass_supply_log <- colSums(wf_biomass_log) # incoming, log1p-biomass weighted
grid$biomass_export_log <- rowSums(wf_biomass_log) # outgoing, log1p-biomass weighted
grid$present_supply <- colSums(wf_present) # incoming presence-weighted supply
grid$present_export <- rowSums(wf_present) # outgoing presence-weighted export

# 03 Structural graph metrics ----
# connectivity metrics from the flow matrix alone (no biomass), so they are
# independent of the SDM covariates. Directed graph with edge i -> j = flow[i, j].
g <- graph_from_adjacency_matrix(flow, weighted = TRUE, mode = "directed")
grid$in_strength <- strength(g, mode = "in") # total incoming export prob (import connectivity)
grid$out_strength <- strength(g, mode = "out") # total outgoing export prob
grid$local_retention <- diag(flow) # self-recruitment probability
grid$eigen_centrality <- suppressWarnings(eigen_centrality(g, directed = TRUE, scale = TRUE)$vector)

# 04 Save ----
saveRDS(grid, here("data", "final", "connectivity_weighted.rds"))
saveRDS(
  list(biomass = wf_biomass, biomass_log = wf_biomass_log, present = wf_present),
  here("data", "intermediate", "weighted_flow_matrix.rds")
)

# 05 Plot weighted supply ----
# coastline cropped to the grid area, projected to UTM metres
land <- ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf") |>
  st_make_valid()
region <- st_bbox(st_transform(st_as_sf(grid, coords = c("x_utm", "y_utm"), crs = 32632), 4326))
region["xmin"] <- region["xmin"] - 0.4
region["ymin"] <- region["ymin"] - 0.4
region["xmax"] <- region["xmax"] + 0.4
region["ymax"] <- region["ymax"] + 0.4
land_region <- suppressWarnings(st_crop(land, region)) |> st_transform(32632)

# tiles for the wet cells, coastline drawn on top; each sink cell coloured by the
# sum of incoming links weighted by source biomass / presence
plot_dat <- grid |> filter(!is.na(depth))
xlim <- range(grid$x_utm)
ylim <- range(grid$y_utm)

p_biomass_supply <- ggplot() +
  geom_tile(data = plot_dat, aes(x_utm, y_utm, fill = biomass_supply)) +
  geom_sf(data = land_region, fill = "grey85", colour = "grey60", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  scale_fill_viridis_c() +
  labs(title = "Biomass-weighted supply", x = NULL, y = NULL, fill = "Supply") +
  theme_light()

p_biomass_supply_log <- ggplot() +
  geom_tile(data = plot_dat, aes(x_utm, y_utm, fill = biomass_supply_log)) +
  geom_sf(data = land_region, fill = "grey85", colour = "grey60", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  scale_fill_viridis_c() +
  labs(title = "log1p(biomass)-weighted supply", x = NULL, y = NULL, fill = "Supply") +
  theme_light()

p_present_supply <- ggplot() +
  geom_tile(data = plot_dat, aes(x_utm, y_utm, fill = present_supply)) +
  geom_sf(data = land_region, fill = "grey85", colour = "grey60", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  scale_fill_viridis_c(option = "mako") +
  labs(title = "Presence-weighted supply", x = NULL, y = NULL, fill = "Supply") +
  theme_light()

p_biomass_supply + p_biomass_supply_log + p_present_supply
