# Weight the raw MIKE connectivity matrix by modeled habitat suitability
# (presence probability and log biomass, from 03_fit_suitability.R) and
# compute graph metrics on the weighted and unweighted matrices: in-degree,
# in-strength, in-closeness, eigenvector centrality, and local transitivity.
# Also computes in-strength with the presence-weighted matrix's edges
# reversed (in_strength_flipped), to check whether direction matters.
#
# Ported from Flemming Thorbjorn Hansen's original connectivity-analysis
# script (2024-06-01): the flow-matrix construction (row i / release[i], as
# in R/old/01_compute_flow_matrix.R), the -log(probability) distance
# transform used for closeness, and the graph-metric definitions all follow
# that script exactly - this is also how ConMetrics_PresAbs.csv's own
# columns were generated. Source-node weighting (habitat coverage there,
# suitability here) is applied the same way: row i of the flow matrix is
# multiplied by weight[i].
#
# Saved as .asc rasters (same grid as data/env/Asc4dtuaqua), one per metric x
# weighting combination, so they can be extracted onto survey points the
# same way as the environmental covariates.

library(tidyverse)
library(here)
library(sf)
library(terra)
library(igraph)

n_col <- 65
n_row <- 36
n_cell <- n_col * n_row

# in-degree, in-strength, in-closeness, eigen, and local transitivity for one
# (weighted or unweighted) flow matrix, mapped back onto the full n_cell grid
# (excluded nodes get 0, matching the legacy script)
compute_metrics <- function(mat) {
  keep <- rowSums(mat) != 0
  sub <- mat[keep, keep, drop = FALSE]

  g <- graph_from_adjacency_matrix(sub, weighted = TRUE, mode = "directed")
  in_degree <- degree(g, mode = "in")
  in_strength <- strength(g, mode = "in")
  in_strength <- in_strength / sum(E(g)$weight)
  eigen <- eigen_centrality(g, directed = TRUE)$vector
  transitivity <- transitivity(g, type = "local")
  transitivity[is.na(transitivity) | is.infinite(transitivity)] <- 0

  # -log(probability) distance transform, as in the legacy script's
  # betweenness/closeness section
  dist <- -log(sub)
  dist[!is.finite(dist)] <- 0
  g_dist <- graph_from_adjacency_matrix(dist, weighted = TRUE, mode = "directed", diag = FALSE)
  in_closeness <- closeness(g_dist, mode = "in", normalized = TRUE, cutoff = -1)

  fill <- function(x) {
    full <- numeric(n_cell)
    full[keep] <- x
    full
  }
  tibble(
    id = 1:n_cell,
    in_degree = fill(in_degree),
    in_strength = fill(in_strength),
    in_closeness = fill(in_closeness),
    eigen = fill(eigen),
    transitivity = fill(transitivity)
  )
}

# in-strength on the transposed matrix (edges reversed)
in_strength_flipped <- function(mat) {
  keep <- rowSums(mat) != 0
  sub <- t(mat)[keep, keep, drop = FALSE]
  g <- graph_from_adjacency_matrix(sub, weighted = TRUE, mode = "directed")
  s <- strength(g, mode = "in")
  s <- s / sum(E(g)$weight)
  full <- numeric(n_cell)
  full[keep] <- s
  full
}

# 01 Load data ----
cmn <- read_csv(
  here("data", "connectivity", "raw", "ABM_cmn_all.csv"),
  col_names = FALSE, na = "1e-35", show_col_types = FALSE
) |>
  as.matrix()
cmn[is.na(cmn)] <- 0
dimnames(cmn) <- NULL

release <- read_csv(
  here("data", "connectivity", "raw", "ABM_rel_all.csv"),
  col_names = FALSE, show_col_types = FALSE
)[[1]]

grid <- readRDS(here("data", "grid", "grid_suitability.rds"))
stopifnot(
  nrow(cmn) == n_cell, ncol(cmn) == n_cell, length(release) == n_cell,
  all(grid$id == 1:n_cell)
)

# 02 Flow matrix ----
# export probability, source -> sink, relative to total release:
# flow[i, j] = cmn[i, j] / release[i]
denom <- matrix(release, nrow = n_cell, ncol = n_cell, byrow = FALSE)
flow <- ifelse(denom > 0, cmn / denom, 0)
flow[is.nan(flow)] <- 0

# 03 Source weights ----
# grid is ordered 1:n_cell by id, matching the matrix's row/column order.
# log(biomass) is min-max normalised to [0, 1] over the predicted cells before
# use, so it stays non-negative (log(biomass) alone goes negative below
# biomass = 1, which produced negative edge weights and complex eigenpairs)
weight_presence <- replace_na(grid$suit_presence, 0)

log_biomass <- log(grid$suit_biomass)
weight_biomass <- (log_biomass - min(log_biomass, na.rm = TRUE)) /
  (max(log_biomass, na.rm = TRUE) - min(log_biomass, na.rm = TRUE))
weight_biomass <- replace_na(weight_biomass, 0)

# 04 Weighted matrices ----
# `flow * w` recycles w down the rows, i.e. multiplies row i by w[i]
flow_presence <- flow * weight_presence
flow_biomass <- flow * weight_biomass

# 05 Compute metrics ----
metrics_unweighted <- compute_metrics(flow) |> rename_with(~ paste0(.x, "_unweighted"), -id)
metrics_presence <- compute_metrics(flow_presence) |> rename_with(~ paste0(.x, "_presence"), -id)
metrics_biomass <- compute_metrics(flow_biomass) |> rename_with(~ paste0(.x, "_biomass"), -id)

metrics <- metrics_unweighted |>
  left_join(metrics_presence, by = "id") |>
  left_join(metrics_biomass, by = "id") |>
  mutate(in_strength_flipped = in_strength_flipped(flow_presence))

# 06 Save as rasters ----
out_dir <- here("data", "connectivity", "derived")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (col in setdiff(names(metrics), "id")) {
  df <- tibble(x = grid$x, y = grid$y, value = metrics[[col]])
  r <- rast(as.data.frame(df), type = "xyz", crs = "EPSG:32632")
  writeRaster(r, file.path(out_dir, paste0(col, ".asc")), filetype = "AAIGrid", overwrite = TRUE)
}

# 07 Plot ----
land <- st_read(
  here("data", "boundaries", "land_small_utm", "land_small_utm.shp"),
  quiet = TRUE
) |>
  st_make_valid()

# nodes dropped from the unweighted subgraph (rowSums == 0) are set to NA
# rather than plotted as a real 0
keep_unweighted <- rowSums(flow) != 0

plot_data <- tibble(
  x = grid$x, y = grid$y,
  in_strength = if_else(keep_unweighted, metrics$in_strength_unweighted, NA_real_)
)

ggplot() +
  geom_sf(data = land, fill = "grey80", colour = NA) +
  geom_tile(data = plot_data, aes(x, y, fill = in_strength), width = 2000, height = 2000) +
  coord_sf(crs = 32632, xlim = range(grid$x), ylim = range(grid$y), expand = FALSE) +
  scale_fill_viridis_c(na.value = NA, trans = "pseudo_log", name = "In-strength\n(unweighted)") +
  theme_light()
