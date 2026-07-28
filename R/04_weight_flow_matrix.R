# Weight the flow matrix by source-cell biomass and probability of presence,
# then derive per-cell connectivity metrics.
#
#   flow[i, j]   = export probability, source i -> sink j   (rows = source)
#   weight[i]    = predicted value at source i (avg biomass or P(presence))
#
#   wf[i, j] = weight[i] * flow[i, j]   (each source's outflow scaled by its own
#                                        value = realized larval export)
#
#   in_strength[j] = colSums(wf) = sum_i weight[i] * flow[i, j]   (incoming)
#                    -> value received at sink j, summed over all incoming links
#
#
# Flow rows/cols and the grid are both in id order (1..2340), so they align.
#
# Runs once per flow matrix (2010-2016 and the pooled "all"), so the year-to-year
# runs can be used as a dispersal sensitivity check. The biomass and presence
# weights are the same in every run - only the flow matrix changes.

library(tidyverse)
library(here)
library(sf)
library(igraph)
library(patchwork)

# closeness on a cost-transformed graph: an export probability p becomes a
# distance log(1 / p), so strong links are short (Costa et al. 2017,
# doi:10.1371/journal.pone.0189021). Cells with no links have no paths -> 0
closeness_cost <- function(flow) {
  cost <- flow
  cost[cost > 0] <- log(1 / cost[cost > 0])
  g_cost <- graph_from_adjacency_matrix(cost, weighted = TRUE, mode = "directed")
  cl <- suppressWarnings(closeness(g_cost, mode = "all", normalized = TRUE, cutoff = -1))
  ifelse(is.finite(cl), cl, 0)
}

# scale each source row by its predicted value, then derive per-cell metrics
weight_flow <- function(flow, grid) {
  # log1p(biomass) compresses the biomass hotspots (assumes saturating larval
  # output) while staying non-negative and keeping empty cells at 0 (raw log
  # would give -Inf for zeros and negatives below 1)
  wf_biomass <- sweep(flow, 1, grid$avg_biomass, "*")
  wf_biomass_log <- sweep(flow, 1, log1p(grid$avg_biomass), "*")
  wf_present <- sweep(flow, 1, grid$prob_present, "*")

  grid$biomass_in_strength <- colSums(wf_biomass) # incoming biomass-weighted supply
  grid$log_biomass_in_strength <- colSums(wf_biomass_log) # incoming, log1p-biomass weighted
  grid$presence_in_strength <- colSums(wf_present) # incoming presence-weighted supply

  # structural metrics from the flow matrix alone (no biomass), so they are
  # independent of the SDM covariates. Directed graph with edge i -> j = flow[i, j]
  g <- graph_from_adjacency_matrix(flow, weighted = TRUE, mode = "directed")
  grid$deg_in <- degree(g, mode = "in") # number of incoming links, unweighted
  grid$in_strength <- strength(g, mode = "in") # total incoming export prob
  grid$local_retention <- diag(flow) # self-recruitment probability
  grid$eigen_centrality <- suppressWarnings(eigen_centrality(g, directed = TRUE, scale = TRUE)$vector)
  grid$closeness_centrality <- closeness_cost(flow)

  list(
    grid = grid,
    wf = list(biomass = wf_biomass, biomass_log = wf_biomass_log, present = wf_present)
  )
}

# 01 Load inputs ----
grid <- readRDS(here("data", "final", "avg_biomass_grid.rds")) # id order, has predictions
labels <- c(as.character(2010:2016), "all")

# 02 Weight each flow matrix ----
# saved per label inside the loop so the large weighted matrices are not all
# held in memory at once
conn <- map(set_names(labels), function(label) {
  flow <- readRDS(here("data", "intermediate", paste0("flow_matrix_", label, ".rds")))

  stopifnot(
    nrow(flow) == ncol(flow),
    nrow(grid) == nrow(flow),
    identical(grid$id, seq_len(nrow(flow)))
  )

  out <- weight_flow(flow, grid)
  saveRDS(out$grid, here("data", "final", paste0("connectivity_weighted_", label, ".rds")))
  saveRDS(out$wf, here("data", "intermediate", paste0("weighted_flow_matrix_", label, ".rds")))
  out$grid
})

# 03 Combined long table for sensitivity ----
# one row per cell per run, so the year-to-year spread can be inspected directly
conn_years <- conn |>
  list_rbind(names_to = "run")
saveRDS(conn_years, here("data", "final", "connectivity_by_year.rds"))

# 04 Report spread across years ----
# how stable is each metric across the year-specific dispersal runs?
conn_years |>
  filter(run != "all", !is.na(depth)) |>
  pivot_longer(
    c(
      biomass_in_strength, log_biomass_in_strength, presence_in_strength,
      deg_in, in_strength, local_retention, eigen_centrality, closeness_centrality
    ),
    names_to = "metric"
  ) |>
  summarise(mean = mean(value), sd = sd(value), .by = c(metric, id)) |>
  summarise(median_cv = median(sd / mean, na.rm = TRUE), .by = metric) |>
  arrange(median_cv) |>
  print()


