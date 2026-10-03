# Build the SPDE mesh used by every spatial model (04_fit_suitability.R,
# 07_fit_sdm.R): resolution varies with the local density of survey stations,
# the coastline is built into the mesh, and land triangles act as a barrier.
#
# Following the fmesher "variable mesh quality" article
# (https://inlabru-org.github.io/fmesher/articles/variable_mesh_quality.html),
# quality.spec gives the maximum triangle edge as a function of location:
#   - at each station, the target edge is the distance to its 20th nearest
#     distinct station location, kept within 0.5-6 km (fine in the densely
#     surveyed beds, coarse on the sparse 2018 fjord-wide grid)
#   - away from stations the target grows by 1 km per km, up to 10 km (land,
#     mesh extension)
# The coastline is an interior constraint, so every triangle is wholly water
# or wholly land, and land triangles get a range of 10% of the water range
# (sdmTMBextra::add_barrier_mesh()).
#
# Saves data/mesh/mesh.rds: the fmesher mesh and the land polygons that define
# the barrier. Model scripts attach it to their data with make_mesh() and
# add_barrier_mesh().

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)

source(here("R", "helpers.R")) # fmesher compatibility shim, needed by add_barrier_mesh() below

# resolution settings (km)
edge_min <- 0.5 # finest target edge, in the most densely sampled beds
edge_max <- 6 # coarsest target edge at a station (sparse 2018 grid)
spacing_k <- 20 # local station spacing = distance to the 20th nearest distinct location
growth <- 1 # target edge grows by 1 km per km away from the nearest station
edge_far <- 10 # coarsest target edge (land, mesh extension)
refine_ratio <- 0.6 # refined edges come out at ~0.6 x the limit, so limit = target / 0.6

# 01 Load data ----
dat <- readRDS(here("data", "cockles", "derived", "cockles_env.rds")) |>
  filter(!is.na(x_utm), !is.na(y_utm))

water <- st_read(here("data", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

land <- st_read(here("data", "boundaries", "land_small_utm", "land_small_utm.shp"), quiet = TRUE) |>
  st_make_valid()
region <- st_bbox(
  c(
    xmin = min(dat$x_utm) * 1000 - 30000, ymin = min(dat$y_utm) * 1000 - 30000,
    xmax = max(dat$x_utm) * 1000 + 30000, ymax = max(dat$y_utm) * 1000 + 30000
  ),
  crs = st_crs(32632)
)
land_region <- suppressWarnings(st_crop(land, region))

# 02 Target edge from local station spacing ----
# distinct station locations (km); repeated visits to a spot count once
stations <- dat |>
  distinct(x = x_utm, y = y_utm) |>
  as.matrix()
station_dist <- as.matrix(dist(stations))
diag(station_dist) <- Inf

spacing <- apply(station_dist, 1, \(r) sort(r, partial = spacing_k)[spacing_k])
station_edge <- pmin(pmax(spacing, edge_min), edge_max)

# target edge anywhere: the nearest station's edge, growing with distance from
# it, capped at edge_far
edge_at <- function(p) {
  d <- sqrt(outer(p[, 1], stations[, 1], "-")^2 + outer(p[, 2], stations[, 2], "-")^2)
  pmin(apply(sweep(d * growth, 2, station_edge, "+"), 1, min), edge_far)
}

# 03 Boundary and coastline ----
# outer boundary: non-convex hull 8 km around the stations
hull <- fmesher::fm_nonconvex_hull(stations, convex = 8)
boundary <- fmesher::fm_as_segm(hull)
hull_inner <- st_buffer(st_sfc(st_geometry(hull)[[1]]), -0.5)

# coastline as an interior constraint: Limfjorden shoreline (km), simplified
# to ~1 km (finer makes fmesher's refinement stall in the narrowest straits),
# clipped 0.5 km inside the hull (shoreline ends lying on the hull otherwise
# create sliver triangles); loops under 3 km (islets) are dropped
coast <- water |>
  st_union() |>
  st_simplify(dTolerance = 1000) |>
  st_boundary()
coast <- st_sfc(st_geometry(coast)[[1]] * (1 / 1000)) |>
  st_intersection(hull_inner) |>
  st_cast("LINESTRING")
coast <- coast[as.numeric(st_length(coast)) >= 3]
interior <- fmesher::fm_as_segm(coast)

# 04 Seed points ----
# the stations plus a 1 km background grid over the mesh area. Points added
# during refinement take the finest size of their neighbours, so with seeds
# only at the stations a fine bed's size spreads along long triangles far into
# sparse areas; the grid anchors the target edge everywhere. Candidates are
# thinned so that no two are closer than their target edge (every seed becomes
# a mesh vertex), keeping the finest first
grid_pts <- st_make_grid(hull_inner, cellsize = 1, what = "centers")
grid_pts <- grid_pts[lengths(st_within(grid_pts, hull_inner)) > 0]
candidates <- rbind(stations, st_coordinates(grid_pts))
candidate_edge <- c(station_edge, edge_at(st_coordinates(grid_pts)))

kept <- matrix(numeric(0), ncol = 2)
for (i in order(candidate_edge)) {
  p <- candidates[i, ]
  if (nrow(kept) == 0 || min((kept[, 1] - p[1])^2 + (kept[, 2] - p[2])^2) >= candidate_edge[i]^2) kept <- rbind(kept, p)
}

# seeds within 300 m of the shore make sliver triangles that break the barrier
# model, so they are dropped (the shoreline constraint covers them)
on_shore <- lengths(st_is_within_distance(st_as_sf(as.data.frame(kept), coords = 1:2), coast, dist = 0.3)) > 0
seeds <- kept[!on_shore, ]

# 05 Mesh ----
# quality.spec$segm lists the boundary points first, then the coastline ones
mesh <- fmesher::fm_rcdt_2d_inla(
  loc = seeds,
  cutoff = 0.2,
  boundary = boundary,
  interior = interior,
  refine = list(max.edge = Inf, max.n.strict = 8000),
  quality.spec = list(
    loc = edge_at(seeds) / refine_ratio,
    segm = edge_at(rbind(boundary$loc, interior$loc)[, 1:2]) / refine_ratio
  )
)

# 06 Barrier ----
# add_barrier_mesh() marks a triangle as land when its centre falls on the
# land map. The mesh follows the ~1 km simplified shoreline, so a few shore
# triangles holding survey stations would count as land; those triangles are
# cut out of the land map so every station lies in water
triangles <- mesh$graph$tv
triangle_sf <- map(seq_len(nrow(triangles)), \(i) st_polygon(list(mesh$loc[triangles[i, c(1:3, 1)], 1:2] * 1000))) |>
  st_sfc(crs = 32632)

station_triangles <- unique(fmesher::fm_bary(mesh, stations)$index)
land_barrier <- st_difference(land_region, st_union(triangle_sf[station_triangles]))

barrier_mesh <- sdmTMBextra::add_barrier_mesh(
  make_mesh(dat, c("x_utm", "y_utm"), mesh = mesh), land_barrier,
  range_fraction = 0.1, proj_scaling = 1000, plot = FALSE
)

stopifnot(!any(station_triangles %in% barrier_mesh$barrier_triangles))

# 07 Save ----
dir.create(here("data", "mesh"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(mesh = mesh, land_barrier = land_barrier), here("data", "mesh", "mesh.rds"))

# 08 Plot ----
mesh_sf <- st_sf(geometry = triangle_sf) |>
  mutate(barrier = seq_len(n()) %in% barrier_mesh$barrier_triangles)

ggplot() +
  geom_sf(data = mesh_sf, aes(fill = barrier), colour = "grey45", linewidth = 0.1) +
  geom_sf(data = water, fill = NA, colour = "black", linewidth = 0.35) +
  geom_point(data = dat, aes(x_utm * 1000, y_utm * 1000), size = 0.12, colour = "#B2182B") +
  scale_fill_manual(values = c(`FALSE` = "#d6e6f2", `TRUE` = "grey85"), labels = c("water", "land (barrier)"), name = NULL) +
  coord_sf(crs = 32632, xlim = range(dat$x_utm * 1000) + c(-12000, 12000), ylim = range(dat$y_utm * 1000) + c(-8000, 8000), expand = FALSE) +
  labs(subtitle = paste(mesh$n, "mesh points")) +
  theme_light()
