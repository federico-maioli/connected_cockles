# Figure 2: the model predictors, so every input to Figure 3 is shown before
# the results - the five continuous environmental covariates and sediment
# type at their native raster resolution (DTU Aqua rasters, reprojected to UTM
# 32N at 100 m), fishing effort of the cockle fishery (swept-area ratio, 250 m
# grid; the models use the season before each survey, the map shows the mean
# over the nine seasons 2016-17 to 2024-25), and presence-weighted in-strength
# on the 2 km connectivity grid. Every panel is clipped to the 2 km grid cells
# holding a survey station, extent the sampled area with a 3 km margin. Each
# environmental covariate has its own muted ramp (dark = "more"); in-strength,
# the only connectivity predictor, uses the multi-hue viridis palette.
# 6.2 x 6.4 in (fits a double-column page) at 600 dpi, 7-8 pt text.

library(tidyverse)
library(here)
library(sf)
library(patchwork)
library(terra)

source(here("R", "helpers.R")) # sediment_folk7(), sar_seasons()

# 01 Data ----
dat <- readRDS(here("data", "cockles", "derived", "cockles_env.rds")) |>
  filter(!is.na(x_utm), !is.na(y_utm))
water <- st_read(here("data", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

# only the cells with at least one survey sample, as in Figure 1a (same cell
# id formula as 01_prepare_grid.R); in-strength is 0 outside the
# presence-weighted network, so it is set to NA and left blank rather than
# drawn as a real 0
cell_id <- function(x, y) (floor((x - 450074) / 2000) + 1) + floor((y - 6258093) / 2000) * 65
sampled_ids <- unique(cell_id(dat$x_utm * 1000, dat$y_utm * 1000))

grid <- readRDS(here("data", "grid", "grid_env_conn.rds")) |>
  filter(id %in% sampled_ids) |>
  mutate(in_strength = na_if(conn_in_strength_presence, 0))

# extent: the sampled area with a 3 km margin, so the maps fill the panels
xlim <- range(dat$x_utm * 1000) + c(-3000, 3000)
ylim <- range(dat$y_utm * 1000) + c(-3000, 3000)

# white mask = the plot frame minus the fjord, drawn over the tiles so the 2 km
# cells are clipped to the coastline instead of spilling onto land
frame <- st_as_sfc(st_bbox(
  c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
  crs = st_crs(32632)
))
land_mask <- st_difference(frame, st_union(water))

# raster cells to a data frame, clipped to the sampled grid cells (so the
# south-eastern basins, outside the surveyed area, are left blank)
fig_extent <- ext(xlim[1], xlim[2], ylim[1], ylim[2])
raster_df <- function(r) {
  as.data.frame(crop(r, fig_extent), xy = TRUE, na.rm = TRUE) |>
    setNames(c("x", "y", "value")) |>
    filter(cell_id(x, y) %in% sampled_ids)
}

# the environmental rasters, as extracted in 02_clean_cockles.R, reprojected to
# UTM 32N at 100 m with nearest neighbour (values stay raw; depth is 50 m and
# shear stress 200 m natively). Depth is negative below sea level in the
# raster, flipped so larger means deeper; sediment codes to Folk 7 classes
env_files <- c(
  depth = "depth.asc",
  temp = "temp_bot_mean.asc",
  sal = "salt_bot_mean.asc",
  oxy = "do4mgl_bot_mean.asc",
  shear_max = "tw_bot_mean_of_max.asc",
  sediment = "sediment.asc"
)
env_rasters <- imap(env_files, \(file, nm) {
  r <- rast(here("data", "env", "Asc4dtuaqua", file)) |>
    project("EPSG:32632", res = 100, method = "near")
  if (nm == "depth") r <- -r
  raster_df(r)
})
env_rasters$sediment <- mutate(env_rasters$sediment, value = sediment_folk7(value)) |>
  filter(!is.na(value))

# fishing effort: mean swept-area ratio over the nine seasons
fishing <- raster_df(mean(sar_seasons()))

# 02 Panels ----
titles <- list(
  depth = expression(Depth ~ (m)),
  temp = expression(Temperature ~ (degree * C)),
  sal = expression(Salinity ~ (psu)),
  oxy = expression("Oxygen (days < 4 mg O"[2] * " l"^{-1} * ")"),
  shear_max = expression("Max. shear stress (N m"^{-2} * ")"),
  sediment = "Sediment (Folk 7)",
  fishing = "Fishing (swept-area ratio)",
  in_strength = "In-strength"
)
# each environmental covariate has its own muted single-hue ramp (none teal or
# green); in-strength uses the vivid multi-hue viridis palette, so the
# connectivity metric stands apart from the environment at a glance
scales <- list(
  depth = scale_fill_gradient(low = "#e6eaf0", high = "#34496a", na.value = NA, name = NULL, breaks = c(0, 10, 20)),
  # temperature and salinity vary little across most of the fjord, so their
  # ramps are capped at the 2nd and 98th percentiles (values beyond take the end
  # colours)
  temp = scale_fill_gradientn(
    colours = c("#fff5f0", "#fb6a4a", "#67000d"), na.value = NA, name = NULL,
    limits = quantile(env_rasters$temp$value, c(0.02, 0.98)), oob = scales::squish,
    breaks = scales::breaks_pretty(n = 3)
  ),
  sal = scale_fill_gradientn(
    colours = c("#f7f4e6", "#c2a83e", "#4a3f0f"), na.value = NA, name = NULL,
    limits = quantile(env_rasters$sal$value, c(0.02, 0.98)), oob = scales::squish,
    breaks = scales::breaks_pretty(n = 3)
  ),
  oxy = scale_fill_gradient(low = "#efefef", high = "#303030", na.value = NA, name = NULL, transform = "pseudo_log", breaks = c(0, 10, 100)),
  shear_max = scale_fill_gradient(low = "#ede6f0", high = "#5a3f6e", na.value = NA, name = NULL, breaks = c(1, 2, 3, 4)),
  sediment = scale_fill_manual(
    values = c(mud = "#7a6650", muddy_sand = "#b29d78", sand = "#e6d396", mixed_sediment = "#b9b9b9"),
    labels = c(mud = "Mud", muddy_sand = "Muddy sand", sand = "Sand", mixed_sediment = "Mixed"),
    name = NULL
  ),
  # most cells are never fished: square-root scale so fished beds stand out
  fishing = scale_fill_gradient(low = "#f4f1ea", high = "#7a3b12", na.value = NA, name = NULL, transform = "sqrt", breaks = c(0, 0.5, 2)),
  in_strength = scale_fill_viridis_c(
    option = "viridis", direction = -1, transform = "sqrt", na.value = NA,
    breaks = c(0.001, 0.01), labels = c("0.001", "0.01"), name = NULL
  )
)

# horizontal legends under each map: a colour bar for continuous predictors,
# keys for sediment
legend_guide <- function(col) {
  if (col == "sediment") {
    guide_legend(nrow = 1, theme = theme(
      legend.key.size = unit(2, "mm"), legend.key.spacing.x = unit(2, "mm"),
      legend.text = element_text(size = 5.5, margin = margin(l = 1))
    ))
  } else {
    guide_colourbar(theme = theme(legend.key.width = unit(16, "mm"), legend.key.height = unit(1.4, "mm")))
  }
}

panels <- map(names(titles), function(col) {
  layer <- switch(col,
    in_strength = geom_tile(data = grid, aes(x, y, fill = in_strength), width = 2000, height = 2000),
    fishing = geom_raster(data = fishing, aes(x, y, fill = value)),
    geom_raster(data = env_rasters[[col]], aes(x, y, fill = value))
  )
  ggplot() +
    layer +
    geom_sf(data = land_mask, fill = "white", colour = NA) +
    geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.15) +
    scales[[col]] +
    guides(fill = legend_guide(col)) +
    coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
    labs(title = titles[[col]]) +
    theme_void(base_size = 7) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = 7, hjust = 0.5, colour = "grey15", margin = margin(b = 1)),
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "mm"),
      legend.text = element_text(size = 5.5, colour = "grey30", margin = margin(t = 0.5)),
      plot.margin = margin(3, 4, 3, 4)
    )
})

# 03 Combine and save ----
fig2 <- wrap_plots(panels, ncol = 3) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 7.5, face = "bold", colour = "grey20"))

dir.create(here("output", "figs", "main"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "main", "fig2_predictors.png"), fig2, width = 6.2, height = 6.4, dpi = 600, bg = "white")
