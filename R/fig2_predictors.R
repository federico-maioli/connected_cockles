# Figure 2: the model predictors - the five environmental covariates at their
# native raster resolution (DTU Aqua rasters, reprojected to UTM 32N at 100 m)
# and presence-weighted in-strength on the 2 km connectivity grid (the
# connectivity metric used in the models), so every input to Figure 3 is shown
# before the results. Same
# extent and coastline mask as Figure 1. Each environmental covariate has its
# own muted ramp (dark = "more"); in-strength, the only non-environmental
# predictor, uses the vivid multi-hue viridis palette.

library(tidyverse)
library(here)
library(sf)
library(patchwork)
library(terra)

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
sampled_ids <- dat |>
  mutate(id = (floor((x_utm * 1000 - 450074) / 2000) + 1) + floor((y_utm * 1000 - 6258093) / 2000) * 65) |>
  distinct(id) |>
  pull(id)

grid <- readRDS(here("data", "grid", "grid_env_conn.rds")) |>
  filter(id %in% sampled_ids) |>
  mutate(in_strength = na_if(conn_in_strength_presence, 0))

# same extent as Figure 1
xlim <- range(dat$x_utm * 1000) + c(-16000, 16000)
ylim <- range(dat$y_utm * 1000) + c(-6000, 6000)

# white mask = the plot frame minus the fjord, drawn over the tiles so the 2 km
# cells are clipped to the coastline instead of spilling onto land
frame <- st_as_sfc(st_bbox(
  c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
  crs = st_crs(32632)
))
land_mask <- st_difference(frame, st_union(water))

# the environmental rasters, as extracted in 02_clean_cockles.R, reprojected to
# UTM 32N at 100 m with nearest neighbour (values stay raw; depth is 50 m and
# shear stress 200 m natively) and cropped to the figure extent. Depth is
# negative below sea level in the raster, flipped so larger means deeper
env_files <- c(
  depth = "depth.asc",
  temp = "temp_bot_mean.asc",
  sal = "salt_bot_mean.asc",
  oxy = "do4mgl_bot_mean.asc",
  shear_max = "tw_bot_mean_of_max.asc"
)
fig_extent <- ext(xlim[1], xlim[2], ylim[1], ylim[2])
env_rasters <- imap(env_files, \(file, nm) {
  r <- rast(here("data", "env", "Asc4dtuaqua", file)) |>
    project("EPSG:32632", res = 100, method = "near") |>
    crop(fig_extent)
  if (nm == "depth") r <- -r
  as.data.frame(r, xy = TRUE, na.rm = TRUE) |>
    setNames(c("x", "y", "value"))
})

# 02 Panels ----
titles <- list(
  depth = expression(Depth ~ (m)),
  temp = expression(Temperature ~ (degree * C)),
  sal = expression(Salinity ~ (psu)),
  oxy = expression("Oxygen (" * "days yr"^{-1} * " < 4 mg O"[2] * " l"^{-1} * ")"),
  shear_max = expression("Max. shear stress (" * N ~ m^{-2} * ")"),
  in_strength = "In-strength"
)
# each environmental covariate has its own muted single-hue ramp (none teal or
# green); in-strength uses the vivid multi-hue viridis palette, so the
# connectivity metric stands apart from the environment at a glance
scales <- list(
  depth = scale_fill_gradient(low = "#e6eaf0", high = "#34496a", na.value = NA, name = NULL),
  # temperature and salinity vary little across most of the fjord, so their
  # ramps are capped at the 2nd and 98th percentiles (values beyond take the end
  # colours)
  temp = scale_fill_gradientn(
    colours = c("#fff5f0", "#fb6a4a", "#67000d"), na.value = NA, name = NULL,
    limits = quantile(env_rasters$temp$value, c(0.02, 0.98)), oob = scales::squish
  ),
  sal = scale_fill_gradientn(
    colours = c("#f7f4e6", "#c2a83e", "#4a3f0f"), na.value = NA, name = NULL,
    limits = quantile(env_rasters$sal$value, c(0.02, 0.98)), oob = scales::squish
  ),
  oxy = scale_fill_gradient(low = "#efefef", high = "#303030", na.value = NA, name = NULL, transform = "pseudo_log"),
  shear_max = scale_fill_gradient(low = "#ede6f0", high = "#5a3f6e", na.value = NA, name = NULL),
  in_strength = scale_fill_viridis_c(
    option = "viridis", direction = -1, transform = "sqrt", na.value = NA,
    breaks = c(0.001, 0.005, 0.01, 0.02), labels = c("0.001", "0.005", "0.01", "0.02"), name = NULL
  )
)

panels <- map(names(titles), function(col) {
  layer <- if (col == "in_strength") {
    geom_tile(data = grid, aes(x, y, fill = in_strength), width = 2000, height = 2000)
  } else {
    geom_raster(data = env_rasters[[col]], aes(x, y, fill = value))
  }
  ggplot() +
    layer +
    geom_sf(data = land_mask, fill = "white", colour = NA) +
    geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.22) +
    scales[[col]] +
    coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
    labs(title = titles[[col]]) +
    theme_void(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = 10, hjust = 0.5, colour = "grey20"),
      legend.position = "right",
      legend.key.width = unit(3, "mm"), legend.key.height = unit(10, "mm"),
      legend.text = element_text(size = 8, colour = "grey30"),
      plot.margin = margin(4, 4, 4, 4)
    )
})

# 03 Combine and save ----
fig2 <- wrap_plots(panels, ncol = 3) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 12, face = "bold", colour = "grey20"))

dir.create(here("output", "figs", "main"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "main", "fig2_predictors.png"), fig2, width = 12, height = 6.6, dpi = 600, bg = "white")
