# Covariate maps: each model covariate at the survey points over the Limfjorden
# outline, styled like Figure 1 (white background, water outline only).

library(tidyverse)
library(here)
library(sf)
library(patchwork)

# 01 Data ----
dat <- readRDS(here("data", "derived", "cockles_connectivity.rds")) |>
  filter(!is.na(x_utm), !is.na(y_utm))
water <- st_read(here("data", "raw", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

pts <- dat |>
  mutate(x_m = x_utm * 1000, y_m = y_utm * 1000) |>
  st_as_sf(coords = c("x_m", "y_m"), crs = 32632)

xlim <- range(dat$x_utm * 1000) + c(-16000, 16000)
ylim <- range(dat$y_utm * 1000) + c(-6000, 6000)

# 02 Covariates to map ----
# each covariate keeps its own colour scheme; dark = "more" (deeper, warmer, more
# saline, more hypoxic days, higher shear). Units in superscript notation.
covars <- c("depth", "temp", "sal", "oxy", "shear_max")
titles <- list(
  depth = expression(Depth ~ (m)),
  temp = expression(Temperature ~ (degree * C)),
  sal = expression(Salinity ~ (psu)),
  oxy = expression("Oxygen deficiency (" * "days yr"^{-1} * " < 4 mg O"[2] * " l"^{-1} * ")"),
  shear_max = expression("Max. shear stress (" * N ~ m^{-2} * ")")
)
trans <- c(depth = "identity", temp = "identity", sal = "identity", oxy = "pseudo_log", shear_max = "identity")

covar_scale <- function(col, tr) {
  switch(col,
    depth = scale_colour_gradient(low = "#d0e1f2", high = "#08306b", name = NULL, trans = tr),
    temp = scale_colour_gradient(low = "#fee0d2", high = "#67000d", name = NULL, trans = tr),
    sal = scale_colour_gradient(low = "#e5f5e0", high = "#005a32", name = NULL, trans = tr),
    oxy = scale_colour_gradient(low = "#efefef", high = "black", name = NULL, trans = tr),
    shear_max = scale_colour_viridis_c(option = "rocket", direction = -1, name = NULL, trans = tr)
  )
}

# 03 Panel builder ----
covar_panel <- function(col) {
  d <- pts |> arrange(.data[[col]])
  ggplot() +
    geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.22) +
    geom_sf(data = d, aes(colour = .data[[col]]), size = 0.45, alpha = 0.85) +
    covar_scale(col, trans[[col]]) +
    coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
    labs(title = titles[[col]]) +
    theme_void(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = 9.5, hjust = 0.5, colour = "grey20"),
      legend.position = "right",
      legend.key.width = unit(3, "mm"), legend.key.height = unit(8, "mm"),
      legend.text = element_text(size = 7, colour = "grey30"),
      plot.tag = element_text(size = 12, face = "bold", colour = "grey20"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

panels <- map(covars, covar_panel)

# 04 Combine and save ----
# 3 columns (3 on top, 2 below); tag panels a-e
fig <- wrap_plots(panels, ncol = 3) +
  plot_annotation(tag_levels = "a")

dir.create(here("output"), showWarnings = FALSE)
ggsave(here("output", "fig_covariates.png"), fig, width = 13, height = 6.5, dpi = 600, bg = "white")
