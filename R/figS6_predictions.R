# Figure S6: cockle biomass predicted by the selected model (Space +
# Environment), on the 2 x 2 km grid cells holding at least one survey
# station: (a) expected biomass = presence probability x biomass where
# present, for the grab survey (Stock) in the most recent year (2025),
# predictions below 1 g/m2 shown as 0; (b) its uncertainty, the SD of 500
# draws from the joint precision matrix, on the log scale.
#
# Reads data/sdm/main/space_env.rds (07_fit_sdm.R) and
# data/grid/grid_env_conn.rds (06_match_connectivity.R).

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
biomass_cols <- carto_pal(7, "BurgYl")

# 01 Data ----
fit <- readRDS(here("data", "sdm", "main", "space_env.rds"))
water <- st_read(here("data", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

xlim <- range(fit$data$x_utm * 1000) + c(-16000, 16000)
ylim <- range(fit$data$y_utm * 1000) + c(-6000, 6000)

# white mask = the plot frame minus the fjord, drawn over the tiles so the 2 km
# cells are clipped to the coastline instead of spilling onto land
frame <- st_as_sfc(st_bbox(
  c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
  crs = st_crs(32632)
))
land_mask <- st_difference(frame, st_union(water))

# 02 Predict ----
# at the centre of every grid cell holding at least one survey station; grid
# covariates are standardised with the model data's own mean and SD, recovered
# from the raw and standardised columns of the fit
sampled_ids <- fit$data |>
  transmute(id = (floor((x_utm * 1000 - 450074) / 2000) + 1) + floor((y_utm * 1000 - 6258093) / 2000) * 65) |>
  distinct(id) |>
  pull(id)

std_vars <- c("depth", "temp", "sal", "oxy", "shear_max")
newdata <- readRDS(here("data", "grid", "grid_env_conn.rds")) |>
  filter(id %in% sampled_ids) |>
  mutate(
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

pred <- predict(fit, newdata = newdata, type = "response")

# uncertainty: 500 draws from the joint precision matrix; for a delta model
# these are the combined prediction on the log scale, log(presence probability
# x biomass where present). The SD across draws is taken on that log scale:
# back-transformed draws are heavy-tailed, and a few extreme ones dominate an
# SD in g/m2
sims <- predict(fit, newdata = newdata, nsim = 500)
pred$sd <- apply(sims, 1, sd)
pred <- mutate(pred, est = if_else(est < biomass_zero, 0, est))
summary(pred$est)
summary(pred$sd)

# 03 Plot ----
map_layers <- list(
  geom_sf(data = land_mask, fill = "white", colour = NA),
  geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.3),
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE),
  theme_void(base_size = 11),
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "bottom",
    legend.box.just = "bottom",
    legend.spacing.x = unit(2, "mm"),
    legend.title = element_text(size = 11, colour = "grey15"),
    legend.text = element_text(size = 10, colour = "grey30")
  )
)

p_est <- ggplot() +
  # zero cells get their own legend key through a dummy linetype mapping (the
  # tile outline itself is not drawn, linewidth = 0)
  geom_tile(data = filter(pred, est == 0), aes(x, y, linetype = "0"), fill = "grey85", linewidth = 0, width = 2000, height = 2000) +
  geom_tile(data = filter(pred, est > 0), aes(x, y, fill = est), width = 2000, height = 2000) +
  map_layers +
  scale_fill_gradientn(
    colours = biomass_cols, trans = "pseudo_log", name = biomass_lab,
    breaks = biomass_breaks, labels = scales::comma(biomass_breaks),
    limits = c(biomass_zero, NA) # bar starts at the zero cut-off, so its first label is 1
  ) +
  scale_linetype_manual(values = "solid", name = NULL) +
  annotation_scale(location = "bl", width_hint = 0.2, height = unit(0.15, "cm"), text_cex = 0.7, line_col = "grey40", text_col = "grey40") +
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
      theme = theme(legend.key.height = unit(3.5, "mm"), legend.key.width = unit(45, "mm"))
    )
  )

p_sd <- ggplot() +
  geom_tile(data = pred, aes(x, y, fill = sd), width = 2000, height = 2000) +
  map_layers +
  scale_fill_gradientn(
    colours = carto_pal(7, "Purp"), name = "SD (log scale)"
  ) +
  annotation_north_arrow(location = "br", which_north = "true", height = unit(0.9, "cm"), width = unit(0.7, "cm"), style = north_arrow_minimal(line_col = "grey40", text_col = "grey40", fill = "grey40")) +
  guides(fill = guide_colourbar(
    title.position = "top", title.hjust = 0.5,
    theme = theme(legend.key.height = unit(3.5, "mm"), legend.key.width = unit(45, "mm"))
  ))

fig_s8 <- p_est + p_sd +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 12, face = "bold", colour = "grey20"))

dir.create(here("output", "figs", "supp"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "supp", "figS6_predictions.png"), fig_s8, width = 13, height = 6.5, dpi = 600, bg = "white")
