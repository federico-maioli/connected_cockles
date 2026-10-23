# Figure S7: predictions of the full two-part model (Space + Environment +
# Connectivity, presence-weighted in-strength) for one grab sample (Stock) in
# the most recent year (2025), at the centre of the 2 x 2 km grid cells:
# (a) occurrence probability, every cell holding a survey station, and (b)
# its uncertainty, the SD of the probability across 500 draws; (c) expected
# biomass = occurrence probability x biomass where present, only in cells
# holding grab stations (the biomass part is fitted to grab beds only), for
# beds not fished in the previous season, predictions below 1 g/m2 shown as
# 0, and (d) its uncertainty, the SD of 500 draws from the joint precision
# matrix of each part, on the log scale.
#
# Reads data/sdm/main/space_env_conn.rds (07_fit_sdm.R: list(occurrence,
# biomass)) and data/grid/grid_env_conn.rds (06_match_connectivity.R).

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
library(rcartocolor)
library(patchwork)

biomass_lab <- expression(Predicted ~ biomass ~ (g ~ m^-2))
biomass_breaks <- c(1, 10, 100, 1000)
# predictions below this (g/m2) are shown as 0
biomass_zero <- 1
biomass_cols <- carto_pal(7, "BurgYl")

# 01 Data ----
fits <- readRDS(here("data", "sdm", "main", "space_env_conn.rds"))
dat <- fits$occurrence$data
water <- st_read(here("data", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)

# extent: the sampled area with a 3 km margin, as in Figure 2
xlim <- range(dat$x_utm * 1000) + c(-3000, 3000)
ylim <- range(dat$y_utm * 1000) + c(-3000, 3000)

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
cell_id <- function(x_utm, y_utm) (floor((x_utm * 1000 - 450074) / 2000) + 1) + floor((y_utm * 1000 - 6258093) / 2000) * 65
sampled_ids <- unique(cell_id(dat$x_utm, dat$y_utm))
grab_ids <- with(filter(dat, survey == "Stock"), unique(cell_id(x_utm, y_utm)))

std_vars <- c("depth", "temp", "sal", "oxy", "shear_max", "conn_in_strength")
newdata <- readRDS(here("data", "grid", "grid_env_conn.rds")) |>
  filter(id %in% sampled_ids) |>
  mutate(
    conn_in_strength = conn_in_strength_presence,
    x_utm = x / 1000, y_utm = y / 1000,
    survey = factor("Stock", levels = levels(dat$survey)),
    year = factor("2025", levels = levels(dat$year))
  )
for (v in std_vars) {
  s <- sd(dat[[v]]) / sd(dat[[paste0(v, "_std")]])
  m <- mean(dat[[v]]) - s * mean(dat[[paste0(v, "_std")]])
  newdata[[paste0(v, "_std")]] <- (newdata[[v]] - m) / s
}
newdata <- filter(newdata, if_all(all_of(c(paste0(std_vars, "_std"), "sediment")), \(x) !is.na(x)))
# biomass where present for a bed not fished in the previous season: the
# standardised log(1 + SAR) of 0, recovered from the raw and standardised
# columns of the fit
s_sar <- sd(log1p(dat$sar_prev)) / sd(dat$sar_prev_std)
m_sar <- mean(log1p(dat$sar_prev)) - s_sar * mean(dat$sar_prev_std)
newdata_grab <- newdata |>
  filter(id %in% grab_ids) |>
  mutate(
    year = factor("2025", levels = levels(fits$biomass$data$year)),
    sar_prev_std = (0 - m_sar) / s_sar
  )

# occurrence probability of one grab sample (offset 0, as for grab stations
# in the fit): inverse cloglog of the prediction
inv_cloglog <- function(x) 1 - exp(-exp(x))
pred_occ <- predict(fits$occurrence, newdata = newdata, offset = rep(0, nrow(newdata))) |>
  mutate(p = inv_cloglog(est))
# its uncertainty: SD of the occurrence probability across 500 draws
set.seed(1)
sims_occ_all <- predict(fits$occurrence, newdata = newdata, offset = rep(0, nrow(newdata)), nsim = 500)
pred_occ$sd <- apply(inv_cloglog(sims_occ_all), 1, sd)

# expected biomass in the grab cells: occurrence probability x biomass where
# present. Uncertainty from 500 draws of each part (the parts share no
# parameters, so their draws are independent), combined on the log scale,
# log(p) + log(mu); the SD is taken on that scale, as back-transformed draws
# are heavy-tailed
set.seed(1)
occ_grab <- predict(fits$occurrence, newdata = newdata_grab, offset = rep(0, nrow(newdata_grab)))
bio_grab <- predict(fits$biomass, newdata = newdata_grab)
sims_occ <- predict(fits$occurrence, newdata = newdata_grab, offset = rep(0, nrow(newdata_grab)), nsim = 500)
sims_bio <- predict(fits$biomass, newdata = newdata_grab, nsim = 500)
pred <- newdata_grab |>
  select(x, y) |>
  mutate(
    est = inv_cloglog(occ_grab$est) * exp(bio_grab$est),
    sd = apply(log(inv_cloglog(sims_occ)) + sims_bio, 1, sd),
    est = if_else(est < biomass_zero, 0, est)
  )
summary(pred_occ$p)
summary(pred$est)
summary(pred$sd)

# 03 Plot ----
# as Figure 2: panel titles above the maps, compact horizontal legends under
# them; 7 in wide (double column) at 600 dpi, 7-8 pt text
colourbar <- guide_colourbar(theme = theme(legend.key.width = unit(16, "mm"), legend.key.height = unit(1.4, "mm")))
map_layers <- list(
  geom_sf(data = land_mask, fill = "white", colour = NA),
  geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.15),
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE),
  theme_void(base_size = 7),
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 7.5, hjust = 0.5, colour = "grey15", margin = margin(b = 2)),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    legend.spacing.x = unit(1, "mm"),
    legend.text = element_text(size = 6, colour = "grey30", margin = margin(t = 0.5)),
    plot.margin = margin(3, 4, 3, 4)
  )
)

p_occ <- ggplot() +
  geom_tile(data = pred_occ, aes(x, y, fill = p), width = 2000, height = 2000) +
  map_layers +
  # square-root scale: probabilities for a 0.1 m2 grab are mostly very low
  scale_fill_gradientn(
    colours = carto_pal(7, "Teal"), name = NULL,
    limits = c(0, NA), transform = "sqrt", breaks = c(0, 0.1, 0.4)
  ) +
  guides(fill = colourbar) +
  labs(title = "Occurrence probability")

p_est <- ggplot() +
  # zero cells get their own legend key through a dummy linetype mapping (the
  # tile outline itself is not drawn, linewidth = 0)
  geom_tile(data = filter(pred, est == 0), aes(x, y, linetype = "0"), fill = "grey85", linewidth = 0, width = 2000, height = 2000) +
  geom_tile(data = filter(pred, est > 0), aes(x, y, fill = est), width = 2000, height = 2000) +
  map_layers +
  scale_fill_gradientn(
    colours = biomass_cols, transform = "pseudo_log", name = NULL,
    breaks = biomass_breaks, labels = scales::comma(biomass_breaks),
    limits = c(biomass_zero, NA) # bar starts at the zero cut-off, so its first label is 1
  ) +
  scale_linetype_manual(values = "solid", name = NULL) +
  # the "0" key sits just left of the colour bar, same height, so it reads as
  # the bar's first step
  guides(
    linetype = guide_legend(
      order = 1, label.position = "bottom",
      override.aes = list(fill = "grey85", linewidth = 0),
      theme = theme(legend.key.height = unit(1.4, "mm"), legend.key.width = unit(3, "mm"))
    ),
    fill = guide_colourbar(order = 2, theme = theme(legend.key.width = unit(16, "mm"), legend.key.height = unit(1.4, "mm")))
  ) +
  labs(title = expression(Expected ~ biomass ~ (g ~ m^-2)))

p_sd <- ggplot() +
  geom_tile(data = pred, aes(x, y, fill = sd), width = 2000, height = 2000) +
  map_layers +
  scale_fill_gradientn(colours = carto_pal(7, "Purp"), name = NULL, breaks = scales::breaks_pretty(n = 3)) +
  guides(fill = colourbar) +
  labs(title = "SD of log biomass")

p_occ_sd <- ggplot() +
  geom_tile(data = pred_occ, aes(x, y, fill = sd), width = 2000, height = 2000) +
  map_layers +
  # square-root scale, as for the probability itself
  scale_fill_gradientn(colours = carto_pal(7, "Purp"), name = NULL, transform = "sqrt", breaks = c(0, 0.05, 0.2)) +
  guides(fill = colourbar) +
  labs(title = "SD of occurrence probability")

# occurrence and its uncertainty on top, biomass and its uncertainty below
fig_s7 <- p_occ + p_occ_sd + p_est + p_sd +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 7.5, face = "bold", colour = "grey20"))

dir.create(here("output", "figs", "supp"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "supp", "figS7_predictions.png"), fig_s7, width = 4.9, height = 4.6, dpi = 600, bg = "white")
