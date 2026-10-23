# Figure 3: the full two-part model (space + environment + in-strength
# connectivity), its two parts side by side:
#   occurrence                binomial (cloglog) - where cockles occur
#   biomass where present     generalised gamma (log) - how much, where they
#                                                       occur
# Panel a: standardized coefficient of every predictor in each part.
# Panel b: partial effects - predicted occurrence probability (one grab
# sample) and biomass where present along each predictor, with the others at
# their mean; sediment type affects occurrence only, previous-season fishing
# biomass only.
#
# Reads data/sdm/main/space_env_conn.rds (07_fit_sdm.R): list(occurrence,
# biomass).

library(tidyverse)
library(here)
library(sdmTMB)
library(ggstats)
library(patchwork)

part_labs <- c(presence = "Occurrence (cloglog link)", biomass = "Biomass where present (log link)")
part_point_cols <- set_names(c("#5D3A9B", "#E66100"), part_labs)

term_labs <- c(
  depth_std = "Depth", `I(depth_std^2)` = "Depth\u00b2", temp_std = "Temperature", oxy_std = "Oxygen",
  sal_std = "Salinity", shear_max_std = "Shear stress",
  sedimentmuddy_sand = "Muddy sand (vs mud)", sedimentsand = "Sand (vs mud)", sedimentmixed_sediment = "Mixed sediment (vs mud)",
  conn_in_strength_std = "In-strength",
  sar_prev_std = "Fishing"
)
term_levels <- c(
  "Fishing", "In-strength", "Mixed sediment (vs mud)", "Sand (vs mud)", "Muddy sand (vs mud)", "Shear stress", "Salinity", "Oxygen", "Temperature",
  "Depth\u00b2", "Depth"
)

fits <- readRDS(here("data", "sdm", "main", "space_env_conn.rds"))

# 01 Panel a: coefficients of both parts ----
# sediment enters occurrence only, previous-season fishing biomass only
coefs <- map(1:2, \(m) tidy(fits[[m]], effects = "fixed", conf.int = TRUE) |> mutate(component = m)) |>
  list_rbind() |>
  filter(term %in% names(term_labs)) |>
  mutate(
    term = factor(term_labs[term], levels = term_levels),
    part = factor(unname(part_labs)[component], levels = unname(part_labs))
  )

p_coeff <- ggplot(coefs, aes(estimate, term, colour = part)) +
  geom_stripped_rows(colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0, linewidth = 0.9, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = part_point_cols, name = NULL) +
  labs(x = "Standardized coefficient", y = NULL) +
  theme_light(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# 02 Panel b: partial effects ----
# predicted occurrence probability (left axis) and log biomass where present
# (right axis, log scale) along each predictor, from the population-level part
# of each fit (re_form = NA, no spatial field), with all other predictors at
# their mean (0, standardised), on mud (the reference sediment), for one grab
# sample in 2021, the year with the most grab stations. Occurrence spans the
# 1st-99th percentile of all stations, biomass that of the grab stations with
# cockles (no extrapolation beyond the data each part is fitted to); depth
# includes its quadratic term. Ribbons are 95% CIs
pe_year <- "2021"
pe_vars <- c(
  depth = "Depth (m)", temp = "Temperature (°C)", sal = "Salinity (psu)",
  oxy = "Oxygen (days < 4 mg/l)", shear_max = "Max. shear stress (N/m²)",
  conn_in_strength = "In-strength"
)
std_cols <- paste0(names(pe_vars), "_std")
dat_occ <- fits$occurrence$data
dat_bio <- fits$biomass$data

# newdata along one predictor over a part's data range; covariates are
# standardised with the scaling of all stations (07_fit_sdm.R)
pe_newdata <- function(var, data) {
  s <- sd(dat_occ[[var]]) / sd(dat_occ[[paste0(var, "_std")]])
  m <- mean(dat_occ[[var]]) - s * mean(dat_occ[[paste0(var, "_std")]])
  nd <- tibble(value = seq(quantile(data[[var]], 0.01), quantile(data[[var]], 0.99), length.out = 100))
  for (col in std_cols) nd[[col]] <- 0
  nd[[paste0(var, "_std")]] <- (nd$value - m) / s
  mutate(
    nd,
    survey = factor("Stock", levels = levels(dat_occ$survey)),
    sediment = factor("mud", levels = levels(dat_occ$sediment)),
    sar_prev_std = 0, # mean fishing effort (biomass part only)
    x_utm = mean(dat_occ$x_utm), y_utm = mean(dat_occ$y_utm)
  )
}

inv_cloglog <- function(x) 1 - exp(-exp(x))

partial_effect <- function(var) {
  nd_occ <- mutate(pe_newdata(var, dat_occ), year = factor(pe_year, levels = levels(dat_occ$year)))
  nd_bio <- mutate(pe_newdata(var, dat_bio), year = factor(pe_year, levels = levels(dat_bio$year)))
  # a grab sample: offset 0, as for every grab station in the fit
  presence <- predict(fits$occurrence, newdata = nd_occ, re_form = NA, se_fit = TRUE, offset = rep(0, nrow(nd_occ)))
  biomass <- predict(fits$biomass, newdata = nd_bio, re_form = NA, se_fit = TRUE)
  bind_rows(
    tibble(
      part = "presence", value = nd_occ$value,
      est = inv_cloglog(presence$est), low = inv_cloglog(presence$est - 1.96 * presence$est_se),
      high = inv_cloglog(presence$est + 1.96 * presence$est_se)
    ),
    tibble(
      part = "biomass", value = nd_bio$value,
      est = biomass$est, low = biomass$est - 1.96 * biomass$est_se, high = biomass$est + 1.96 * biomass$est_se
    )
  )
}

# one plot per predictor: log biomass (the link scale) is rescaled onto the
# probability axis and labelled on a secondary axis
pe_plot <- function(var) {
  d <- partial_effect(var)
  d_p <- filter(d, part == "presence")
  d_b <- filter(d, part == "biomass")
  p_max <- max(d_p$high)
  b_range <- range(c(d_b$low, d_b$high))
  to_p <- function(b) (b - b_range[1]) / diff(b_range) * p_max
  ggplot() +
    geom_ribbon(data = d_p, aes(value, ymin = low, ymax = high), fill = part_point_cols[[1]], alpha = 0.2) +
    geom_line(data = d_p, aes(value, est), colour = part_point_cols[[1]], linewidth = 0.8) +
    geom_ribbon(data = d_b, aes(value, ymin = to_p(low), ymax = to_p(high)), fill = part_point_cols[[2]], alpha = 0.2) +
    geom_line(data = d_b, aes(value, to_p(est)), colour = part_point_cols[[2]], linewidth = 0.8) +
    # axis titles are drawn once for the whole panel b (see 03 below)
    scale_y_continuous(
      name = NULL,
      sec.axis = sec_axis(\(y) b_range[1] + y / p_max * diff(b_range), name = NULL)
    ) +
    scale_x_continuous(expand = c(0, 0), breaks = scales::breaks_pretty(n = 4)) +
    labs(x = pe_vars[[var]]) +
    theme_light(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title.y.left = element_text(colour = part_point_cols[[1]]),
      axis.title.y.right = element_text(colour = part_point_cols[[2]]),
      axis.text.y.left = element_text(colour = part_point_cols[[1]]),
      axis.text.y.right = element_text(colour = part_point_cols[[2]])
    )
}

# the y-axis titles are identical in every plot, so patchwork shows each once;
# the panel tag sits on the first plot
# previous-season fishing enters biomass only: one curve, on the right
# (biomass) axis like the other panels, along the raw swept-area ratio from 0
# to its 99th percentile at the grab stations with cockles (standardised as
# log(1 + SAR), 07_fit_sdm.R)
pe_fishing <- function() {
  s <- sd(log1p(dat_occ$sar_prev)) / sd(dat_occ$sar_prev_std)
  m <- mean(log1p(dat_occ$sar_prev)) - s * mean(dat_occ$sar_prev_std)
  nd <- tibble(value = seq(0, quantile(dat_bio$sar_prev, 0.99), length.out = 100))
  for (col in std_cols) nd[[col]] <- 0
  nd <- mutate(
    nd,
    sar_prev_std = (log1p(value) - m) / s,
    year = factor(pe_year, levels = levels(dat_bio$year)),
    x_utm = mean(dat_occ$x_utm), y_utm = mean(dat_occ$y_utm)
  )
  biomass <- predict(fits$biomass, newdata = nd, re_form = NA, se_fit = TRUE)
  d <- tibble(value = nd$value, est = biomass$est, low = est - 1.96 * biomass$est_se, high = est + 1.96 * biomass$est_se)
  ggplot(d, aes(value)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = part_point_cols[[2]], alpha = 0.2) +
    geom_line(aes(y = est), colour = part_point_cols[[2]], linewidth = 0.8) +
    scale_y_continuous(name = NULL, position = "right") + # title as in the panels above
    scale_x_continuous(expand = c(0, 0), breaks = scales::breaks_pretty(n = 4)) +
    labs(x = "Fishing (swept-area ratio)") +
    theme_light(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title.y.right = element_text(colour = part_point_cols[[2]]),
      axis.text.y.right = element_text(colour = part_point_cols[[2]])
    )
}

# sediment type enters occurrence only: predicted occurrence probability of
# each class with 95% CI, for one grab sample in 2021 with the other
# predictors at their mean (the marginal means emmeans would give at these
# fixed values; emmeans' default would average over survey and year instead)
pe_sediment <- function() {
  sed_labs <- c(mud = "Mud", muddy_sand = "Muddy sand", sand = "Sand", mixed_sediment = "Mixed")
  nd <- tibble(sediment = factor(levels(dat_occ$sediment), levels = levels(dat_occ$sediment)))
  for (col in std_cols) nd[[col]] <- 0
  nd <- mutate(
    nd,
    survey = factor("Stock", levels = levels(dat_occ$survey)),
    year = factor(pe_year, levels = levels(dat_occ$year)),
    x_utm = mean(dat_occ$x_utm), y_utm = mean(dat_occ$y_utm)
  )
  presence <- predict(fits$occurrence, newdata = nd, re_form = NA, se_fit = TRUE, offset = rep(0, nrow(nd)))
  d <- tibble(
    sediment = factor(sed_labs[as.character(nd$sediment)], levels = sed_labs),
    est = inv_cloglog(presence$est),
    low = inv_cloglog(presence$est - 1.96 * presence$est_se),
    high = inv_cloglog(presence$est + 1.96 * presence$est_se)
  )
  ggplot(d, aes(sediment, est)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, linewidth = 0.8, colour = part_point_cols[[1]]) +
    geom_point(size = 2, colour = part_point_cols[[1]]) +
    # class names inside the panel, written vertically down from the top above
    # each interval, so the x axis lines up with the other panels of the row
    geom_text(aes(y = max(high) * 1.7, label = sediment), angle = 90, hjust = 1, size = 2.5, colour = "grey25") +
    scale_y_continuous(name = NULL, limits = c(0, max(d$high) * 1.72)) +
    labs(x = "Sediment") +
    theme_light(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y.left = element_text(colour = part_point_cols[[1]]),
      axis.text.y.left = element_text(colour = part_point_cols[[1]])
    )
}

pe_plots <- c(map(names(pe_vars), pe_plot), list(pe_sediment(), pe_fishing()))
# the two y-axis titles, drawn once at the left and right of the whole grid,
# so they sit level with its middle row
axis_label <- function(label, colour, rot) {
  wrap_elements(full = grid::textGrob(label, rot = rot, gp = grid::gpar(col = colour, fontsize = 10)))
}
p_partial <- axis_label("Occurrence probability", part_point_cols[[1]], 90) +
  wrap_plots(pe_plots, ncol = 3) +
  axis_label("Log biomass where present", part_point_cols[[2]], -90) +
  plot_layout(widths = unit(c(0.5, 1, 0.5), c("cm", "null", "cm")))

# 03 Assemble and save ----
# both panels are wrapped as single elements and tagged at the figure level,
# so the tags sit at the same left edge and panel b's axis titles sit next to
# its plots rather than being pushed out by panel a's wide labels
fig3 <- wrap_elements(full = p_coeff) / wrap_elements(full = p_partial) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = list(c("a)", "b)"))) &
  theme(plot.tag = element_text(size = 12, face = "bold", colour = "grey20"))

dir.create(here("output", "figs", "main"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "main", "fig3_coeff.png"), fig3, width = 7.5, height = 9.5, dpi = 600, bg = "white")
