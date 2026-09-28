# Figure 2: the full model (space + environment + in-strength connectivity),
# with and without the spatial field, two panels stacked with patchwork.
#
# Panel A: standardized coefficient of every predictor, with vs without the
# field - shows how much each one moves when the field is dropped. In-strength
# is the only connectivity metric whose coefficient survives this (the other
# four are checked the same way in figS5_space_confounding.R, not repeated
# here).
# Panel B: variance partitioning of the same two models (Nakagawa R2, split
# HMSC-style among Environment/Connectivity - see 02_fit_sdm.R and
# nakagawa_sdmtmb() in R/helpers.R) - shows how much of panel A's stability
# comes from the field soaking up variance that would otherwise land on the
# fixed effects.
#
# Both panels read data/sdm/main/<response>_space_env_conn(.rds/_variance.rds)
# and .../<response>_env_conn(.rds/_variance.rds) - the full in-strength model,
# fitted in 02_fit_sdm.R.
#
# Writes two files from the same builder:
#   fig2_coeff.png             biomass (Tweedie)   - main text
#   figS3_coeff_presence.png   presence (binomial) - supplementary

library(tidyverse)
library(here)
library(sdmTMB)
library(ggstats)
library(patchwork)

# the two models behind both panels: full model, spatial field on vs off
structure_labs <- c(space_env_conn = "Spatial field: on", env_conn = "Spatial field: off")
field_cols <- c("Spatial field: on" = "#0072B2", "Spatial field: off" = "#D55E00")

term_labs <- c(
  depth_std = "Depth", temp_std = "Temperature", oxy_std = "Oxygen",
  sal_std = "Salinity", shear_max_std = "Shear stress",
  conn_in_strength_std = "In-strength"
)
term_levels <- c("In-strength", "Shear stress", "Salinity", "Oxygen", "Temperature", "Depth")

# Okabe-Ito, colourblind-safe; unexplained variance stays neutral grey
part_cols <- c(
  "Spatial field" = "#0072B2", "Environment" = "#009E73",
  "In-strength" = "#E69F00", "Unexplained" = "#E6E6E6"
)
label_cols <- c(
  "Spatial field" = "white", "Environment" = "white",
  "In-strength" = "grey15", "Unexplained" = "grey15"
)

# 01 Panel A: coefficients, with vs without the spatial field ----
coeff_panel <- function(resp, xlab) {
  d <- map(names(structure_labs), function(structure) {
    fit <- readRDS(here("data", "sdm", "main", paste0(resp, "_", structure, ".rds")))
    tidy(fit, effects = "fixed", conf.int = TRUE) |>
      mutate(field = structure_labs[[structure]])
  }) |>
    list_rbind() |>
    filter(term %in% names(term_labs)) |>
    mutate(
      term = factor(term_labs[term], levels = term_levels),
      field = factor(field, levels = unname(structure_labs))
    )

  dodge <- position_dodge(width = 0.55)
  ggplot(d, aes(estimate, term, colour = field)) +
    geom_stripped_rows(colour = NA) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0, linewidth = 0.9, position = dodge) +
    geom_point(size = 2.4, position = dodge) +
    scale_colour_manual(values = field_cols, name = NULL) +
    labs(x = xlab, y = NULL) +
    theme_light(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

# 02 Panel B: variance partitioning, with vs without the spatial field ----
# one stacked bar per level of `group`; segments >= 5% are labelled in place,
# thinner ones (e.g. Connectivity) get their value placed just above the slice
partition_plot <- function(d) {
  d <- d |> filter(share > 0)

  labs_d <- d |>
    arrange(group, component) |>
    mutate(x_mid = cumsum(share) - share / 2, .by = group)

  wide_labs <- labs_d |> filter(share >= 0.05)
  # just the value - the fill colour already says which component this is
  narrow_labs <- labs_d |>
    filter(share < 0.05) |>
    summarise(
      label = paste(scales::percent(share, accuracy = 0.1), collapse = "\n"),
      x_mid = mean(x_mid),
      .by = group
    )

  ggplot(d, aes(x = share, y = group, fill = component)) +
    geom_col(width = 0.55, colour = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
    geom_text(
      data = wide_labs,
      aes(x = x_mid, y = group, label = scales::percent(share, accuracy = 0.1), colour = component),
      inherit.aes = FALSE, size = 3, show.legend = FALSE
    ) +
    geom_text(
      data = narrow_labs,
      aes(x = x_mid, y = group, label = label),
      inherit.aes = FALSE, size = 2.6, colour = "grey25", lineheight = 0.85,
      position = position_nudge(y = 0.35)
    ) +
    scale_fill_manual(values = part_cols, name = NULL, drop = TRUE) +
    scale_colour_manual(values = label_cols, guide = "none") +
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = scales::percent, expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
    labs(x = expression("Share of total variance (Nakagawa " * R^2 * ")"), y = NULL) +
    theme_light(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 1))
}

variance_panel <- function(resp) {
  d <- map(names(structure_labs), function(structure) {
    readRDS(here("data", "sdm", "main", paste0(resp, "_", structure, "_variance.rds"))) |>
      mutate(
        # block name from 02_fit_sdm.R's blocks_env_conn_in_strength is
        # "Connectivity" (shared across all 5 metrics there); this figure is
        # in-strength only, so it gets the more specific label here
        component = recode(component,
          spatial = "Spatial field", distribution = "Unexplained", Connectivity = "In-strength"
        ),
        group = structure_labs[[structure]]
      )
  }) |>
    list_rbind() |>
    filter(share > 0) |>
    mutate(
      component = factor(component, levels = names(part_cols)),
      # "on" first, at the top of the bar chart (ggplot draws the first
      # discrete-axis level at the bottom, so the level order is reversed)
      group = factor(group, levels = rev(unname(structure_labs)))
    )
  partition_plot(d)
}

# 03 Assemble both panels ----
build_fig <- function(resp, xlab) {
  coeff_panel(resp, xlab) / variance_panel(resp) +
    plot_annotation(tag_levels = "a", tag_suffix = ")")
}

# 04 Biomass (main text) ----
p_biomass <- build_fig("biomass", "Standardized coefficient (log link)")

# 05 Presence (supplementary) ----
p_present <- build_fig("present", "Standardized coefficient (logit link)")

# 06 Save ----
dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "fig2_coeff.png"), p_biomass,
  width = 7, height = 8, dpi = 600, bg = "white"
)
ggsave(here("output", "figs", "figS3_coeff_presence.png"), p_present,
  width = 7, height = 8, dpi = 600, bg = "white"
)
