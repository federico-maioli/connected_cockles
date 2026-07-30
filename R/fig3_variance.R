# Figure 3: variance partitioning of cockle biomass, one bar per connectivity
# metric, from the full spatial model.
#
# Also writes the supplementary counterpart without the spatial random field
# (figS6), which shares the same builder: with the field removed, whatever it was
# absorbing is redistributed onto the fixed blocks.
#
# Shares are computed in 07_variance_partitioning.R.

library(tidyverse)
library(here)
library(ggrepel)

# Okabe-Ito, colourblind-safe; unexplained variance stays neutral grey
part_cols <- c(
  "Spatial field" = "#0072B2", "Environment" = "#009E73",
  "Survey (gear)" = "#CC79A7", "Connectivity" = "#E69F00", "Unexplained" = "#E6E6E6"
)

# text colour for labels that sit inside their segment
label_cols <- c(
  "Spatial field" = "white", "Environment" = "white",
  "Survey (gear)" = "grey15", "Connectivity" = "grey15", "Unexplained" = "grey15"
)

partition_plot <- function(d) {
  # drop components that are structurally absent (no field when spatial = off)
  d <- d |> filter(share > 0)

  # segment midpoints; wide segments are labelled in place, narrow ones (< 5%)
  # are lifted above the bar and repelled sideways so they stay readable
  labs_d <- d |>
    arrange(metric_lab, component) |>
    mutate(x_mid = cumsum(share) - share / 2, .by = metric_lab) |>
    mutate(narrow = share < 0.05)

  ggplot(d, aes(x = share, y = metric_lab, fill = component)) +
    geom_col(width = 0.62, colour = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
    geom_text(
      data = filter(labs_d, !narrow),
      aes(x = x_mid, y = metric_lab, label = scales::percent(share, accuracy = 0.1), colour = component),
      inherit.aes = FALSE, size = 2.9, show.legend = FALSE
    ) +
    geom_text_repel(
      data = filter(labs_d, narrow),
      aes(x = x_mid, y = metric_lab, label = scales::percent(share, accuracy = 0.1)),
      inherit.aes = FALSE, size = 2.9, colour = "grey25",
      nudge_y = 0.55, direction = "x", min.segment.length = 0,
      segment.size = 0.3, segment.colour = "black", box.padding = 0.12, seed = 1
    ) +
    scale_fill_manual(values = part_cols, name = NULL, drop = TRUE) +
    scale_colour_manual(values = label_cols, guide = "none") +
    scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
    scale_y_discrete(expand = expansion(add = c(0.45, 0.9))) +
    labs(x = expression("Share of total variance (Nakagawa " * R^2 * ")"), y = NULL) +
    theme_light(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 1))
}

# 01 Load shares ----
metric_labs <- c(
  log_biomass_in_strength = "Biomass (log) in-strength",
  presence_in_strength = "Presence in-strength",
  deg_in = "In-degree",
  in_strength = "In-strength",
  eigen_centrality = "Eigenvector centrality",
  closeness_centrality = "Closeness centrality"
)

parts <- readRDS(here("data", "derived", "variance_partition.rds")) |>
  mutate(
    metric_lab = factor(metric_labs[metric], levels = rev(unname(metric_labs))),
    component = factor(component, levels = c(
      "Spatial field", "Environment", "Survey (gear)", "Connectivity", "Unexplained"
    ))
  )

# 02 Build and save ----
p_space <- partition_plot(filter(parts, structure == "space_env_conn"))
p_nospace <- partition_plot(filter(parts, structure == "env_conn"))

dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "fig3_variance.png"), p_space,
  width = 9.5, height = 4.8, dpi = 600, bg = "white"
)
ggsave(here("output", "figs", "figS6_variance_nospace.png"), p_nospace,
  width = 9.5, height = 4.8, dpi = 600, bg = "white"
)
