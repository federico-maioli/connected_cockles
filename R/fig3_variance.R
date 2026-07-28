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

part_cols <- c(
  "Spatial field" = "#4C72B0", "Environment" = "#2E8B57",
  "Survey (gear)" = "#7F7F7F", "Connectivity" = "#E58606", "Unexplained" = "#ECECEC"
)

partition_plot <- function(d) {
  # drop components that are structurally absent (no field when spatial = off)
  d <- d |> filter(share > 0)
  ggplot(d, aes(x = share, y = metric_lab, fill = component)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = if_else(share >= 0.05, scales::percent(share, accuracy = 0.1), "")),
      position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.9, colour = "grey10"
    ) +
    scale_fill_manual(values = part_cols, name = NULL, drop = TRUE) +
    scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
    labs(x = expression("Share of total variance (Nakagawa " * R^2 * ")"), y = NULL) +
    theme_light(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 1))
}

# 01 Load shares ----
metric_labs <- c(
  log_biomass_in_strength = "Biomass in-strength (log)",
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
  width = 9.5, height = 4.2, dpi = 600, bg = "white"
)
ggsave(here("output", "figs", "figS6_variance_nospace.png"), p_nospace,
  width = 9.5, height = 4.2, dpi = 600, bg = "white"
)
