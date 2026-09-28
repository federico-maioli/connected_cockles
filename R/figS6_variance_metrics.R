# Supplementary figure: variance partitioning across all 5 connectivity
# metrics (spatial field on, both responses) - shows the small connectivity
# share seen for in-strength (fig2_coeff.R panel b) is not particular to that
# one metric.
#
# Reads the per-model variance tables written by 02_fit_sdm.R:
#   data/sdm/main/<response>_space_env_conn_variance.rds             (in-strength)
#   data/sdm/sensitivity/<response>_space_env_conn_<metric>_variance.rds  (the other 4)

library(tidyverse)
library(here)

# Okabe-Ito, colourblind-safe; unexplained variance stays neutral grey
part_cols <- c(
  "Spatial field" = "#0072B2", "Environment" = "#009E73",
  "Connectivity" = "#E69F00", "Unexplained" = "#E6E6E6"
)
label_cols <- c(
  "Spatial field" = "white", "Environment" = "white",
  "Connectivity" = "grey15", "Unexplained" = "grey15"
)

# one stacked bar per level of `group`, faceted on `facet_var`; d must have
# share, component, group. Segments >= 5% are labelled in place; segments too
# thin to label are bundled into one line of text to the right of the bar
# (avoids the label crowding/collision a repel-based approach runs into with
# many thin segments across many bars - see fig2_coeff.R for the simpler,
# two-bar version of this same idea)
partition_plot <- function(d, facet_var) {
  d <- d |> filter(share > 0)
  by_vars <- c("group", facet_var)

  labs_d <- d |>
    arrange(pick(all_of(by_vars)), component) |>
    mutate(x_mid = cumsum(share) - share / 2, .by = all_of(by_vars))

  wide_labs <- labs_d |> filter(share >= 0.05)
  narrow_labs <- labs_d |>
    filter(share < 0.05) |>
    summarise(
      label = paste(component, scales::percent(share, accuracy = 0.1), collapse = "\n"),
      .by = all_of(by_vars)
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
      aes(x = 1.02, y = group, label = label),
      inherit.aes = FALSE, size = 2.6, colour = "grey25", hjust = 0, vjust = 0.5, lineheight = 0.85
    ) +
    scale_fill_manual(values = part_cols, name = NULL, drop = TRUE) +
    scale_colour_manual(values = label_cols, guide = "none") +
    scale_x_continuous(
      breaks = seq(0, 1, 0.25), labels = scales::percent,
      expand = expansion(mult = c(0, 0.4))
    ) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
    facet_wrap(vars(.data[[facet_var]]), scales = "free_x") +
    labs(x = expression("Share of total variance (Nakagawa " * R^2 * ")"), y = NULL) +
    theme_light(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(colour = "black")
    ) +
    guides(fill = guide_legend(nrow = 1))
}

# 01 Load the space_env_conn variance table for every metric, both responses ----
metric_labs <- c(
  in_degree = "In-degree", in_strength = "In-strength", in_closeness = "In-closeness",
  eigen = "Eigenvector centrality", transitivity = "Transitivity"
)

variance_file <- function(resp, metric) {
  if (metric == "in_strength") {
    here("data", "sdm", "main", paste0(resp, "_space_env_conn_variance.rds"))
  } else {
    here("data", "sdm", "sensitivity", paste0(resp, "_space_env_conn_", metric, "_variance.rds"))
  }
}

sens <- expand_grid(response = c("biomass", "present"), metric = names(metric_labs)) |>
  mutate(variance = map2(response, metric, \(resp, m) readRDS(variance_file(resp, m)))) |>
  unnest(variance) |>
  mutate(
    component = recode(component, spatial = "Spatial field", distribution = "Unexplained"),
    component = factor(component, levels = names(part_cols)),
    group = factor(metric_labs[metric], levels = rev(unname(metric_labs))),
    response_lab = recode(response, biomass = "Biomass (log link)", present = "Presence (logit link)")
  )

# 02 Plot ----
p_sens <- partition_plot(sens, facet_var = "response_lab")

# 03 Save ----
dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "figS6_variance_metrics.png"), p_sens,
  width = 9.5, height = 4.4, dpi = 600, bg = "white"
)
