# Supplementary figure: does the biomass-connectivity relationship depend on
# which dispersal year the connectivity was computed from?
#
# Each connectivity metric is extracted at the survey points once per flow matrix
# (2010-2016) and once for the pooled run, then correlated with observed cockle
# biomass. If the correlations sit in a tight band across years, the pooled run
# used in the models is representative and the choice of dispersal year is not
# driving the result.
#
# Spearman correlation, because biomass is strongly right-skewed (Tweedie) and
# the metrics are on very different scales.

library(tidyverse)
library(here)

# 01 Load the per-run extraction ----
sens <- readRDS(here("data", "derived", "cockles_connectivity_sensitivity.rds")) |>
  filter(!is.na(biomass))

metrics <- c(
  "log_biomass_in_strength", "presence_in_strength", "deg_in", "in_strength",
  "eigen_centrality", "closeness_centrality"
)
metric_labs <- c(
  log_biomass_in_strength = "Biomass (log) in-strength",
  presence_in_strength = "Presence in-strength",
  deg_in = "In-degree",
  in_strength = "In-strength",
  eigen_centrality = "Eigenvector centrality",
  closeness_centrality = "Closeness centrality"
)

# 02 Correlation with biomass, per run ----
cors <- sens |>
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") |>
  summarise(
    rho = cor(biomass, value, method = "spearman"),
    .by = c(run, metric)
  ) |>
  mutate(metric = factor(metric_labs[metric], levels = unname(metric_labs)))

# pooled run drawn as a reference line, the years as points
pooled <- cors |> filter(run == "all")
yearly <- cors |> filter(run != "all") |> mutate(year = as.integer(run))

# 03 Report the spread ----
yearly |>
  summarise(
    min_rho = round(min(rho), 3),
    max_rho = round(max(rho), 3),
    range = round(max(rho) - min(rho), 3),
    .by = metric
  ) |>
  left_join(pooled |> transmute(metric, pooled_rho = round(rho, 3)), by = "metric") |>
  arrange(range) |>
  print()

# 04 Plot ----
p <- ggplot(yearly, aes(year, rho)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_hline(
    data = pooled, aes(yintercept = rho),
    colour = "#D55E00", linewidth = 0.5
  ) +
  geom_line(colour = "grey40", linewidth = 0.4) +
  geom_point(size = 2, colour = "grey20") +
  facet_wrap(~metric, nrow = 2) +
  scale_x_continuous(breaks = c(2011, 2013, 2015), expand = expansion(mult = 0.12)) +
  labs(
    x = "Dispersal year of the flow matrix",
    y = "Spearman correlation with biomass"
  ) +
  theme_light(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "grey20", face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# 05 Save ----
dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "figS4_year_sensitivity.png"), p,
  width = 8.5, height = 5.6, dpi = 600, bg = "white"
)
