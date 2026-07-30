# Supplementary figure: how strongly does connectivity covary with the
# environmental predictors?
#
# This is the evidence behind the interpretation of the negative connectivity
# coefficient - that well-connected cells are also the deep, brackish, low-energy
# parts of the fjord, so connectivity and habitat are not independent predictors.
#
# Computed at the survey stations, which is what the fitted models see and the
# level at which the VIFs in 06 are computed. Note that survey effort is uneven
# across grid cells (7468 stations in 324 cells) and connectivity is constant
# within a cell, so these correlations are weighted towards heavily sampled cells.

library(tidyverse)
library(here)

conn_labs <- c(
  log_biomass_in_strength = "Biomass (log) in-strength",
  presence_in_strength = "Presence in-strength",
  deg_in = "In-degree",
  in_strength = "In-strength",
  eigen_centrality = "Eigenvector centrality",
  closeness_centrality = "Closeness centrality"
)
env_labs <- c(
  depth = "Depth",
  temp = "Temperature",
  oxy = "Oxygen",
  sal = "Salinity",
  shear_max = "Shear stress"
)

# 01 Load ----
dat <- readRDS(here("data", "derived", "cockles_connectivity.rds"))

# 02 Correlations ----
# Spearman throughout: the metrics are strongly right-skewed
cors <- cor(
  dat[, names(conn_labs)], dat[, names(env_labs)],
  method = "spearman", use = "pairwise.complete.obs"
) |>
  as.data.frame() |>
  rownames_to_column("metric") |>
  pivot_longer(-metric, names_to = "covariate", values_to = "rho") |>
  mutate(
    metric = factor(conn_labs[metric], levels = rev(unname(conn_labs))),
    covariate = factor(env_labs[covariate], levels = unname(env_labs))
  )

print(cors |> arrange(desc(abs(rho))) |> mutate(rho = round(rho, 2)), n = 8)

# 03 Plot ----
p <- ggplot(cors, aes(covariate, metric, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.1, colour = "black") +
  scale_fill_distiller(
    palette = "RdBu", limits = c(-1, 1), direction = 1,
    name = "Spearman correlation", guide = guide_colourbar(title.position = "top")
  ) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_light(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    legend.key.height = unit(3, "mm"), legend.key.width = unit(14, "mm")
  )

# 04 Save ----
dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "figS7_connectivity_environment.png"), p,
  width = 5.6, height = 5.2, dpi = 600, bg = "white"
)
