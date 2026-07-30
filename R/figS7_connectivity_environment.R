# Supplementary figure: what does connectivity covary with?
#
# Two questions in one panel:
#   - are the connectivity metrics collinear with the environmental predictors?
#   - are they simply standing in for residual spatial structure?
# The last column is the spatial random field estimated from the Space +
# Environment model, i.e. the structure connectivity would have to explain if it
# were acting as a spatial smoother. It is not an environmental covariate, hence
# the separator.
#
# Everything is computed at the survey stations, the level at which the models
# are fitted and the VIFs in 06 are computed. Survey effort is uneven across grid
# cells (7468 stations in 324 cells) and connectivity is constant within a cell,
# so these correlations are weighted towards heavily sampled cells.

library(tidyverse)
library(here)
library(sdmTMB)

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
  shear_max = "Shear stress",
  spatial_field = "Spatial field"
)

# 01 Load ----
# the fitted data carries both the covariates and the model rows, so the field
# and the covariates are guaranteed to align
fit <- readRDS(here("data", "derived", "sdm_fits.rds"))[["biomass_space_env"]]
dat <- fit$data
dat$spatial_field <- predict(fit)$est_rf

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

cat("correlation with the residual spatial field:\n")
cors |>
  filter(covariate == "Spatial field") |>
  mutate(rho = round(rho, 2)) |>
  arrange(rho) |>
  print()

# 03 Plot ----
p <- ggplot(cors, aes(covariate, metric, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.1, colour = "black") +
  # the spatial field is not an environmental covariate
  geom_vline(xintercept = 5.5, colour = "grey30", linewidth = 0.5) +
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
  width = 6.4, height = 5.2, dpi = 600, bg = "white"
)
