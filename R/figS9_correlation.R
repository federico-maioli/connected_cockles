# Supplementary figure: correlation among all predictors used in 02_fit_sdm.R,
# environment and connectivity together (not split into separate blocks as in
# the VIF check there), so any cross-block collinearity is visible too.

library(tidyverse)
library(here)

# 01 Load data ----
dat <- readRDS(here("data", "derived", "cockles_clean.rds"))

pred_vars <- c(
  depth = "depth", temp = "temp", sal = "sal", oxy = "oxy", shear_max = "shear_max",
  in_degree = "conn_in_degree", in_strength = "conn_in_strength",
  in_closeness = "conn_in_closeness", eigen = "conn_eigen", transitivity = "conn_transitivity"
)

# 02 Correlation matrix, long format for plotting ----
cor_mat <- dat |>
  select(all_of(pred_vars)) |>
  cor(use = "complete.obs")

cor_long <- cor_mat |>
  as_tibble(rownames = "var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "r") |>
  mutate(
    var1 = factor(var1, levels = names(pred_vars)),
    var2 = factor(var2, levels = rev(names(pred_vars)))
  )

# 03 Plot ----
p <- ggplot(cor_long, aes(var1, var2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.8) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), name = "Pearson r") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_light(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# 04 Save ----
dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "figS9_correlation.png"), p,
  width = 7, height = 6, dpi = 600, bg = "white"
)
