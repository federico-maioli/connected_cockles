# Table 2: variance partitioning (LaTeX), components as rows and connectivity
# metrics as columns, for the full spatial model.
#
# Shares are computed in 07_variance_partitioning.R. Writes a booktabs table to
# output/tables for \input into the manuscript.

library(tidyverse)
library(here)

# 01 Load shares ----
metric_labs <- c(
  log_biomass_in_strength = "Biomass in-strength (log)",
  deg_in = "In-degree",
  in_strength = "In-strength",
  eigen_centrality = "Eigenvector centrality",
  closeness_centrality = "Closeness centrality"
)
component_order <- c(
  "Spatial field", "Environment", "Survey (gear)", "Connectivity", "Unexplained"
)

parts <- readRDS(here("data", "derived", "variance_partition.rds")) |>
  filter(structure == "space_env_conn")

# 02 Wide layout ----
metric_order <- names(metric_labs)
wide <- parts |>
  mutate(
    pct = sprintf("%.1f", 100 * share),
    component = factor(component, levels = component_order)
  ) |>
  select(component, metric, pct) |>
  pivot_wider(names_from = metric, values_from = pct) |>
  arrange(component)

header <- paste0("Component & ", paste(metric_labs[metric_order], collapse = " & "), " \\\\")
rows <- wide |>
  select(component, all_of(metric_order)) |>
  as.matrix() |>
  apply(1, \(r) paste0(paste(r, collapse = " & "), " \\\\"))

# mean R2 across metrics, quoted in the caption
r2 <- parts |>
  filter(component == "Connectivity") |>
  summarise(cond = mean(conditional_r2), marg = mean(marginal_r2))

# 03 Build LaTeX ----
latex <- c(
  "\\begin{table}[ht]",
  "\\centering",
  sprintf("\\caption{Variance partitioning of cockle biomass from the spatial Tweedie SDM (Nakagawa \\& Schielzeth marginal / conditional $R^2$), with the connectivity slot rotated over five metrics. The spatial field and the three fixed-effect blocks together give the conditional $R^2$ (mean %.1f\\%% across metrics); the fixed blocks alone give the marginal $R^2$ (mean %.1f\\%%). Shares are percentages of total variance and sum to 100\\%% within each column.}", 100 * r2$cond, 100 * r2$marg),
  "\\label{tab:variance-partitioning}",
  paste0("\\begin{tabular}{l", strrep("r", length(metric_order)), "}"),
  "\\toprule",
  header,
  "\\midrule",
  rows[1:4],
  "\\midrule",
  rows[5],
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

# 04 Write ----
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
writeLines(latex, here("output", "tables", "tab2_variance_partitioning.tex"))
