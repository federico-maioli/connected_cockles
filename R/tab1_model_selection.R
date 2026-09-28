# Model-selection table (LaTeX) for the paper.
#
# The six main structures (space, environment, and connectivity in-strength,
# crossed) for the biomass (Tweedie) response only, ordered by delta AIC.
# Writes a booktabs table to output/tables for \input into Overleaf.
#
# Reads data/sdm/main/model_comparison.rds, written by 02_fit_sdm.R.

library(tidyverse)
library(here)

# 01 Load comparison ----
comp <- readRDS(here("data", "sdm", "main", "model_comparison.rds")) |>
  filter(response == "biomass")

# 02 Delta AIC, ordered best first ----
structure_labels <- c(
  space = "Space",
  env = "Environment",
  space_env = "Space + Environment",
  space_conn = "Space + Connectivity",
  env_conn = "Environment + Connectivity",
  space_env_conn = "Space + Environment + Connectivity"
)

tab <- comp |>
  mutate(
    label = structure_labels[model],
    delta_aic = aic - min(aic, na.rm = TRUE)
  ) |>
  arrange(delta_aic) |>
  select(label, aic, delta_aic)

# 03 Build LaTeX ----
# one decimal place; the best model (delta = 0) is bold; a model that did not
# converge has no AIC and is shown as a dash
fmt <- function(x) {
  best <- !is.na(x) & x == min(x, na.rm = TRUE)
  s <- formatC(x, format = "f", digits = 1)
  s <- if_else(is.na(x), "---", s)
  if_else(best, paste0("\\textbf{", s, "}"), s)
}
tab <- tab |>
  mutate(aic_f = fmt(aic), delta_aic_f = fmt(delta_aic))

rows <- tab |>
  transmute(line = paste0(label, " & ", aic_f, " & ", delta_aic_f, " \\\\")) |>
  pull(line)

latex <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Model selection for cockle biomass (Tweedie), ordered by $\\Delta$AIC (lower is better; the best model is shown in bold). Connectivity is in-strength throughout; environment comprises depth, temperature, oxygen, salinity, and maximum shear stress; space is a spatial random field. A dash marks a model that did not converge.}",
  "\\label{tab:model-selection}",
  "\\begin{tabular}{lrr}",
  "\\toprule",
  "Model & AIC & $\\Delta$AIC \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

# 04 Write ----
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
writeLines(latex, here("output", "tables", "tab1_model_selection.tex"))
