# Model-selection table (LaTeX) for the paper.
#
# The six main structures (space, environment, and connectivity in-strength,
# crossed), all delta-gamma models, ordered by delta AIC, with the marginal
# and conditional R2 of each component (presence; biomass where present).
# Writes a booktabs table to output/tables for \input into Overleaf.
#
# Reads data/sdm/main/model_comparison.rds, written by 07_fit_sdm.R.

library(tidyverse)
library(here)

# 01 Load comparison ----
comp <- readRDS(here("data", "sdm", "main", "model_comparison.rds"))

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
  arrange(delta_aic)

# 03 Build LaTeX ----
# AIC to one decimal place, the best model (delta = 0) in bold; R2 to two
# decimals; a model that did not converge is shown as a dash
fmt_aic <- function(x) {
  best <- !is.na(x) & x == min(x, na.rm = TRUE)
  s <- if_else(is.na(x), "---", formatC(x, format = "f", digits = 1))
  if_else(best, paste0("\\textbf{", s, "}"), s)
}
fmt_r2 <- function(x) if_else(is.na(x), "---", formatC(x, format = "f", digits = 2))

rows <- tab |>
  transmute(line = paste0(
    label, " & ", fmt_aic(aic), " & ", fmt_aic(delta_aic), " & ",
    fmt_r2(r2_marginal_presence), " & ", fmt_r2(r2_conditional_presence), " & ",
    fmt_r2(r2_marginal_biomass), " & ", fmt_r2(r2_conditional_biomass), " \\\\"
  )) |>
  pull(line)

latex <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Model selection for cockle biomass (delta-gamma hurdle models), ordered by $\\Delta$AIC (lower is better; the best model is shown in bold). $R^2_m$ and $R^2_c$ are the marginal (fixed effects) and conditional (fixed effects and spatial field) $R^2$ of each component: presence (binomial) and biomass where present (gamma). All models include survey (gear) and year as factors. Connectivity is presence-weighted in-strength throughout; environment comprises depth (linear and quadratic), temperature, oxygen, salinity, and maximum shear stress; space is a spatial random field. A dash marks a model that did not converge.}",
  "\\label{tab:model-selection}",
  "\\begin{tabular}{lrrrrrr}",
  "\\toprule",
  " & & & \\multicolumn{2}{c}{Presence} & \\multicolumn{2}{c}{Biomass where present} \\\\",
  "\\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
  "Model & AIC & $\\Delta$AIC & $R^2_m$ & $R^2_c$ & $R^2_m$ & $R^2_c$ \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

# 04 Write ----
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
writeLines(latex, here("output", "tables", "tab1_model_selection.tex"))
