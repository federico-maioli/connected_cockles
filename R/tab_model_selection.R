# Model-selection table (LaTeX) for the paper.
#
# All twelve model structures (the connectivity slot rotated over three
# predictors), ordered by biomass Delta AIC. Writes a booktabs table to
# output/tables for \input into Overleaf.

library(tidyverse)
library(here)

# 01 Load comparison ----
comp <- readRDS(here("data", "final", "sdm_model_comparison.rds"))

# 02 Build meaningful labels for every model ----
# base structure name + connectivity type (in parentheses) where present
structure_labels <- c(
  space = "Space",
  env = "Environment",
  space_env = "Space + Environment",
  space_conn = "Space + Connectivity",
  env_conn = "Environment + Connectivity",
  space_env_conn = "Space + Environment + Connectivity"
)
conn_labels <- c(
  in_strength = "in-strength",
  biomass_supply_log = "biomass supply, log",
  eigen = "eigenvector centrality"
)

wide <- comp |>
  select(response, model, structure, delta_aic) |>
  pivot_wider(names_from = response, values_from = delta_aic) |>
  mutate(
    conn = str_remove(model, paste0("^", structure, "_?")),
    label = if_else(
      conn == "",
      structure_labels[structure],
      paste0(structure_labels[structure], " (", conn_labels[conn], ")")
    )
  )

# order by biomass Delta AIC (best first)
tab <- wide |>
  arrange(biomass) |>
  select(label, present, biomass)

# 03 Build LaTeX ----
# bold the best model in each response (delta = 0); one decimal place
fmt <- function(x, best) {
  s <- formatC(x, format = "f", digits = 1)
  if_else(best, paste0("\\textbf{", s, "}"), s)
}
tab <- tab |>
  mutate(
    present_f = fmt(present, present == min(present)),
    biomass_f = fmt(biomass, biomass == min(biomass))
  )

rows <- tab |>
  transmute(line = paste0(label, " & ", present_f, " & ", biomass_f, " \\\\")) |>
  pull(line)

latex <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Model selection for cockle presence and biomass, ordered by biomass $\\Delta$AIC. Values are $\\Delta$AIC relative to the best model within each response (lower is better; the best model in each column is shown in bold). Every model includes a survey (gear) intercept; environment comprises depth, temperature, oxygen, salinity, and maximum shear stress; space is a spatial random field; the connectivity metric is given in parentheses.}",
  "\\label{tab:model-selection}",
  "\\begin{tabular}{lrr}",
  "\\toprule",
  "Model & Presence $\\Delta$AIC & Biomass $\\Delta$AIC \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

# 04 Write ----
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
writeLines(latex, here("output", "tables", "model_selection.tex"))
