# Model-selection tables (LaTeX) for the paper.
#
# Table 1: the four main models (space, + environment, + in-strength, +
# both), all delta-gamma with survey, year and a spatial field, ordered by
# AIC, with the out-of-sample ELPD difference to the best model from the
# spatially blocked cross-validation.
# Table S1: the full model with each of the five connectivity metrics.
# Written to output/tables for \input into Overleaf.
#
# Reads data/sdm/main/model_comparison.rds and
# data/sdm/sensitivity/model_comparison.rds (07_fit_sdm.R) and
# data/sdm/cv/elpd_compare.rds (08_cross_validation.R).

library(tidyverse)
library(here)

structure_labels <- c(
  space = "Space",
  space_env = "Space + Environment",
  space_conn = "Space + Connectivity",
  space_env_conn = "Space + Environment + Connectivity"
)
metric_labels <- c(
  space_env_conn = "In-strength",
  space_env_conn_in_degree = "In-degree",
  space_env_conn_in_closeness = "In-closeness",
  space_env_conn_eigen = "Eigenvector centrality",
  space_env_conn_transitivity = "Transitivity"
)

# AIC to one decimal place, the best model (delta = 0) in bold; a model that
# did not converge is shown as a dash
fmt <- function(x, best = FALSE) {
  s <- if_else(is.na(x), "---", formatC(x, format = "f", digits = 1))
  if_else(rep_len(best, length(s)), paste0("\\textbf{", s, "}"), s)
}

dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)

# 01 Table 1: main models ----
comp <- readRDS(here("data", "sdm", "main", "model_comparison.rds"))
elpd <- readRDS(here("data", "sdm", "cv", "elpd_compare.rds"))

tab1 <- comp |>
  left_join(select(elpd, model, elpd_diff, se_diff), by = "model") |>
  mutate(
    label = structure_labels[model],
    delta_aic = aic - min(aic, na.rm = TRUE),
    best = !is.na(delta_aic) & delta_aic == 0
  ) |>
  arrange(delta_aic)

rows1 <- tab1 |>
  transmute(line = paste0(
    label, " & ", fmt(aic, best), " & ", fmt(delta_aic, best), " & ",
    fmt(elpd_diff, elpd_diff == 0), " (", fmt(se_diff), ") \\\\"
  )) |>
  pull(line)

writeLines(c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Model selection for cockle biomass (delta-gamma hurdle models), ordered by $\\Delta$AIC (lower is better; the best model is shown in bold). $\\Delta$ELPD is the difference in expected log predictive density to the best model in spatially blocked 10-fold cross-validation (higher is better), with its standard error in parentheses. All models include survey (gear) and year as factors and a spatial random field. Connectivity is presence-weighted in-strength; environment comprises depth (linear and quadratic), temperature, oxygen, salinity, and maximum shear stress. A dash marks a model that did not converge.}",
  "\\label{tab:model-selection}",
  "\\begin{tabular}{lrrr}",
  "\\toprule",
  "Model & AIC & $\\Delta$AIC & $\\Delta$ELPD (SE) \\\\",
  "\\midrule",
  rows1,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
), here("output", "tables", "tab1_model_selection.tex"))

# 02 Table S1: full model with each connectivity metric ----
tab_s1 <- bind_rows(
  filter(comp, model == "space_env_conn"),
  readRDS(here("data", "sdm", "sensitivity", "model_comparison.rds"))
) |>
  mutate(
    label = metric_labels[model],
    delta_aic = aic - min(aic, na.rm = TRUE),
    best = !is.na(delta_aic) & delta_aic == 0
  ) |>
  arrange(delta_aic)

rows_s1 <- tab_s1 |>
  transmute(line = paste0(label, " & ", fmt(aic, best), " & ", fmt(delta_aic, best), " \\\\")) |>
  pull(line)

writeLines(c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{AIC of the full Space + Environment + Connectivity model with each of the five connectivity metrics, ordered by $\\Delta$AIC (the best model is shown in bold). All metrics are presence-weighted; all models include survey, year, the environmental predictors and a spatial random field.}",
  "\\label{tab:connectivity-metrics}",
  "\\begin{tabular}{lrr}",
  "\\toprule",
  "Connectivity metric & AIC & $\\Delta$AIC \\\\",
  "\\midrule",
  rows_s1,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
), here("output", "tables", "tabS1_connectivity_metrics.tex"))
