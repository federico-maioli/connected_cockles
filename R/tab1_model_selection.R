# Model-selection tables (LaTeX) for the paper.
#
# Table 1: the four main models (space, + environment, + in-strength, +
# both), all with survey, year and a spatial field. For each part of the
# delta-gamma model (presence; biomass where present) and in total: the AIC
# difference to the best model and, from the spatially blocked
# cross-validation, the difference in expected log predictive density (ELPD)
# to the best model with its standard error.
# Table S1: the full model with each of the five connectivity metrics, AIC
# difference per part and in total.
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

# one decimal place, the best value (0) in bold; a model that did not converge
# is shown as a dash
fmt <- function(x) {
  s <- if_else(is.na(x), "---", formatC(x, format = "f", digits = 1))
  if_else(!is.na(x) & x == 0, paste0("\\textbf{", s, "}"), s)
}
fmt_elpd <- function(diff, se) if_else(diff == 0, "\\textbf{0.0}", paste0(fmt(diff), " (", formatC(se, format = "f", digits = 1), ")"))

delta_aic <- function(comp) {
  mutate(
    comp,
    d_presence = aic_presence - min(aic_presence, na.rm = TRUE),
    d_biomass = aic_biomass - min(aic_biomass, na.rm = TRUE),
    d_total = aic - min(aic, na.rm = TRUE)
  )
}

dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)

# 01 Table 1: main models ----
elpd <- readRDS(here("data", "sdm", "cv", "elpd_compare.rds")) |>
  select(model, part, elpd_diff, se_diff) |>
  pivot_wider(names_from = part, values_from = c(elpd_diff, se_diff))

tab1 <- readRDS(here("data", "sdm", "main", "model_comparison.rds")) |>
  delta_aic() |>
  left_join(elpd, by = "model") |>
  mutate(label = structure_labels[model]) |>
  arrange(d_total)

rows1 <- tab1 |>
  transmute(line = paste0(
    label, " & ",
    fmt(d_presence), " & ", fmt(d_biomass), " & ", fmt(d_total), " & ",
    fmt_elpd(elpd_diff_presence, se_diff_presence), " & ",
    fmt_elpd(elpd_diff_biomass, se_diff_biomass), " & ",
    fmt_elpd(elpd_diff_total, se_diff_total), " \\\\"
  )) |>
  pull(line)

writeLines(c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\small",
  "\\caption{Model selection for cockle presence, biomass where present and both together (the delta-gamma model). $\\Delta$AIC is the difference in AIC to the best model (lower AIC is better). $\\Delta$ELPD is the difference in expected log predictive density to the best model in spatially blocked 10-fold cross-validation (higher ELPD is better), with its standard error in parentheses. The best model in each column is shown in bold. The total AIC and ELPD are the sums of the two parts. All models include survey (gear) and year as factors and a spatial random field. Connectivity is presence-weighted in-strength; environment comprises depth (linear and quadratic), temperature, oxygen, salinity, and maximum shear stress. A dash marks a model that did not converge.}",
  "\\label{tab:model-selection}",
  "\\begin{tabular}{lrrrrrr}",
  "\\toprule",
  " & \\multicolumn{3}{c}{$\\Delta$AIC} & \\multicolumn{3}{c}{$\\Delta$ELPD (SE)} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
  "Model & Presence & Biomass & Total & Presence & Biomass & Total \\\\",
  "\\midrule",
  rows1,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
), here("output", "tables", "tab1_model_selection.tex"))

# 02 Table S1: full model with each connectivity metric ----
tab_s1 <- bind_rows(
  filter(readRDS(here("data", "sdm", "main", "model_comparison.rds")), model == "space_env_conn"),
  readRDS(here("data", "sdm", "sensitivity", "model_comparison.rds"))
) |>
  delta_aic() |>
  mutate(label = metric_labels[model]) |>
  arrange(d_total)

rows_s1 <- tab_s1 |>
  transmute(line = paste0(label, " & ", fmt(d_presence), " & ", fmt(d_biomass), " & ", fmt(d_total), " \\\\")) |>
  pull(line)

writeLines(c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{$\\Delta$AIC of the full Space + Environment + Connectivity model with each of the five connectivity metrics, for presence, biomass where present and in total (the best model in each column is shown in bold). All metrics are presence-weighted; all models include survey, year, the environmental predictors and a spatial random field.}",
  "\\label{tab:connectivity-metrics}",
  "\\begin{tabular}{lrrr}",
  "\\toprule",
  " & \\multicolumn{3}{c}{$\\Delta$AIC} \\\\",
  "\\cmidrule(lr){2-4}",
  "Connectivity metric & Presence & Biomass & Total \\\\",
  "\\midrule",
  rows_s1,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
), here("output", "tables", "tabS1_connectivity_metrics.tex"))
