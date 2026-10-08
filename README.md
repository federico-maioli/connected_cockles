# connected_cockles

Does larval connectivity from a biophysical dispersal model improve species
distribution models for cockles (*Cerastoderma* spp.) in the Limfjorden,
Denmark, beyond the environment and space?

Data are kept local (see `.gitignore`); this repository tracks the analysis
scripts, figures, tables, and the manuscript only. Exploratory analyses live in
`R/exploration/`, which is also not tracked.

## Layout

```
R/
  01_..08_     analysis, run in order; write only to data/
  fig*/tab*    outputs; write only to output/
  helpers.R    fmesher compatibility shim
data/          local only, never tracked
output/
  figs/main/   fig1..fig3 (main text)
  figs/supp/   figS1..figS9 (supplementary)
  tables/      tab1, tabS1 (LaTeX, \input directly by the manuscript)
manuscript/    main.tex, supp.tex (Supporting Information), bibliography, and style files
```

## Analysis - run in order

| Script | Purpose |
|--------|---------|
| `01_prepare_grid.R` | 2 x 2 km connectivity grid with environmental covariates |
| `02_clean_cockles.R` | Cockle surveys (Stock 2018 suction dredge, Stock/HighRes 2021-2025 grab, KSKV dredge 2018-2023) with environmental covariates |
| `03_build_mesh.R` | Variable-resolution SPDE mesh with the coastline built in and a land barrier |
| `04_fit_suitability.R` | Presence suitability (survey + year + spatial field) predicted onto the grid |
| `05_weight_connectivity.R` | Presence-weighted connectivity metrics (in-degree, in-strength, in-closeness, eigenvector centrality, transitivity) |
| `06_match_connectivity.R` | Attach the metrics to the survey stations and the grid |
| `07_fit_sdm.R` | Delta-gamma SDMs: space, + environment, + in-strength, + both; the full model with each other metric |
| `08_cross_validation.R` | Spatially blocked 10-fold cross-validation of the four main models (ELPD) |

## Outputs - run once the analysis has run

| Script | Writes |
|--------|--------|
| `fig1_map.R` | `figs/main/fig1_map.png` - survey stations and settlement footprints |
| `fig2_predictors.R` | `figs/main/fig2_predictors.png` - environmental predictors and in-strength |
| `fig3_coeff.R` | `figs/main/fig3_coeff.png` - coefficients and partial effects |
| `fig_supp.R` | `figs/supp/figS1..figS8` (all but S7), numbered in order of citation in main.tex |
| `figS7_predictions.R` | `figs/supp/figS7_predictions.png` - biomass predicted by the full model (Space + Environment + Connectivity) and its uncertainty |
| `tab1_model_selection.R` | `tables/tab1_model_selection.tex`, `tables/tabS1_connectivity_metrics.tex` |

## Manuscript

`manuscript/main.tex` pulls figures and tables directly from `output/` via
`\graphicspath` and `\input`, so they cannot drift from the fitted models.

## Package versions

- `sdmTMB` >= 1.1.0.9024 (development version, `remotes::install_github("pbs-assess/sdmTMB")`)
  is needed for `predictive = "mle-mvn"` and `loo::elpd()` in `08_cross_validation.R`.
- `sdmTMBextra` 0.0.5 calls `fmesher::fm_identical_CRS`, renamed in fmesher >= 0.7;
  `helpers.R` registers the old name as an alias. Remove it once `sdmTMBextra` is updated.
