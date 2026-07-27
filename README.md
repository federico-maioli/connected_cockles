# connected_cockles

Integrating larval connectivity from a biophysical dispersal model into species
distribution models for cockles (*Cerastoderma edule*) in the Limfjord, Denmark.

Data are kept local (see `.gitignore`); this repository tracks the analysis
scripts, figures, and tables only.

## Pipeline (`R/`)

| Script | Purpose |
|--------|---------|
| `00_clean_cockles.R` | Clean the cockle survey data, attach depth |
| `01_compute_flow_matrix.R` | Flow (settlement-probability) matrix from the dispersal model |
| `02_make_grid.R` | Build the spatial grid |
| `03_predict_avg_biomass.R` | Predict mean biomass / presence to the grid (depth-only) |
| `04_weight_flow_matrix.R` | Biomass-weighted supply and structural graph metrics |
| `05_extract_connectivity.R` | Attach connectivity metrics to survey observations |
| `06_fit_sdm.R` | Fit and compare the spatial SDMs (presence, biomass) |
| `07_variance_partitioning.R` | Nakagawa marginal / conditional variance partition |
| `fig1_study_area.R` | Figure 1 — study area and connectivity clusters |
| `fig2_coefficients.R` | Figure 2 — standardized coefficients |
| `fig_covariates.R` | Covariate maps |
| `tab_model_selection.R` | Model-selection table (LaTeX) |

## Outputs (`output/`)

Figures (`.png`), tables (`tables/*.tex`), and methods text (`methods_*.tex`).

## Models

Spatial GLMMs (`sdmTMB`) with a coastline barrier mesh: biomass (Tweedie) and
presence (Bernoulli), with environmental, gear, and connectivity covariates and
a spatial random field.
