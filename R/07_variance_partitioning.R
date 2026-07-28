# Variance partitioning for the cockle biomass (Tweedie) SDMs.
#
# Nakagawa & Schielzeth (2013) marginal / conditional R2, all on the link scale:
#   var_f  variance of the fixed-effects linear predictor
#   var_l  variance of the spatial random field (0 when the field is off)
#   var_d  Tweedie observation-level (distribution) variance
# marginal R2    = var_f / (var_f + var_l + var_d)             fixed effects
# conditional R2 = (var_f + var_l) / (var_f + var_l + var_d)   fixed + field
# Total variance is split into: spatial field, environment, survey, connectivity,
# unexplained. The fixed (marginal) share is divided among the three fixed blocks
# by each block's covariance with the fixed predictor, so the blocks sum exactly
# to the marginal R2.
#
# Run for every connectivity metric, twice:
#   space_env_conn  with the spatial field    -> fig3 (main text) + LaTeX table
#   env_conn        without the spatial field -> supplementary sensitivity
# Without the field, whatever the field was absorbing has to go somewhere, so the
# comparison shows how much of each block's share depends on it being there.
#
# Fits are reused from 06 so the partitioning matches the AIC table and the
# coefficient figures exactly (same data, same mesh, same estimates).

library(tidyverse)
library(here)
library(sdmTMB)

env_cols <- c("depth_std", "temp_std", "oxy_std", "sal_std", "shear_max_std")

# split one fitted model's total variance into components
partition_fit <- function(fit, conn_term) {
  # with a field, predict() splits est into est_non_rf + est_rf; with the field
  # off it returns est only, which is already the fixed-effects predictor
  pr <- predict(fit)
  var_f <- if (is.null(pr$est_non_rf)) var(pr$est) else var(pr$est_non_rf)
  var_l <- if (is.null(pr$est_rf)) 0 else var(pr$est_rf)
  rp <- tidy(fit, effects = "ran_pars")
  phi <- rp$estimate[rp$term == "phi"]
  p_tw <- rp$estimate[rp$term == "tweedie_p"]
  var_d <- mean(log1p(phi * exp(pr$est)^(p_tw - 2))) # Tweedie observation variance
  tot <- var_f + var_l + var_d

  # each fixed block's share = cov(block, fixed predictor) / tot; the three sum
  # exactly to the marginal (fixed) R2 (block covariances absorbed by this split)
  coefs <- tidy(fit, effects = "fixed") |>
    select(term, estimate) |>
    deframe()
  X <- model.matrix(delete.response(terms(formula(fit))), data = fit$data)
  lp_term <- sweep(X, 2, coefs[colnames(X)], "*")
  b_env <- rowSums(lp_term[, env_cols, drop = FALSE])
  b_survey <- lp_term[, "surveyStock"]
  b_conn <- lp_term[, conn_term]
  eta_fixed <- b_env + b_survey + b_conn

  tibble(
    component = c("Spatial field", "Environment", "Survey (gear)", "Connectivity", "Unexplained"),
    share = c(
      var_l / tot,
      cov(b_env, eta_fixed) / tot,
      cov(b_survey, eta_fixed) / tot,
      cov(b_conn, eta_fixed) / tot,
      var_d / tot
    ),
    marginal_r2 = var_f / tot,
    conditional_r2 = (var_f + var_l) / tot
  )
}

# 01 Load fits ----
# only converged models are partitioned
fits <- readRDS(here("data", "derived", "sdm_fits.rds"))
ok_ids <- readRDS(here("data", "derived", "sdm_model_comparison.rds")) |>
  filter(converged) |>
  transmute(id = paste(response, model, sep = "_")) |>
  pull(id)

metric_labs <- c(
  log_biomass_in_strength = "Biomass in-strength (log)",
  deg_in = "In-degree",
  in_strength = "In-strength",
  eigen_centrality = "Eigenvector centrality",
  closeness_centrality = "Closeness centrality"
)

# 02 Partition every metric, with and without the spatial field ----
grid <- expand_grid(
  metric = names(metric_labs),
  structure = c("space_env_conn", "env_conn")
) |>
  mutate(id = paste("biomass", structure, metric, sep = "_")) |>
  filter(id %in% ok_ids)

parts <- grid |>
  mutate(part = map2(id, metric, \(i, m) partition_fit(fits[[i]], paste0(m, "_std")))) |>
  unnest(part) |>
  mutate(
    metric_lab = factor(metric_labs[metric], levels = rev(unname(metric_labs))),
    component = factor(component, levels = c(
      "Spatial field", "Environment", "Survey (gear)", "Connectivity", "Unexplained"
    ))
  )

# 03 Report ----
parts |>
  filter(component == "Connectivity") |>
  select(metric_lab, structure, connectivity_share = share, marginal_r2, conditional_r2) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(metric_lab, structure) |>
  print(n = Inf)

# 04 Plot ----
part_cols <- c(
  "Spatial field" = "#4C72B0", "Environment" = "#2E8B57",
  "Survey (gear)" = "#7F7F7F", "Connectivity" = "#E58606", "Unexplained" = "#ECECEC"
)

partition_plot <- function(d) {
  # drop components that are structurally absent (no field when spatial = off)
  d <- d |> filter(share > 0)
  ggplot(d, aes(x = share, y = metric_lab, fill = component)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = if_else(share >= 0.05, scales::percent(share, accuracy = 0.1), "")),
      position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.9, colour = "grey10"
    ) +
    scale_fill_manual(values = part_cols, name = NULL, drop = TRUE) +
    scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
    labs(x = expression("Share of total variance (Nakagawa " * R^2 * ")"), y = NULL) +
    theme_light(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 1))
}

p_space <- partition_plot(filter(parts, structure == "space_env_conn"))
p_nospace <- partition_plot(filter(parts, structure == "env_conn"))

dir.create(here("output"), showWarnings = FALSE)
ggsave(here("output", "fig3_variance_partitioning.png"), p_space,
  width = 9.5, height = 4.2, dpi = 600, bg = "white"
)
ggsave(here("output", "figS_variance_partitioning_nospace.png"), p_nospace,
  width = 9.5, height = 4.2, dpi = 600, bg = "white"
)

# 05 LaTeX table (main model, with the spatial field) ----
# components as rows, connectivity metrics as columns
wide <- parts |>
  filter(structure == "space_env_conn") |>
  mutate(pct = sprintf("%.1f", 100 * share)) |>
  select(component, metric, pct) |>
  pivot_wider(names_from = metric, values_from = pct) |>
  arrange(component)

metric_order <- names(metric_labs)
header <- paste0("Component & ", paste(metric_labs[metric_order], collapse = " & "), " \\\\")
rows <- wide |>
  select(component, all_of(metric_order)) |>
  as.matrix() |>
  apply(1, \(r) paste0(paste(r, collapse = " & "), " \\\\"))

r2 <- parts |>
  filter(structure == "space_env_conn", component == "Connectivity") |>
  summarise(cond = mean(conditional_r2), marg = mean(marginal_r2))

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

dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
writeLines(latex, here("output", "tables", "variance_partitioning.tex"))
