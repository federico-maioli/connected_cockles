# Variance partitioning for the cockle biomass (Tweedie) SDM.
#
# Nakagawa & Schielzeth (2013) marginal / conditional R2 for the full model
# (survey + environment + connectivity + spatial field), all on the link scale:
#   var_f  variance of the fixed-effects linear predictor
#   var_l  variance of the spatial random field
#   var_d  Tweedie observation-level (distribution) variance
# marginal R2    = var_f / (var_f + var_l + var_d)     fixed effects
# conditional R2 = (var_f + var_l) / (var_f + var_l + var_d)  fixed + field
# Reported as a table splitting total variance into: spatial field, environment,
# survey, connectivity, unexplained. The fixed (marginal) share is divided among
# the three fixed blocks by each block's covariance with the fixed predictor, so
# the blocks sum exactly to the marginal R2.

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
library(rnaturalearth)

# compatibility shim: sdmTMBextra 0.0.5 calls fmesher::fm_identical_CRS, renamed
# to fm_crs_is_identical in fmesher >= 0.7; register the old name as an alias
local({
  exps <- .getNamespaceInfo(asNamespace("fmesher"), "exports")
  if (!exists("fm_identical_CRS", envir = exps, inherits = FALSE)) {
    assign("fm_identical_CRS", "fm_crs_is_identical", envir = exps)
  }
})

# 01 Load and standardise ----
std_vars <- c("depth", "temp", "oxy", "sal", "shear_max", "in_strength")
dat <- readRDS(here("data", "final", "cockles_connectivity.rds")) |>
  mutate(survey = factor(survey)) |>
  mutate(across(all_of(std_vars), ~ as.numeric(scale(.x)), .names = "{.col}_std")) |>
  filter(
    !is.na(biomass), !is.na(present), !is.na(x_utm), !is.na(y_utm),
    if_all(all_of(paste0(std_vars, "_std")), ~ !is.na(.x))
  )

# 02 Coastline barrier mesh ----
land <- ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf") |>
  st_make_valid()
region <- st_bbox(st_transform(
  st_as_sf(dat |> transmute(x = x_utm * 1000, y = y_utm * 1000), coords = c("x", "y"), crs = 32632),
  4326
))
region["xmin"] <- region["xmin"] - 0.4
region["ymin"] <- region["ymin"] - 0.4
region["xmax"] <- region["xmax"] + 0.4
region["ymax"] <- region["ymax"] + 0.4
land_region <- suppressWarnings(st_crop(land, region)) |> st_transform(32632)

mesh <- make_mesh(dat, c("x_utm", "y_utm"), cutoff = 1)
barrier_mesh <- sdmTMBextra::add_barrier_mesh(
  mesh, land_region,
  range_fraction = 0.1, proj_scaling = 1000, plot = FALSE
)

# 03 Fit the full biomass model ----
fixed_form <- ~ survey + depth_std + temp_std + oxy_std + sal_std + shear_max_std + in_strength_std
fit <- sdmTMB(
  update(fixed_form, biomass ~ .),
  data = dat, mesh = barrier_mesh, spatial = "on", family = tweedie(link = "log")
)

# 04 Nakagawa marginal / conditional R2 ----
pr <- predict(fit)
var_f <- var(pr$est_non_rf) # fixed-effects linear predictor
var_l <- var(pr$est_rf) # spatial random field
rp <- tidy(fit, effects = "ran_pars")
phi <- rp$estimate[rp$term == "phi"]
p_tw <- rp$estimate[rp$term == "tweedie_p"]
var_d <- mean(log1p(phi * exp(pr$est)^(p_tw - 2))) # Tweedie observation variance

tot <- var_f + var_l + var_d
marginal_r2 <- var_f / tot
conditional_r2 <- (var_f + var_l) / tot
cat(sprintf(
  "marginal (fixed) R2 = %.3f, conditional (fixed + field) R2 = %.3f\n",
  marginal_r2, conditional_r2
))

# 05 Split the fixed share between environment, survey and connectivity ----
# each fixed block's share = cov(block, fixed predictor) / tot; the three sum
# exactly to the marginal (fixed) R2 (block covariances absorbed by this split)
coefs <- tidy(fit, effects = "fixed") |>
  select(term, estimate) |>
  deframe()
X <- model.matrix(fixed_form, data = dat)
lp_term <- sweep(X, 2, coefs[colnames(X)], "*")
env_cols <- c("depth_std", "temp_std", "oxy_std", "sal_std", "shear_max_std")
b_env <- rowSums(lp_term[, env_cols, drop = FALSE])
b_survey <- lp_term[, "surveyStock"]
b_conn <- lp_term[, "in_strength_std"]
eta_fixed <- b_env + b_survey + b_conn

partition <- tibble(
  component = c("Spatial field", "Environment", "Survey (gear)", "Connectivity", "Unexplained"),
  share = c(
    var_l / tot,
    cov(b_env, eta_fixed) / tot,
    cov(b_survey, eta_fixed) / tot,
    cov(b_conn, eta_fixed) / tot,
    var_d / tot
  )
)
print(partition)
cat(sprintf("sum of shares = %.3f\n", sum(partition$share)))

# 06 LaTeX table ----
tab <- partition |> mutate(pct = sprintf("%.1f", 100 * share))
rows <- paste0(tab$component, " & ", tab$pct, " \\\\")

latex <- c(
  "\\begin{table}[ht]",
  "\\centering",
  sprintf("\\caption{Variance partitioning of cockle biomass from the spatial Tweedie SDM (Nakagawa \\& Schielzeth marginal / conditional $R^2$). The spatial field and the three fixed-effect blocks together give the conditional $R^2$ (%.1f\\%%); the fixed blocks alone give the marginal $R^2$ (%.1f\\%%). Shares are percentages of total variance and sum to 100\\%%.}", 100 * conditional_r2, 100 * marginal_r2),
  "\\label{tab:variance-partitioning}",
  "\\begin{tabular}{lr}",
  "\\toprule",
  "Component & Variance explained (\\%) \\\\",
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

# 07 Plot (optional companion to the table) ----
part_cols <- c(
  "Spatial field" = "#4C72B0", "Environment" = "#2E8B57",
  "Survey (gear)" = "#7F7F7F", "Connectivity" = "#E58606", "Unexplained" = "#ECECEC"
)
plot_dat <- partition |>
  mutate(component = factor(component, levels = names(part_cols)))
comp_labs <- plot_dat |>
  mutate(lab = paste0(component, " (", scales::percent(share, accuracy = 0.1), ")")) |>
  select(component, lab) |>
  deframe()

p_var <- ggplot(plot_dat, aes(x = share, y = "Biomass", fill = component)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = if_else(share >= 0.04, scales::percent(share, accuracy = 0.1), "")),
    position = position_stack(vjust = 0.5), size = 3.2, colour = "grey10"
  ) +
  scale_fill_manual(values = part_cols, labels = comp_labs, name = NULL) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
  labs(x = expression("Share of total variance (Nakagawa " * R^2 * ")"), y = NULL) +
  theme_light(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 1))

dir.create(here("output"), showWarnings = FALSE)
ggsave(here("output", "fig3_variance_partitioning.png"), p_var, width = 9, height = 3, dpi = 600, bg = "white")
