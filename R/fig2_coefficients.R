# Figure 2: standardized coefficients of the best-supported model (biomass).
#
# The full Space + Environment + Connectivity model, with the connectivity slot
# rotated over the five metrics. Because the metrics give nearly-tied AICs, the
# point is that the effect sizes agree: each predictor's estimates cluster and
# the connectivity coefficient keeps its sign across metrics. Intercept and
# survey (gear) terms are dropped as nuisance; estimates are on the standardized
# link scale, so effects are comparable within a response.
#
# Writes two files from the same builder:
#   fig2_coefficients.png            biomass (Tweedie)  - main text
#   figS_coefficients_presence.png   presence (binomial) - supplementary

library(tidyverse)
library(here)
library(sdmTMB)
library(rcartocolor)
library(ggstats)

# tidy the fixed effects of every converged space_env_conn fit for one response
coef_data <- function(fits, resp, ok_ids) {
  ids <- names(fits)[startsWith(names(fits), paste0(resp, "_space_env_conn_"))]
  ids <- intersect(ids, ok_ids)
  map(ids, function(id) {
    f <- fits[[id]]
    if (is.null(f)) {
      return(NULL)
    }
    tidy(f, effects = "fixed", conf.int = TRUE) |>
      mutate(conn = str_remove(id, paste0("^", resp, "_space_env_conn_")))
  }) |>
    list_rbind()
}

coef_plot <- function(d, xlab) {
  dodge <- position_dodge(width = 0.6)
  ggplot(d, aes(estimate, term, colour = metric)) +
    geom_stripped_rows(colour = NA) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
      orientation = "y", width = 0, position = dodge, linewidth = 0.5
    ) +
    geom_point(position = dodge, size = 2.4) +
    scale_colour_carto_d(palette = "Vivid", name = "Connectivity metric") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(x = xlab, y = NULL) +
    theme_light(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

# 01 Load fits ----
# non-converged fits are dropped: their coefficients are not trustworthy
fits <- readRDS(here("data", "derived", "sdm_fits.rds"))
ok_ids <- readRDS(here("data", "derived", "sdm_model_comparison.rds")) |>
  filter(converged) |>
  transmute(id = paste(response, model, sep = "_")) |>
  pull(id)

# 02 Labels ----
# every connectivity metric shares one "Connectivity" row, so the metrics line
# up against each other in the same slot
term_labs <- c(
  depth_std = "Depth",
  temp_std = "Temperature",
  oxy_std = "Oxygen",
  sal_std = "Salinity",
  shear_max_std = "Shear stress",
  log_biomass_in_strength_std = "Connectivity",
  deg_in_std = "Connectivity",
  in_strength_std = "Connectivity",
  eigen_centrality_std = "Connectivity",
  closeness_centrality_std = "Connectivity"
)
term_levels <- c("Connectivity", "Shear stress", "Salinity", "Oxygen", "Temperature", "Depth")

conn_labs <- c(
  log_biomass_in_strength = "Biomass in-strength (log)",
  deg_in = "In-degree",
  in_strength = "In-strength",
  eigen_centrality = "Eigenvector centrality",
  closeness_centrality = "Closeness centrality"
)
conn_levels <- unname(conn_labs)

prep <- function(resp) {
  coef_data(fits, resp, ok_ids) |>
    filter(term %in% names(term_labs)) |>
    mutate(
      term = factor(term_labs[term], levels = term_levels),
      metric = factor(conn_labs[conn], levels = conn_levels)
    )
}

# 03 Biomass (main text) ----
p_biomass <- prep("biomass") |>
  coef_plot("Standardized coefficient (log link)")

# 04 Presence (supplementary) ----
p_present <- prep("present") |>
  coef_plot("Standardized coefficient (logit link)")

# 05 Save ----
dir.create(here("output"), showWarnings = FALSE)
ggsave(here("output", "fig2_coefficients.png"), p_biomass,
  width = 8, height = 5.4, dpi = 600, bg = "white"
)
ggsave(here("output", "figS_coefficients_presence.png"), p_present,
  width = 8, height = 5.4, dpi = 600, bg = "white"
)
