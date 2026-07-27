# Figure 2: standardized coefficients of the best-supported model.
#
# The full Space + Environment + Connectivity model, with the connectivity slot
# rotated over the three metrics, for both responses. Because the three metrics
# give nearly-tied AICs, the point is that the effect sizes agree: each
# predictor's estimates cluster and the connectivity coefficient is negative for
# every metric. Intercept and survey (gear) terms are dropped as nuisance;
# estimates are on the standardized link scale, so effects are comparable within
# a response.

library(tidyverse)
library(here)
library(sdmTMB)
library(rcartocolor)
library(ggstats)

# 01 Load fits ----
fits <- readRDS(here("data", "intermediate", "sdm_fits.rds"))
conn_ids <- names(fits)[grepl("space_env_conn", names(fits))]

# 02 Tidy fixed effects ----
coefs <- map_dfr(conn_ids, function(id) {
  f <- fits[[id]]
  if (is.null(f)) {
    return(NULL)
  }
  tidy(f, effects = "fixed", conf.int = TRUE) |>
    mutate(fit = id)
}) |>
  separate(fit, into = c("response", "model_id"), sep = "_", extra = "merge")

# 03 Labels ----
term_labs <- c(
  depth_std = "Depth",
  temp_std = "Temperature",
  oxy_std = "Oxygen",
  sal_std = "Salinity",
  shear_max_std = "Shear stress",
  eigen_centrality_std = "Connectivity",
  in_strength_std = "Connectivity",
  biomass_supply_log_std = "Connectivity"
)
term_levels <- c("Connectivity", "Shear stress", "Salinity", "Oxygen", "Temperature", "Depth")

conn_labs <- c(
  in_strength = "In-strength",
  biomass_supply_log = "Biomass supply (log)",
  eigen = "Eigenvector centrality"
)
conn_levels <- c("In-strength", "Biomass supply (log)", "Eigenvector centrality")

coefs <- coefs |>
  filter(term %in% names(term_labs)) |>
  mutate(
    term = factor(term_labs[term], levels = term_levels),
    conn = str_extract(model_id, "(eigen|in_strength|biomass_supply_log)$"),
    metric = factor(conn_labs[conn], levels = conn_levels),
    response = factor(response, levels = c("present", "biomass"), labels = c("Presence", "Biomass"))
  )

# 04 Plot ----
dodge <- position_dodge(width = 0.55)
p <- ggplot(coefs, aes(estimate, term, colour = metric)) +
  geom_stripped_rows(colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
    height = 0, position = dodge, linewidth = 0.5
  ) +
  geom_point(position = dodge, size = 2.4) +
  facet_wrap(~response, scales = "free_x") +
  scale_colour_carto_d(palette = "Vivid", name = "Connectivity metric") +
  labs(x = "Standardized coefficient (link scale)", y = NULL) +
  theme_light(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "grey20", face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# 05 Save ----
dir.create(here("output"), showWarnings = FALSE)
ggsave(here("output", "fig2_coefficients.png"), p, width = 10, height = 6, dpi = 600, bg = "white")
