# Supplementary figure: is the connectivity effect confounded with space?
#
# The same Environment + Connectivity model is fitted with and without a spatial
# random field (space_env_conn vs env_conn). If connectivity were merely standing
# in for broad spatial structure, dropping the field would inflate its
# coefficient; if the two estimates agree, the connectivity effect is not simply
# spatial autocorrelation in disguise.
#
# Both responses are shown; estimates are on the standardized link scale, so the
# with / without pair is directly comparable within a response.

library(tidyverse)
library(here)
library(sdmTMB)
library(ggstats)

# tidy the fixed effects of one structure for one response, tagged by metric
coef_data <- function(fits, resp, structure, ok_ids) {
  prefix <- paste0(resp, "_", structure, "_")
  ids <- names(fits)[startsWith(names(fits), prefix)]
  ids <- intersect(ids, ok_ids)
  map(ids, function(id) {
    f <- fits[[id]]
    if (is.null(f)) {
      return(NULL)
    }
    tidy(f, effects = "fixed", conf.int = TRUE) |>
      mutate(conn = str_remove(id, paste0("^", prefix)), structure = structure)
  }) |>
    list_rbind()
}

# 01 Load fits ----
# non-converged fits are dropped: their coefficients are not trustworthy
fits <- readRDS(here("data", "intermediate", "sdm_fits.rds"))
ok_ids <- readRDS(here("data", "final", "sdm_model_comparison.rds")) |>
  filter(converged) |>
  transmute(id = paste(response, model, sep = "_")) |>
  pull(id)

# 02 Labels ----
conn_labs <- c(
  log_biomass_in_strength = "Biomass in-strength (log)",
  deg_in = "In-degree",
  in_strength = "In-strength",
  eigen_centrality = "Eigenvector centrality",
  closeness_centrality = "Closeness centrality"
)
conn_terms <- paste0(names(conn_labs), "_std")
structure_labs <- c(space_env_conn = "With spatial field", env_conn = "Without spatial field")

# 03 Connectivity coefficient, with and without the field ----
# keep only the connectivity term of each fit - that is the one at issue
coefs <- c("biomass", "present") |>
  set_names() |>
  map(function(resp) {
    bind_rows(
      coef_data(fits, resp, "space_env_conn", ok_ids),
      coef_data(fits, resp, "env_conn", ok_ids)
    )
  }) |>
  list_rbind(names_to = "response") |>
  filter(term %in% conn_terms) |>
  mutate(
    metric = factor(conn_labs[conn], levels = rev(unname(conn_labs))),
    structure = factor(structure_labs[structure], levels = unname(structure_labs)),
    response = factor(response,
      levels = c("present", "biomass"),
      labels = c("Presence (logit)", "Biomass (log)")
    )
  )

# 04 Plot ----
dodge <- position_dodge(width = 0.55)
p <- ggplot(coefs, aes(estimate, metric, colour = structure)) +
  geom_stripped_rows(colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
    orientation = "y", width = 0, position = dodge, linewidth = 0.5
  ) +
  geom_point(position = dodge, size = 2.4) +
  facet_wrap(~response, scales = "free_x") +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
  labs(
    x = "Standardized connectivity coefficient", y = NULL,
    caption = "A missing point means that model did not converge."
  ) +
  theme_light(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "grey20", face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.caption = element_text(colour = "grey40", hjust = 0)
  )

# 05 Save ----
dir.create(here("output"), showWarnings = FALSE)
ggsave(here("output", "figS_space_confounding.png"), p,
  width = 9, height = 4.5, dpi = 600, bg = "white"
)

# 06 Report the shift ----
# how much does each connectivity coefficient move when the field is dropped?
coefs |>
  select(response, metric, structure, estimate) |>
  pivot_wider(names_from = structure, values_from = estimate) |>
  rename(with_field = `With spatial field`, without_field = `Without spatial field`) |>
  mutate(
    across(c(with_field, without_field), \(x) round(x, 3)),
    ratio = round(without_field / with_field, 2)
  ) |>
  arrange(response, metric) |>
  print(n = Inf)
