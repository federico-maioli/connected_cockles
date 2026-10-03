# Figure 3: the full delta-gamma model (space + environment + in-strength
# connectivity), its two components side by side:
#   presence                  binomial (logit) - where cockles occur
#   biomass where present     gamma (log)      - how much, where they occur
# Panel a: standardized coefficient of every predictor in each component.
# Panel b: variance partitioning of each component (Nakagawa R2, split
# HMSC-style among Environment/Connectivity - see 07_fit_sdm.R and
# nakagawa_sdmtmb() in R/helpers.R); the biomass component is partitioned over
# the positive observations only.
# The same model without the spatial field is compared in fig_supp.R (Figure S2).
#
# Reads data/sdm/main/space_env_conn(.rds/_variance.rds), fitted in 07_fit_sdm.R.

library(tidyverse)
library(here)
library(sdmTMB)
library(ggstats)
library(patchwork)

part_labs <- c(presence = "Presence (logit link)", biomass = "Biomass where present (log link)")
part_point_cols <- c("Presence (logit link)" = "#5D3A9B", "Biomass where present (log link)" = "#E66100")

term_labs <- c(
  depth_std = "Depth", `I(depth_std^2)` = "Depth\u00b2", temp_std = "Temperature", oxy_std = "Oxygen",
  sal_std = "Salinity", shear_max_std = "Shear stress",
  conn_in_strength_std = "In-strength"
)
term_levels <- c("In-strength", "Shear stress", "Salinity", "Oxygen", "Temperature", "Depth\u00b2", "Depth")

# Okabe-Ito, colourblind-safe; unexplained (distribution) variance is not drawn -
# it is the empty remainder of each bar up to 100%
var_cols <- c("Spatial field" = "#0072B2", "Environment" = "#009E73", "In-strength" = "#E69F00", "Survey" = "#999999", "Year" = "#CC79A7")
label_cols <- c("Spatial field" = "white", "Environment" = "white", "In-strength" = "grey15", "Survey" = "grey15", "Year" = "grey15")

fit <- readRDS(here("data", "sdm", "main", "space_env_conn.rds"))

# 01 Panel a: coefficients of both components ----
coefs <- map(1:2, \(m) tidy(fit, effects = "fixed", model = m, conf.int = TRUE) |> mutate(component = m)) |>
  list_rbind() |>
  filter(term %in% names(term_labs)) |>
  mutate(
    term = factor(term_labs[term], levels = term_levels),
    part = factor(unname(part_labs)[component], levels = unname(part_labs))
  )

p_coeff <- ggplot(coefs, aes(estimate, term, colour = part)) +
  geom_stripped_rows(colour = NA) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0, linewidth = 0.9, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = part_point_cols, name = NULL) +
  labs(x = "Standardized coefficient", y = NULL) +
  theme_light(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# 02 Panel b: variance partitioning of both components ----
# segments >= 5% are labelled in place, thinner ones (e.g. in-strength) get
# their value just above the slice
variance <- readRDS(here("data", "sdm", "main", "space_env_conn_variance.rds")) |>
  filter(share > 0, component != "distribution") |>
  mutate(
    # 07_fit_sdm.R names the block "Connectivity"; this figure is in-strength only
    component = factor(recode(component, spatial = "Spatial field", Connectivity = "In-strength"), levels = names(var_cols)),
    # presence first, at the top of the bar chart (ggplot draws the first
    # discrete-axis level at the bottom, so the level order is reversed)
    group = factor(part_labs[part], levels = rev(unname(part_labs)))
  )

variance_labs <- variance |>
  arrange(group, component) |>
  mutate(x_mid = cumsum(share) - share / 2, .by = group)
narrow_labs <- variance_labs |>
  filter(share < 0.05) |>
  summarise(label = paste(scales::percent(share, accuracy = 0.1), collapse = "\n"), x_mid = mean(x_mid), .by = group)

p_variance <- ggplot(variance, aes(x = share, y = group, fill = component)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
  geom_text(
    data = filter(variance_labs, share >= 0.05),
    aes(x = x_mid, y = group, label = scales::percent(share, accuracy = 0.1), colour = component),
    inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  geom_text(
    data = narrow_labs, aes(x = x_mid, y = group, label = label),
    inherit.aes = FALSE, size = 2.6, colour = "grey25", lineheight = 0.85,
    position = position_nudge(y = 0.5)
  ) +
  scale_fill_manual(values = var_cols, name = NULL, drop = TRUE) +
  scale_colour_manual(values = label_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent, expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.45, 0.75))) +
  labs(x = expression("Share of total variance (Nakagawa " * R^2 * ")"), y = NULL) +
  theme_light(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))

# 03 Assemble and save ----
fig3 <- p_coeff / p_variance +
  plot_layout(heights = c(1.8, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")

dir.create(here("output", "figs"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "fig3_coeff.png"), fig3, width = 7.5, height = 7, dpi = 600, bg = "white")
