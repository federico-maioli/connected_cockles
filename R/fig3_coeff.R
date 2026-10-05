# Figure 3: the full delta-gamma model (space + environment + in-strength
# connectivity), its two components side by side:
#   presence                  binomial (logit) - where cockles occur
#   biomass where present     gamma (log)      - how much, where they occur
# Panel a: standardized coefficient of every predictor in each component.
# Panel b: partial effects - predicted presence probability and biomass where
# present along each predictor, with the others at their mean.
#
# Reads data/sdm/main/space_env_conn.rds (07_fit_sdm.R).

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

# 02 Panel b: partial effects ----
# predicted presence probability (left axis) and log biomass where present (right
# axis, log scale) along each predictor, from the population-level part of
# the full model (re_form = NA, no spatial field), with all other predictors
# at their mean (0, standardised) and for the grab survey in 2021, the year
# with the most grab stations. Each predictor spans the 1st-99th percentile of
# the model data; depth includes its quadratic term. Ribbons are 95% CIs
pe_year <- "2021"
pe_vars <- c(
  depth = "Depth (m)", temp = "Temperature (°C)", sal = "Salinity (psu)",
  oxy = "Oxygen (days < 4 mg/l)", shear_max = "Max. shear stress (N/m²)",
  conn_in_strength = "In-strength"
)
std_cols <- paste0(names(pe_vars), "_std")

partial_effect <- function(var) {
  raw <- fit$data[[var]]
  std <- fit$data[[paste0(var, "_std")]]
  s <- sd(raw) / sd(std)
  m <- mean(raw) - s * mean(std)
  nd <- tibble(value = seq(quantile(raw, 0.01), quantile(raw, 0.99), length.out = 100))
  for (col in std_cols) nd[[col]] <- 0
  nd[[paste0(var, "_std")]] <- (nd$value - m) / s
  nd <- mutate(
    nd,
    survey = factor("Stock", levels = levels(fit$data$survey)),
    year = factor(pe_year, levels = levels(fit$data$year)),
    x_utm = mean(fit$data$x_utm), y_utm = mean(fit$data$y_utm)
  )
  presence <- predict(fit, newdata = nd, re_form = NA, se_fit = TRUE, model = 1)
  biomass <- predict(fit, newdata = nd, re_form = NA, se_fit = TRUE, model = 2)
  tibble(
    value = nd$value,
    p = plogis(presence$est), p_low = plogis(presence$est - 1.96 * presence$est_se), p_high = plogis(presence$est + 1.96 * presence$est_se),
    b = biomass$est, b_low = biomass$est - 1.96 * biomass$est_se, b_high = biomass$est + 1.96 * biomass$est_se
  )
}

# one plot per predictor: log biomass (the link scale) is rescaled onto the
# probability axis and labelled on a secondary axis
pe_plot <- function(var) {
  d <- partial_effect(var)
  p_max <- max(d$p_high)
  b_range <- range(c(d$b_low, d$b_high))
  to_p <- function(b) (b - b_range[1]) / diff(b_range) * p_max
  ggplot(d, aes(value)) +
    geom_ribbon(aes(ymin = p_low, ymax = p_high), fill = part_point_cols[[1]], alpha = 0.2) +
    geom_line(aes(y = p), colour = part_point_cols[[1]], linewidth = 0.8) +
    geom_ribbon(aes(ymin = to_p(b_low), ymax = to_p(b_high)), fill = part_point_cols[[2]], alpha = 0.2) +
    geom_line(aes(y = to_p(b)), colour = part_point_cols[[2]], linewidth = 0.8) +
    scale_y_continuous(
      name = "Presence probability",
      sec.axis = sec_axis(
        \(y) b_range[1] + y / p_max * diff(b_range),
        name = "Log biomass where present"
      )
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = pe_vars[[var]]) +
    theme_light(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title.y.left = element_text(colour = part_point_cols[[1]]),
      axis.title.y.right = element_text(colour = part_point_cols[[2]]),
      axis.text.y.left = element_text(colour = part_point_cols[[1]]),
      axis.text.y.right = element_text(colour = part_point_cols[[2]])
    )
}

# the y-axis titles are identical in every plot, so patchwork shows each once;
# the panel tag sits on the first plot
pe_plots <- map(names(pe_vars), pe_plot)
pe_plots[[1]] <- pe_plots[[1]] + labs(tag = "b)")
p_partial <- wrap_plots(pe_plots, ncol = 3) +
  plot_layout(axis_titles = "collect")

# 03 Assemble and save ----
# nested (not wrapped) so the plot areas of a and b line up; the axis title on
# the left of b is freed from that alignment so it sits next to its plots
# rather than out at the edge of panel a's wide labels
fig3 <- (p_coeff + labs(tag = "a)")) / free(p_partial, type = "label", side = "l") +
  plot_layout(heights = c(1, 1.1)) &
  theme(plot.tag = element_text(size = 12, face = "bold", colour = "grey20"))

dir.create(here("output", "figs", "main"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "main", "fig3_coeff.png"), fig3, width = 7.5, height = 8, dpi = 600, bg = "white")
