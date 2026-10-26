# Figure S8: is suitable habitat that holds no cockles out of reach of larvae?
# The sampled 2 x 2 km cells are classed as
#   occupied          predicted present by the Space + Environment occurrence
#                     model (07_fit_sdm.R), with its spatial field
#   suitable, empty   predicted absent, but environmentally suitable
#   unsuitable, empty
# Environmental suitability comes from the same occurrence model refitted
# without the spatial field (survey, year, sediment and environment), so that
# only the environment defines it. Both are predicted at the cell centres for a
# grab sample in 2025 and cut where they best separate cells with and without
# recorded cockles (maximum true skill statistic, TSS). Connectivity is left
# out of both, so the classes do not depend on the supply compared across them.
# Supply is presence-weighted in-strength as used in the SDMs (raw MIKE flow
# matrix, each source cell weighted by its predicted presence from
# 04_fit_suitability.R, self-supply excluded; 05_weight_connectivity.R),
# unnormalised. If suitable empty cells receive as much supply as occupied
# cells, larval access does not explain why they are empty.
# (a) map of the classes, (b) supply by class.
#
# Reads data/sdm/main/space_env.rds, data/grid/grid_env.rds,
# data/grid/grid_suitability.rds and data/connectivity/raw.

library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
library(patchwork)

n_col <- 65
n_cell <- 65 * 36
cell_id <- function(x_utm, y_utm) (floor((x_utm * 1000 - 450074) / 2000) + 1) + floor((y_utm * 1000 - 6258093) / 2000) * n_col

# 01 Flow matrix ----
# flow[i, j] = larvae released in i that end in j / larvae released in i
cmn <- read_csv(here("data", "connectivity", "raw", "ABM_cmn_all.csv"), col_names = FALSE, na = "1e-35", show_col_types = FALSE) |>
  as.matrix()
cmn[is.na(cmn)] <- 0
dimnames(cmn) <- NULL
release <- read_csv(here("data", "connectivity", "raw", "ABM_rel_all.csv"), col_names = FALSE, show_col_types = FALSE)[[1]]
flow <- ifelse(matrix(release, n_cell, n_cell) > 0, cmn / release, 0)
diag(flow) <- 0

# 02 Cells and predictors ----
fit <- readRDS(here("data", "sdm", "main", "space_env.rds"))$occurrence
dat <- fit$data |>
  mutate(id = cell_id(x_utm, y_utm))
cells <- dat |>
  summarise(stations = n(), observed = any(present == 1), .by = id)

# grid covariates standardised with the model data's own mean and SD, as in
# figS4_predictions.R
std_vars <- c("depth", "temp", "sal", "oxy", "shear_max")
newdata <- readRDS(here("data", "grid", "grid_env.rds")) |>
  filter(id %in% cells$id) |>
  mutate(
    x_utm = x / 1000, y_utm = y / 1000,
    survey = factor("Stock", levels = levels(dat$survey)),
    year = factor("2025", levels = levels(dat$year))
  )
for (v in std_vars) {
  s <- sd(dat[[v]]) / sd(dat[[paste0(v, "_std")]])
  m <- mean(dat[[v]]) - s * mean(dat[[paste0(v, "_std")]])
  newdata[[paste0(v, "_std")]] <- (newdata[[v]] - m) / s
}
newdata <- filter(newdata, if_all(all_of(c(paste0(std_vars, "_std"), "sediment")), \(x) !is.na(x)))

# 03 Occupancy and environmental suitability ----
fit_env <- sdmTMB(formula(fit$formula[[1]]), data = dat, spatial = "off", family = binomial(link = "logit"))
print(tidy(fit_env, conf.int = TRUE))

cells <- cells |>
  inner_join(
    newdata |>
      select(id, x, y) |>
      mutate(
        p_occ = plogis(predict(fit, newdata = newdata)$est),
        suitability = predict(fit_env, newdata = newdata, type = "response")$est
      ),
    by = "id"
  )

# true skill statistic of a cut-off on score s against observed cell occupancy,
# the cut-off that maximises it, and the AUC
tss <- function(s, k) mean(s[cells$observed] >= k) + mean(s[!cells$observed] < k) - 1
best_cut <- function(s) (\(k) k[which.max(map_dbl(k, \(x) tss(s, x)))])(sort(unique(s)))
auc <- function(s, y) (sum(rank(s)[y]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
cut_occ <- best_cut(cells$p_occ)
cut_suit <- best_cut(cells$suitability)
cat("occupancy: AUC", round(auc(cells$p_occ, cells$observed), 2), "| TSS", round(tss(cells$p_occ, cut_occ), 2), "\n")
cat("suitability: AUC", round(auc(cells$suitability, cells$observed), 2), "| TSS", round(tss(cells$suitability, cut_suit), 2), "\n")

# 04 Supply by class ----
w <- replace_na(readRDS(here("data", "grid", "grid_suitability.rds"))$suit_presence, 0)
cells <- cells |>
  mutate(
    occupied = p_occ >= cut_occ,
    supply = colSums(flow * w)[id],
    class = case_when(
      occupied ~ "Occupied",
      suitability >= cut_suit ~ "Suitable, empty",
      .default = "Unsuitable, empty"
    ) |> factor(levels = c("Occupied", "Suitable, empty", "Unsuitable, empty"))
  )
print(table(observed = cells$observed, predicted = cells$occupied))

# reference: the 10th percentile of supply among occupied cells
supply_ref <- quantile(cells$supply[cells$occupied], 0.1)
summarise_cells <- function(d) {
  d |>
    summarise(
      cells = n(),
      reached = round(100 * mean(supply > 0), 1),
      as_supplied_as_beds = round(100 * mean(supply >= supply_ref), 1),
      median_supply = signif(median(supply), 3),
      .by = class
    ) |>
    arrange(class)
}
print(summarise_cells(cells))
# empty cells with few stations are weak absences
print(summarise_cells(filter(cells, occupied | stations >= 5)))
print(wilcox.test(supply ~ class, data = filter(cells, class != "Unsuitable, empty")))

# 05 Plot ----
# styled as Figure S4: white land mask over the 2 km tiles, sampled extent with
# a 3 km margin, 7 pt text
water <- st_read(here("data", "boundaries", "limfjorden", "Limfjorden.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(32632)
xlim <- range(dat$x_utm * 1000) + c(-3000, 3000)
ylim <- range(dat$y_utm * 1000) + c(-3000, 3000)
frame <- st_as_sfc(st_bbox(c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]), crs = st_crs(32632)))
land_mask <- st_difference(frame, st_union(water))

class_cols <- setNames(c(rcartocolor::carto_pal(7, "Teal")[c(6, 3)], "grey80"), levels(cells$class))
# ridge lines a shade darker than the map fills, so the grey class stays visible
line_cols <- setNames(c(class_cols[1:2], "grey55"), names(class_cols))

p_map <- ggplot() +
  geom_tile(data = cells, aes(x, y, fill = class), width = 2000, height = 2000) +
  geom_sf(data = land_mask, fill = "white", colour = NA) +
  geom_sf(data = water, fill = NA, colour = "grey55", linewidth = 0.15) +
  coord_sf(xlim = xlim, ylim = ylim, crs = 32632, expand = FALSE) +
  scale_fill_manual(values = class_cols, name = NULL) +
  guides(fill = guide_legend(theme = theme(legend.key.size = unit(2.5, "mm")))) +
  theme_void(base_size = 7) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    legend.text = element_text(size = 6, colour = "grey30"),
    plot.margin = margin(3, 4, 3, 4)
  )

# one density ridge per class on a log scale, each cell a tick under its ridge
p_supply <- ggplot(cells, aes(supply, fct_rev(class), fill = class, colour = class)) +
  geom_vline(xintercept = supply_ref, linetype = 2, linewidth = 0.3, colour = "grey50") +
  ggridges::geom_density_ridges(
    alpha = 0.6, linewidth = 0.4, scale = 0.9,
    jittered_points = TRUE, position = ggridges::position_points_jitter(width = 0, height = 0),
    point_shape = "|", point_size = 1.5, point_alpha = 0.8
  ) +
  scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_y_discrete(labels = \(x) str_replace(x, ", ", ",\n"), expand = expansion(add = c(0.2, 0.9))) +
  scale_colour_manual(values = line_cols, guide = "none") +
  scale_fill_manual(values = class_cols, guide = "none") +
  labs(x = "Presence-weighted larval supply", y = NULL) +
  theme_light(base_size = 7) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

fig_s8 <- p_map + p_supply +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 7.5, face = "bold", colour = "grey20"))

dir.create(here("output", "figs", "supp"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("output", "figs", "supp", "figS8_accessibility.png"), fig_s8, width = 5.5, height = 3, dpi = 600, bg = "white")
