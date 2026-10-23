# Shared helpers for the cockle SDM pipeline.

# compatibility shim: sdmTMBextra 0.0.5 calls fmesher::fm_identical_CRS, renamed
# to fm_crs_is_identical in fmesher >= 0.7; register the old name as an alias.
# required by sdmTMBextra::add_barrier_mesh() in 03_build_mesh.R,
# 04_fit_suitability.R and 07_fit_sdm.R - do not remove unless sdmTMBextra is
# updated
local({
  exps <- .getNamespaceInfo(asNamespace("fmesher"), "exports")
  if (!exists("fm_identical_CRS", envir = exps, inherits = FALSE)) {
    assign("fm_identical_CRS", "fm_crs_is_identical", envir = exps)
  }
})

# sediment.asc codes to Folk 7 seabed substrate classes. The raster is
# EMODnet's 1:250k seabed substrate (GEUS), rasterised in Folk 7 coding:
# checked against the EMODnet WFS layer gtk:seabed_substrate_250k across the
# fjord (86-98% area-wise agreement per code). 1 mud, 6 muddy sand, 5 sand;
# codes 3 (6 km2) and 4 are both mixed sediment. Coarse sediment (1 km2) and
# rock do not occur in the raster. Folk 7 is used rather than Folk 5, which
# merges mud and muddy sand, the main sediment contrast where cockles occur
# (occurrence AIC -15; R/exploration/check_sediment_folk7.R). The most common
# class comes first, as the reference level.
sediment_folk7 <- function(code) {
  classes <- c(
    `1` = "mud", `6` = "muddy_sand", `5` = "sand",
    `3` = "mixed_sediment", `4` = "mixed_sediment"
  )
  factor(unname(classes[as.character(code)]), levels = c("mud", "muddy_sand", "sand", "mixed_sediment"))
}

# swept-area ratio (SAR) of the cockle fishery (HFisk, data/sar): swept area
# per 25 x 25 m cell and fishing season, 2016-17 to 2024-25; only fished
# cells are listed (absent = unfished). Summed onto a 250 m grid (UTM 32N):
# SAR = swept area / 62,500 m2. Returns a terra SpatRaster, one layer per
# season named "16-17", ..., "24-25". Used by 02_clean_cockles.R and
# fig2_predictors.R
sar_seasons <- function() {
  sar_dir <- here::here("data", "sar", "HFisk_SAR_CPUE_2016-2025")
  if (!dir.exists(sar_dir)) unzip(here::here("data", "sar", "HFisk_SAR_CPUE_2016-2025.zip"), exdir = here::here("data", "sar"))
  template <- terra::rast(xmin = 440000, xmax = 590000, ymin = 6250000, ymax = 6335000, resolution = 250, crs = "EPSG:32632")
  files <- list.files(sar_dir, pattern = "\\.shp$", full.names = TRUE)
  stack <- terra::rast(lapply(files, \(f) {
    cells <- sf::st_read(f, quiet = TRUE) |>
      sf::st_transform(32632) |>
      sf::st_centroid() |>
      suppressWarnings() # attributes are constant over each 25 m cell
    terra::rasterize(terra::vect(cells), template, field = "AreaSum", fun = "sum", background = 0) / 250^2
  }))
  names(stack) <- regmatches(basename(files), regexpr("\\d{2}-\\d{2}", basename(files)))
  stack
}
