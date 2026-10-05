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
