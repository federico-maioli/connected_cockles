# Compute the FLOW matrix (export probability, source -> sink) from ABM_cmn_all.
#
#   flow[i, j] = agents released at i that settle at j  /  total agents released at i
#              = cmn[i, j] / release[i]
#
# This is the "downstream / export probability relative to total release" matrix
# (conmatprobrel in the collaborator's script, the default relsel = 1 case).
# Rows = source node, columns = sink node, on the 65 x 36 = 2340 grid.

library(tidyverse)
library(here)

# 01 Load inputs ----
# raw connectivity matrix (agent counts), no header
cmn_path <- here("data", "raw", "cockles_matrices", "ABM_cmn_all.csv")
cmn <- as.matrix(read_csv(cmn_path, col_names = FALSE, show_col_types = FALSE))
dimnames(cmn) <- NULL
cmn[is.na(cmn)] <- 0

# total release per source node (single column, no header)
rel_path <- here("data", "extra", "cockles_release", "ABM_rel_all.csv")
release <- read_csv(rel_path, col_names = FALSE, show_col_types = FALSE)[[1]]
release <- as.numeric(release)

stopifnot(nrow(cmn) == ncol(cmn)) # square
stopifnot(length(release) == nrow(cmn)) # aligns with sources (rows)

# 02 Compute flow matrix ----
# divide each row i by release[i]; leave rows with zero release as 0
denom <- matrix(release, nrow = nrow(cmn), ncol = ncol(cmn), byrow = FALSE)
flow <- ifelse(denom > 0, cmn / denom, 0)
flow <- round(flow, 7)
flow[is.nan(flow)] <- 0
dimnames(flow) <- NULL

# 03 Checks ----
row_export <- rowSums(flow) # total export prob per source (should be in [0, 1])
cat("Flow dimensions:", nrow(flow), "x", ncol(flow), "\n")
cat("Value range:", range(flow), "\n")
cat("Nodes with non-zero release:", sum(release > 0), "/", length(release), "\n")
cat("Row export prob range:", round(range(row_export), 4), "\n")
cat("Self-retention (mean of diagonal, non-zero release):",
  round(mean(diag(flow)[release > 0]), 6), "\n")

# 04 Save ----
out_path <- here("data", "intermediate", "flow_matrix.rds")
saveRDS(flow, out_path)
cat("Saved flow matrix to:", out_path, "\n")
