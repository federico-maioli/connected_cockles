# Compute the FLOW matrix (export probability, source -> sink) per year and pooled.
#
#   flow[i, j] = agents released at i that settle at j  /  total agents released at i
#              = cmn[i, j] / release[i]
#
# This is the "downstream / export probability relative to total release" matrix
# (conmatprobrel in the collaborator's script, the default relsel = 1 case).
# Rows = source node, columns = sink node, on the 65 x 36 = 2340 grid.
#
# Writes data/intermediate/flow_matrix_<year>.rds for 2010-2016 and
# data/intermediate/flow_matrix_all.rds for the pooled run.

library(tidyverse)
library(here)

# read a raw connectivity matrix (agent counts, no header)
read_cmn <- function(label) {
  path <- here("data", "raw", "connectivity", "cockles_matrices", paste0("ABM_cmn_", label, ".csv"))
  cmn <- as.matrix(read_csv(path, col_names = FALSE, show_col_types = FALSE))
  dimnames(cmn) <- NULL
  cmn[is.na(cmn)] <- 0
  cmn
}

# read total release per source node (single column, no header)
read_release <- function(label) {
  path <- here("data", "raw", "connectivity", "cockles_release", paste0("ABM_rel_", label, ".csv"))
  release <- read_csv(path, col_names = FALSE, show_col_types = FALSE)[[1]]
  as.numeric(release)
}

# divide each row i by release[i]; rows with zero release stay 0
compute_flow <- function(cmn, release) {
  denom <- matrix(release, nrow = nrow(cmn), ncol = ncol(cmn), byrow = FALSE)
  flow <- ifelse(denom > 0, cmn / denom, 0)
  flow[is.nan(flow)] <- 0
  flow <- round(flow, 7)
  dimnames(flow) <- NULL
  flow
}

# 01 Set up ----
labels <- c(as.character(2010:2016), "all")
out_dir <- here("data", "intermediate")

# 02 Compute and save one flow matrix per label ----
summaries <- map(labels, function(label) {
  cmn <- read_cmn(label)
  release <- read_release(label)

  stopifnot(nrow(cmn) == ncol(cmn)) # square
  stopifnot(length(release) == nrow(cmn)) # aligns with sources (rows)

  flow <- compute_flow(cmn, release)

  out_path <- file.path(out_dir, paste0("flow_matrix_", label, ".rds"))
  saveRDS(flow, out_path)

  row_export <- rowSums(flow) # total export prob per source
  tibble(
    label = label,
    n_nodes = nrow(flow),
    n_released = sum(release > 0),
    max_flow = max(flow),
    max_export = max(row_export),
    mean_self_retention = mean(diag(flow)[release > 0])
  )
})

# 03 Report ----
summaries |>
  list_rbind() |>
  mutate(across(c(max_flow, max_export, mean_self_retention), \(x) round(x, 5))) |>
  print(n = Inf)

