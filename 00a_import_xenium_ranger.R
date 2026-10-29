#!/usr/bin/env Rscript
# =============================================================================
# 00a_import_xenium_ranger.R -- standard input from Xenium Ranger outputs
#
# Inputs (config.R)
#   XR_DIRS     named vector, sample_id = Xenium Ranger / Xenium Onboard
#               Analysis output folder
#   XR_LABELS   cell labels: one CSV (cell_id, label, sample_id) for all
#               samples, or a named vector of one CSV (cell_id, label) per
#               sample. Cells missing from it are unlabelled (NA).
#   XR_IMPORT_COUNTS  FALSE skips the counts matrix (modules 04 and 05 need it)
#
# What is read from each output folder
#   cells.parquet (if the arrow package is installed) or cells.csv.gz:
#     cell_id, x_centroid, y_centroid (micrometres), transcript_counts.
#   cell_feature_matrix.h5 (hdf5r) or the cell_feature_matrix/ MEX folder:
#     only features of type "Gene Expression" are kept; negative-control
#     probes and codewords, unassigned and deprecated codewords are dropped.
#   n_count = transcript_counts (gene transcripts); if that column is absent,
#   the column sums of the gene-expression matrix.
#
# Writes (INPUT_DIR)
#   <sample>/cells.rds, <sample>/counts.rds, samples.csv (created, or
#   extended with new sample ids), label_map_template.csv (every label found,
#   compartment blank: fill it in and save it as label_map.csv),
#   00_import_report.xlsx (per_sample, labels).
#   An existing label_map.csv is never overwritten.
#
# Scale: Xenium centroids are already in micrometres. The report still prints
#   section width and height (mm) and the median nearest-neighbour distance,
#   which should be of the order of a cell diameter (about 5-20 um).
#
# Usage   Rscript 00a_import_xenium_ranger.R [config.R]
# ASCII only.
# =============================================================================

source("helpers.R")
load_config()
if (length(XR_DIRS) == 0 || is.null(names(XR_DIRS)) || any(names(XR_DIRS) == ""))
  stop("config: XR_DIRS must be a named vector, sample_id = output folder")
require_param("XR_LABELS")
dir.create(INPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
read_cells <- function(d) {
  pq <- file.path(d, "cells.parquet"); cz <- file.path(d, "cells.csv.gz")
  if (file.exists(pq) && requireNamespace("arrow", quietly = TRUE)) {
    x <- as.data.frame(arrow::read_parquet(pq))
  } else if (file.exists(cz)) {
    x <- utils::read.csv(gzfile(cz), stringsAsFactors = FALSE)
  } else if (file.exists(pq)) {
    stop(d, ": only cells.parquet found; install the arrow package to read it")
  } else stop(d, ": neither cells.parquet nor cells.csv.gz found")
  need <- c("cell_id", "x_centroid", "y_centroid")
  if (!all(need %in% names(x))) stop(d, ": cells table lacks ", paste(setdiff(need, names(x)), collapse = ", "))
  x$cell_id <- as.character(x$cell_id)
  x
}

read_h5_matrix <- function(f) {
  h <- hdf5r::H5File$new(f, mode = "r"); on.exit(h$close_all())
  g <- h[["matrix"]]
  shp <- g[["shape"]][]
  m <- Matrix::sparseMatrix(i = g[["indices"]][] + 1L, p = g[["indptr"]][],
                            x = as.numeric(g[["data"]][]), dims = shp)
  ft <- g[["features"]]
  rownames(m) <- ft[["name"]][]
  colnames(m) <- g[["barcodes"]][]
  attr(m, "feature_type") <- ft[["feature_type"]][]
  m
}

read_mex_matrix <- function(d) {
  m  <- Matrix::readMM(gzfile(file.path(d, "matrix.mtx.gz")))
  fe <- utils::read.delim(gzfile(file.path(d, "features.tsv.gz")), header = FALSE,
                          stringsAsFactors = FALSE)
  bc <- readLines(gzfile(file.path(d, "barcodes.tsv.gz")))
  m <- methods::as(m, "CsparseMatrix")
  rownames(m) <- fe[[2]]; colnames(m) <- bc
  attr(m, "feature_type") <- if (ncol(fe) >= 3) fe[[3]] else rep("Gene Expression", nrow(fe))
  m
}

read_counts <- function(d) {
  h5 <- file.path(d, "cell_feature_matrix.h5"); mx <- file.path(d, "cell_feature_matrix")
  m <- if (file.exists(h5)) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) stop("install hdf5r to read ", h5)
    read_h5_matrix(h5)
  } else if (dir.exists(mx)) read_mex_matrix(mx) else
    stop(d, ": no cell_feature_matrix.h5 or cell_feature_matrix/ folder")
  ftype <- attr(m, "feature_type")
  keep <- ftype == "Gene Expression"
  msg("      matrix: %d features, %d gene expression kept", nrow(m), sum(keep))
  m <- m[keep, , drop = FALSE]
  if (anyDuplicated(rownames(m))) rownames(m) <- make.unique(rownames(m))
  m
}

read_labels <- function(sid) {
  if (length(XR_LABELS) > 1 || !is.null(names(XR_LABELS))) {
    f <- XR_LABELS[[sid]]
    if (is.null(f)) stop("XR_LABELS has no file for sample ", sid)
    l <- utils::read.csv(f, stringsAsFactors = FALSE)
  } else {
    l <- utils::read.csv(XR_LABELS, stringsAsFactors = FALSE)
    if (!"sample_id" %in% names(l))
      stop(XR_LABELS, ": a single label file needs a sample_id column")
    l <- l[as.character(l$sample_id) == sid, , drop = FALSE]
  }
  if (!all(c("cell_id", "label") %in% names(l))) stop("label file needs cell_id and label")
  l$cell_id <- as.character(l$cell_id)
  if (anyDuplicated(l$cell_id)) stop(sid, ": a cell_id appears twice in the labels")
  l
}

# -----------------------------------------------------------------------------
report <- list(); lab_counts <- list()
for (sid in names(XR_DIRS)) {
  d <- XR_DIRS[[sid]]
  msg("[%s] %s", sid, d)
  cx <- read_cells(d)
  lb <- read_labels(sid)
  cells <- data.frame(cell_id = cx$cell_id, x_um = as.numeric(cx$x_centroid),
                      y_um = as.numeric(cx$y_centroid),
                      label = lb$label[match(cx$cell_id, lb$cell_id)],
                      n_count = NA_real_, stringsAsFactors = FALSE)
  unmatched <- sum(!lb$cell_id %in% cx$cell_id)
  if (unmatched > 0) msg("      %d labelled cell ids not found in the cells table", unmatched)

  counts <- NULL
  if (XR_IMPORT_COUNTS) {
    m <- read_counts(d)
    j <- match(cells$cell_id, colnames(m))
    if (anyNA(j)) stop(sid, ": ", sum(is.na(j)), " cells have no column in the matrix")
    counts <- m[, j, drop = FALSE]
  }
  if ("transcript_counts" %in% names(cx)) {
    cells$n_count <- as.numeric(cx$transcript_counts)
  } else if (!is.null(counts)) {
    cells$n_count <- as.numeric(Matrix::colSums(counts))
    msg("      no transcript_counts column; n_count = gene-expression column sums")
  } else stop(sid, ": no transcript_counts column and counts not imported")

  od <- file.path(INPUT_DIR, sid); dir.create(od, showWarnings = FALSE)
  saveRDS(cells, file.path(od, "cells.rds"))
  if (!is.null(counts)) saveRDS(counts, file.path(od, "counts.rds"))

  ok <- is.finite(cells$x_um) & is.finite(cells$y_um)
  report[[sid]] <- data.frame(
    sample_id = sid, source = d, n_cells = nrow(cells),
    n_labelled = sum(!is.na(cells$label)), pct_labelled = 100 * mean(!is.na(cells$label)),
    labelled_ids_not_found = unmatched,
    median_n_count = stats::median(cells$n_count, na.rm = TRUE),
    width_mm = diff(range(cells$x_um[ok])) / 1000,
    height_mm = diff(range(cells$y_um[ok])) / 1000,
    median_nn_distance_um = median_nn_distance(cells$x_um[ok], cells$y_um[ok]),
    n_genes = if (is.null(counts)) NA_integer_ else nrow(counts),
    stringsAsFactors = FALSE)
  tl <- table(cells$label)
  lab_counts[[sid]] <- data.frame(sample_id = sid, label = names(tl),
                                  n_cells = as.integer(tl), stringsAsFactors = FALSE)
  msg("      cells=%d labelled=%.1f%% width=%.2f mm height=%.2f mm median NN=%.1f um",
      nrow(cells), report[[sid]]$pct_labelled, report[[sid]]$width_mm,
      report[[sid]]$height_mm, report[[sid]]$median_nn_distance_um)
}

finish_import(report, lab_counts, INPUT_DIR)
write_session_info(INPUT_DIR)
msg("Done.")
