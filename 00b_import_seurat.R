#!/usr/bin/env Rscript
# =============================================================================
# 00b_import_seurat.R -- standard input from a Seurat v5 object
#
# Inputs (config.R)
#   SEURAT_RDS         saved Seurat object (saveRDS), one or many samples
#   SEURAT_LABEL_COL   meta.data column with the cell labels
#   SEURAT_SAMPLE_COL  meta.data column with the sample id
#   SEURAT_COUNT_COL   meta.data column with total counts per cell; NULL =
#                      nCount_<assay>
#   SEURAT_ASSAY       NULL = the default assay
#   SEURAT_IMPORT_COUNTS  FALSE skips the counts matrix (04 and 05 need it)
#   COORD_SCALE        REQUIRED. Micrometres per coordinate unit. Objects
#                      built with LoadXenium() store centroids in
#                      micrometres: COORD_SCALE <- 1. See "Checking the
#                      scale" below.
#
# What is read
#   Centroids from every field of view (Images(obj)) with
#   GetTissueCoordinates(); in a merged object each sample has its own field
#   of view and cell names carry the prefixes added by merge(), so cells are
#   matched to meta.data by name. Every cell must appear in exactly one field
#   of view. Coordinates are multiplied by COORD_SCALE. Count layers split by
#   merge() are joined (JoinLayers) before the counts are extracted. Only
#   SeuratObject is needed, not Seurat.
#
# Checking the scale
#   After scaling, the script prints per sample the section width and height
#   in mm and the median nearest-neighbour distance in um (also in
#   00_import_report.xlsx). Width and height should match the physical
#   section; the nearest-neighbour distance should be of the order of a cell
#   diameter (roughly 5-20 um in solid tissue). A value in the hundreds or
#   thousands means the coordinates are pixels: set COORD_SCALE to the image
#   pixel size in um per pixel (for Xenium morphology images it is recorded
#   as pixel_size in experiment.xenium).
#
# Writes (INPUT_DIR)
#   <sample>/cells.rds, <sample>/counts.rds, samples.csv,
#   label_map_template.csv (if no label_map.csv exists), 00_import_report.xlsx
#
# Usage   Rscript 00b_import_seurat.R [config.R]
# ASCII only.
# =============================================================================

source("helpers.R")
load_config()
require_param("SEURAT_RDS", "SEURAT_LABEL_COL", "SEURAT_SAMPLE_COL")
if (is.null(COORD_SCALE) || !is.numeric(COORD_SCALE) || length(COORD_SCALE) != 1 ||
    !is.finite(COORD_SCALE) || COORD_SCALE <= 0)
  stop("config: COORD_SCALE (micrometres per coordinate unit) must be set; ",
       "1 for objects built with LoadXenium(). See the header of this script.")
suppressPackageStartupMessages(library(SeuratObject))
dir.create(INPUT_DIR, showWarnings = FALSE, recursive = TRUE)

obj   <- readRDS(SEURAT_RDS)
assay <- if (is.null(SEURAT_ASSAY)) DefaultAssay(obj) else SEURAT_ASSAY
meta  <- obj[[]]
for (cc in c(SEURAT_LABEL_COL, SEURAT_SAMPLE_COL))
  if (!cc %in% names(meta)) stop("meta.data has no column ", cc)
count_col <- if (is.null(SEURAT_COUNT_COL)) paste0("nCount_", assay) else SEURAT_COUNT_COL
if (!count_col %in% names(meta)) stop("meta.data has no column ", count_col)
msg("[1] %d cells, assay %s, %d field(s) of view", nrow(meta), assay, length(Images(obj)))

# -----------------------------------------------------------------------------
# Centroids from every field of view
# -----------------------------------------------------------------------------
if (length(Images(obj)) == 0) stop("the object has no field of view (Images() is empty)")
xy <- do.call(rbind, lapply(Images(obj), function(im) {
  tc <- GetTissueCoordinates(obj[[im]])
  if (!"cell" %in% names(tc)) tc$cell <- rownames(tc)
  data.frame(cell = as.character(tc$cell), x = as.numeric(tc$x), y = as.numeric(tc$y),
             fov = im, stringsAsFactors = FALSE)
}))
dup <- unique(xy$cell[duplicated(xy$cell)])
if (length(dup)) stop(length(dup), " cells appear in more than one field of view, e.g. ",
                      paste(head(dup, 3), collapse = ", "))
j <- match(rownames(meta), xy$cell)
if (anyNA(j)) stop(sum(is.na(j)), " cells have no centroid in any field of view, e.g. ",
                   paste(head(rownames(meta)[is.na(j)], 3), collapse = ", "))
xy <- xy[j, ]
xy$x <- xy$x * COORD_SCALE; xy$y <- xy$y * COORD_SCALE

# -----------------------------------------------------------------------------
# Counts
# -----------------------------------------------------------------------------
counts <- NULL
if (SEURAT_IMPORT_COUNTS) {
  a <- obj[[assay]]
  if (inherits(a, "Assay5")) {
    if (length(grep("^counts", Layers(a))) > 1) a <- JoinLayers(a)
    counts <- LayerData(a, layer = "counts")
  } else {
    counts <- GetAssayData(a, slot = "counts")
  }
  counts <- methods::as(counts, "CsparseMatrix")
  counts <- counts[, rownames(meta), drop = FALSE]
}

# -----------------------------------------------------------------------------
# One folder per sample
# -----------------------------------------------------------------------------
sids <- as.character(meta[[SEURAT_SAMPLE_COL]])
report <- list(); lab_counts <- list()
for (sid in unique(sids)) {
  sel <- which(sids == sid)
  cells <- data.frame(cell_id = rownames(meta)[sel], x_um = xy$x[sel], y_um = xy$y[sel],
                      label = as.character(meta[[SEURAT_LABEL_COL]][sel]),
                      n_count = as.numeric(meta[[count_col]][sel]),
                      stringsAsFactors = FALSE)
  od <- file.path(INPUT_DIR, sid); dir.create(od, showWarnings = FALSE)
  saveRDS(cells, file.path(od, "cells.rds"))
  if (!is.null(counts)) saveRDS(counts[, sel, drop = FALSE], file.path(od, "counts.rds"))
  ok <- is.finite(cells$x_um) & is.finite(cells$y_um)
  report[[sid]] <- data.frame(
    sample_id = sid, source = SEURAT_RDS, fields_of_view = paste(unique(xy$fov[sel]), collapse = ";"),
    n_cells = nrow(cells), n_labelled = sum(!is.na(cells$label)),
    pct_labelled = 100 * mean(!is.na(cells$label)), coord_scale = COORD_SCALE,
    median_n_count = stats::median(cells$n_count, na.rm = TRUE),
    width_mm = diff(range(cells$x_um[ok])) / 1000,
    height_mm = diff(range(cells$y_um[ok])) / 1000,
    median_nn_distance_um = median_nn_distance(cells$x_um[ok], cells$y_um[ok]),
    n_genes = if (is.null(counts)) NA_integer_ else nrow(counts),
    stringsAsFactors = FALSE)
  tl <- table(cells$label)
  lab_counts[[sid]] <- data.frame(sample_id = sid, label = names(tl),
                                  n_cells = as.integer(tl), stringsAsFactors = FALSE)
  msg("[%s] cells=%d labelled=%.1f%% width=%.2f mm height=%.2f mm median NN=%.1f um",
      sid, nrow(cells), report[[sid]]$pct_labelled, report[[sid]]$width_mm,
      report[[sid]]$height_mm, report[[sid]]$median_nn_distance_um)
  nn <- report[[sid]]$median_nn_distance_um
  if (is.finite(nn) && (nn < 2 || nn > 50))
    msg("[!] %s: median nearest-neighbour distance %.1f um is implausible for cells; check COORD_SCALE", sid, nn)
}
finish_import(report, lab_counts, INPUT_DIR)
write_session_info(INPUT_DIR)
msg("Done.")
