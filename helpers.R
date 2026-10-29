# =============================================================================
# helpers.R -- functions shared by every script in xenium-spatial-niche
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(igraph)
  library(dbscan)
  library(ggplot2)
  library(openxlsx)
})

# =============================================================================
# Configuration and messages
# =============================================================================

# Sources the configuration file into the global environment. The file is the
# first command-line argument if one is given, otherwise config.R in the
# working directory (the repository root).
load_config <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  cfg  <- if (length(args) >= 1) args[1] else "config.R"
  if (!file.exists(cfg)) stop("configuration file not found: ", cfg)
  sys.source(cfg, envir = globalenv())
  msg("configuration: %s", normalizePath(cfg))
  invisible(cfg)
}

msg <- function(fmt, ...) message(sprintf(fmt, ...))

# Stops with a clear message if a required configuration value is missing.
require_param <- function(...) {
  nm <- c(...)
  for (n in nm) {
    if (!exists(n, envir = globalenv()) || is.null(get(n, envir = globalenv())))
      stop("config: ", n, " must be set for this script")
  }
  invisible(TRUE)
}

module_dir <- function(results_dir, name) {
  d <- file.path(results_dir, name)
  dir.create(file.path(d, "cache"), showWarnings = FALSE, recursive = TRUE)
  d
}

# =============================================================================
# Standard input
# =============================================================================
# input/samples.csv        sample_id [, any sample-level covariates]
# input/label_map.csv      label, compartment
# input/<sample>/cells.rds data.frame: cell_id, x_um, y_um, label, n_count
# input/<sample>/counts.rds dgCMatrix genes x cells, columns in cells.rds order

read_samples <- function(input_dir, keep = NULL) {
  f <- file.path(input_dir, "samples.csv")
  if (!file.exists(f)) stop("missing ", f, "; run an import script first")
  s <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample_id" %in% names(s)) stop(f, " needs a sample_id column")
  s$sample_id <- as.character(s$sample_id)
  if (!is.null(keep)) {
    miss <- setdiff(keep, s$sample_id)
    if (length(miss)) stop("SAMPLES not in samples.csv: ", paste(miss, collapse = ", "))
    s <- s[match(keep, s$sample_id), , drop = FALSE]
  }
  s
}

read_label_map <- function(input_dir) {
  f <- file.path(input_dir, "label_map.csv")
  if (!file.exists(f)) stop("missing ", f, " (columns: label, compartment)")
  m <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("label", "compartment") %in% names(m)))
    stop(f, " needs columns label and compartment")
  m$label <- as.character(m$label)
  m$compartment <- as.character(m$compartment)
  if (anyDuplicated(m$label)) stop(f, ": a label appears more than once")
  if (any(is.na(m$compartment) | m$compartment == ""))
    stop(f, ": every label needs a compartment")
  m
}

# Reads one sample and attaches its compartment. Labels absent from the label
# map are treated as unlabelled (NA) and reported once.
load_sample <- function(sample_id, input_dir, label_map, counts = FALSE) {
  d <- file.path(input_dir, sample_id)
  cells <- readRDS(file.path(d, "cells.rds"))
  need <- c("cell_id", "x_um", "y_um", "label", "n_count")
  miss <- setdiff(need, names(cells))
  if (length(miss)) stop(sample_id, ": cells.rds lacks ", paste(miss, collapse = ", "))
  cells$cell_id <- as.character(cells$cell_id)
  cells$label   <- as.character(cells$label)
  unk <- setdiff(unique(cells$label[!is.na(cells$label)]), label_map$label)
  if (length(unk)) {
    msg("    %s: %d label(s) not in label_map.csv, treated as unlabelled: %s",
        sample_id, length(unk), paste(head(unk, 10), collapse = ", "))
    cells$label[cells$label %in% unk] <- NA_character_
  }
  cells$compartment <- label_map$compartment[match(cells$label, label_map$label)]
  out <- list(cells = cells, counts = NULL)
  if (counts) {
    f <- file.path(d, "counts.rds")
    if (!file.exists(f)) stop(sample_id, ": counts.rds is needed by this script")
    m <- readRDS(f)
    m <- methods::as(m, "CsparseMatrix")
    if (ncol(m) != nrow(cells)) stop(sample_id, ": counts.rds columns do not match cells.rds")
    if (!is.null(colnames(m)) && !identical(colnames(m), cells$cell_id))
      stop(sample_id, ": counts.rds column names differ from cells.rds cell_id")
    out$counts <- m
  }
  out
}

# Shared last step of both import scripts: samples.csv (created or extended),
# label_map_template.csv when no label_map.csv exists yet, and the import
# report. An existing label_map.csv is never overwritten.
finish_import <- function(report, lab_counts, input_dir) {
  RP <- do.call(rbind, report); LC <- do.call(rbind, lab_counts)
  sf <- file.path(input_dir, "samples.csv")
  if (file.exists(sf)) {
    s <- utils::read.csv(sf, stringsAsFactors = FALSE, check.names = FALSE)
    new <- setdiff(RP$sample_id, s$sample_id)
    if (length(new)) {
      add <- s[0, , drop = FALSE][seq_along(new), , drop = FALSE]
      add$sample_id <- new
      s <- rbind(s, add)
    }
  } else s <- data.frame(sample_id = RP$sample_id)
  utils::write.csv(s, sf, row.names = FALSE)

  labs <- sort(unique(LC$label))
  lm_f <- file.path(input_dir, "label_map.csv")
  if (file.exists(lm_f)) {
    lm <- utils::read.csv(lm_f, stringsAsFactors = FALSE)
    miss <- setdiff(labs, lm$label)
    if (length(miss)) msg("[!] labels not in label_map.csv (will be treated as unlabelled): %s",
                          paste(miss, collapse = ", "))
  } else {
    utils::write.csv(data.frame(label = labs, compartment = ""),
                     file.path(input_dir, "label_map_template.csv"), row.names = FALSE)
    msg("[!] no label_map.csv yet: fill in the compartment column of %s and save it as %s",
        file.path(input_dir, "label_map_template.csv"), lm_f)
  }
  write_xlsx(list(per_sample = RP, labels = LC), file.path(input_dir, "00_import_report.xlsx"))
}

# Tissue component per cell, as written by 01_geometry_qc.R (0 = debris).
read_components <- function(results_dir, sample_id, cell_id) {
  f <- file.path(results_dir, "01_geometry_qc", "components.rds")
  if (!file.exists(f)) stop("missing ", f, "; run 01_geometry_qc.R first")
  cmp <- readRDS(f)
  cmp <- cmp[cmp$sample_id == sample_id, , drop = FALSE]
  if (nrow(cmp) == 0) stop(sample_id, ": not in components.rds; rerun 01_geometry_qc.R")
  out <- cmp$component[match(cell_id, cmp$cell_id)]
  out[is.na(out)] <- 0L
  as.integer(out)
}

# =============================================================================
# Geometry
# =============================================================================

# Connected tissue components on a binned occupancy grid. A bin is tissue if it
# holds at least min_cells_bin cells; occupied bins are joined by rook (4) or
# queen (8) adjacency; components smaller than min_region_bins bins are debris.
# Returns an integer per cell: component id, or 0 for debris / sparse bins.
tissue_components <- function(x, y, bin_um, min_cells_bin, connectivity,
                              min_region_bins) {
  if (!connectivity %in% c(4, 8)) stop("connectivity must be 4 or 8")
  bx <- floor(x / bin_um); by <- floor(y / bin_um)
  bkey <- paste(bx, by, sep = "_")
  tb   <- table(bkey)
  occ  <- names(tb)[tb >= min_cells_bin]
  if (length(occ) == 0) return(rep(0L, length(x)))

  obx <- as.integer(sub("_.*$", "", occ))
  oby <- as.integer(sub("^.*_", "", occ))
  n   <- length(occ); idx <- seq_len(n)

  if (connectivity == 4) { dx <- c(1, 0);       dy <- c(0, 1) }
  else                   { dx <- c(1, 0, 1, 1); dy <- c(0, 1, 1, -1) }

  ei <- integer(0); ej <- integer(0)
  for (k in seq_along(dx)) {
    hit <- match(paste(obx + dx[k], oby + dy[k], sep = "_"), occ)
    ok  <- !is.na(hit)
    ei  <- c(ei, idx[ok]); ej <- c(ej, hit[ok])
  }
  g <- igraph::make_empty_graph(n = n, directed = FALSE)
  if (length(ei) > 0) g <- igraph::add_edges(g, rbind(ei, ej))
  memb <- as.integer(igraph::components(g)$membership)

  sz <- table(memb)
  memb[memb %in% as.integer(names(sz)[sz < min_region_bins])] <- 0L
  # Component ids are igraph's, not renumbered after debris removal, so that
  # permutation strata keep the same order as in the original analysis.

  out <- memb[match(bkey, occ)]
  out[is.na(out)] <- 0L
  as.integer(out)
}

# Undirected edges between points within r of each other, as a two-column
# integer matrix with i < j. Exact (kd-tree fixed-radius search).
radius_edges <- function(x, y, r) {
  n <- length(x)
  if (n < 2) return(matrix(integer(0), ncol = 2))
  nn  <- dbscan::frNN(cbind(x, y), eps = r, sort = FALSE)
  len <- lengths(nn$id)
  if (sum(len) == 0) return(matrix(integer(0), ncol = 2))
  i <- rep.int(seq_len(n), len)
  j <- unlist(nn$id, use.names = FALSE)
  keep <- i < j
  cbind(i[keep], j[keep])
}

# For each query point, the target points within r: a data.frame with one row
# per (query, target) pair. Query and target sets are separate.
bipartite_pairs <- function(qx, qy, tx, ty, r) {
  if (length(qx) == 0 || length(tx) == 0)
    return(data.frame(q = integer(0), t = integer(0)))
  nn  <- dbscan::frNN(cbind(tx, ty), eps = r, query = cbind(qx, qy), sort = FALSE)
  len <- lengths(nn$id)
  data.frame(q = rep.int(seq_along(qx), len),
             t = unlist(nn$id, use.names = FALSE))
}

# Number of target points within r of each query point.
count_within <- function(qx, qy, tx, ty, r) {
  if (length(qx) == 0) return(integer(0))
  if (length(tx) == 0) return(integer(length(qx)))
  nn <- dbscan::frNN(cbind(tx, ty), eps = r, query = cbind(qx, qy), sort = FALSE)
  as.integer(lengths(nn$id))
}

# Distance from each query point to its nearest target point; NA beyond max_d.
nearest_dist <- function(qx, qy, tx, ty, max_d) {
  out <- rep(NA_real_, length(qx))
  if (length(qx) == 0 || length(tx) == 0) return(out)
  k <- dbscan::kNN(cbind(tx, ty), k = 1, query = cbind(qx, qy))
  d <- as.numeric(k$dist[, 1])
  d[!is.finite(d) | d > max_d] <- NA_real_
  d
}

# Median distance from each point to its nearest other point.
median_nn_distance <- function(x, y) {
  if (length(x) < 2) return(NA_real_)
  stats::median(dbscan::kNN(cbind(x, y), k = 1)$dist[, 1])
}

# =============================================================================
# Permutation
# =============================================================================

# Splits indices into shuffle blocks once, so each permutation is a cheap
# relabel. strata is any vector; labels move only within equal strata.
strata_index <- function(strata) split(seq_along(strata), strata)

shuffle_within <- function(labels, idx_by_stratum) {
  out <- labels
  for (s in idx_by_stratum) if (length(s) > 1L) out[s] <- labels[sample(s)]
  out
}

# Counts edges joining A state a to B state b. codes: 0 = not tested,
# 1..nA = A states, nA+1..nA+nB = B states. Each edge is counted once.
cross_counts <- function(ci, cj, nA, nB) {
  a1 <- ci >= 1L & ci <= nA; b2 <- cj > nA
  a2 <- cj >= 1L & cj <= nA; b1 <- ci > nA
  u <- a1 & b2; v <- a2 & b1
  acode <- c(ci[u], cj[v])
  bcode <- c(cj[u], ci[v]) - nA
  if (length(acode) == 0) return(matrix(0L, nA, nB))
  matrix(tabulate((acode - 1L) * nB + bcode, nbins = nA * nB),
         nrow = nA, ncol = nB, byrow = TRUE)
}

# Observed A x B edge counts and their permutation null. Returns matrices
# (A states x B states): obs, exp, sd, log2 enrichment, z and two-sided p.
edge_permutation <- function(code, idx_by_stratum, ei, ej, nA, nB, nperm) {
  obs   <- cross_counts(code[ei], code[ej], nA, nB)
  sum_x <- matrix(0, nA, nB); sum_x2 <- matrix(0, nA, nB)
  ge    <- matrix(0L, nA, nB); le <- matrix(0L, nA, nB)
  for (b in seq_len(nperm)) {
    perm <- shuffle_within(code, idx_by_stratum)
    cnt  <- cross_counts(perm[ei], perm[ej], nA, nB)
    sum_x <- sum_x + cnt; sum_x2 <- sum_x2 + cnt^2
    ge <- ge + (cnt >= obs); le <- le + (cnt <= obs)
  }
  mu  <- sum_x / nperm
  sdv <- sqrt(pmax(sum_x2 / nperm - mu^2, 0))
  list(obs = obs, exp = mu, sd = sdv,
       log2fc = log2((obs + 1) / (mu + 1)),
       z = ifelse(sdv > 0, (obs - mu) / sdv, NA_real_),
       p = pmin(2 * pmin((1 + ge) / (nperm + 1), (1 + le) / (nperm + 1)), 1))
}

# =============================================================================
# Statistics
# =============================================================================

# Wilcoxon signed-rank test against zero across samples, with the
# Hodges-Lehmann location estimate and its 95% CI. NULL below min_n samples.
hl_test <- function(v, min_n) {
  v <- v[is.finite(v)]
  if (length(v) < min_n) return(NULL)
  wt <- suppressWarnings(stats::wilcox.test(v, mu = 0, exact = FALSE,
                                            conf.int = TRUE))
  data.frame(n_samples = length(v), estimate = as.numeric(wt$estimate),
             ci_low = wt$conf.int[1], ci_high = wt$conf.int[2],
             median = stats::median(v), frac_positive = mean(v > 0),
             p_raw = wt$p.value, stringsAsFactors = FALSE)
}

stars <- function(q) ifelse(is.na(q), "", ifelse(q < 0.001, "***",
                     ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", ""))))

# Stars only where q < 0.05 AND the estimate lies in the chosen direction.
directional_stars <- function(estimate, q, direction) {
  if (!direction %in% c("positive", "negative", "both"))
    stop("STAR_DIRECTION must be positive, negative or both")
  ok <- switch(direction,
               positive = estimate > 0,
               negative = estimate < 0,
               both     = rep(TRUE, length(estimate)))
  ifelse(ok %in% TRUE, stars(q), "")
}

# BH correction within families: rows whose pair is pre-specified form one
# family, every other row the exploratory family. Extra grouping columns (for
# example the null model or the metric) split both families further.
bh_by_family <- function(p, is_primary, groups = NULL) {
  key <- paste(if (is.null(groups)) "" else groups, is_primary, sep = "|")
  q <- rep(NA_real_, length(p))
  for (k in unique(key)) {
    s <- key == k & !is.na(p)
    if (any(s)) q[s] <- stats::p.adjust(p[s], method = "BH")
  }
  q
}

is_primary_pair <- function(a, b, pairs) {
  if (length(pairs) == 0) return(rep(FALSE, length(a)))
  key <- vapply(pairs, function(z) paste(z[1], z[2], sep = "|"), character(1))
  paste(a, b, sep = "|") %in% key
}

# Interpolated distance at which a curve first crosses zero from above. NA
# unless the curve starts positive.
crossing <- function(x, y) {
  o <- order(x); x <- x[o]; y <- y[o]
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 2 || y[1] <= 0) return(NA_real_)
  for (i in seq_len(length(x) - 1)) {
    if (y[i] > 0 && y[i + 1] <= 0)
      return(x[i] + (x[i + 1] - x[i]) * y[i] / (y[i] - y[i + 1]))
  }
  NA_real_
}

# Within-sample depth strata: n_strata quantile bins of n_count.
depth_strata <- function(n_count, n_strata) {
  br <- unique(stats::quantile(n_count, probs = seq(0, 1, length.out = n_strata + 1),
                               na.rm = TRUE))
  if (length(br) < 3) return(rep(1L, length(n_count)))
  as.integer(cut(n_count, breaks = br, include.lowest = TRUE))
}

# Size-weighted mean over strata of (mean in arm) - (mean outside arm). Strata
# with fewer than min_per_cell cells in either arm are skipped.
strat_diff <- function(value, arm, strata, min_per_cell = 3) {
  num <- 0; den <- 0
  for (s in unique(strata)) {
    a <- arm & strata == s; b <- (!arm) & strata == s
    if (sum(a) < min_per_cell || sum(b) < min_per_cell) next
    w <- sum(a) + sum(b)
    num <- num + w * (mean(value[a]) - mean(value[b])); den <- den + w
  }
  if (den > 0) num / den else NA_real_
}

# Partial Pearson correlation of every column of Y (cells x genes, sparse,
# binary detection) with xv, given covariates Z (cells x k, with intercept).
# All genes at once by residualisation (Frisch-Waugh).
partial_cor <- function(Y, xv, Z) {
  ZtZi <- solve(crossprod(Z))
  rx   <- as.numeric(xv - Z %*% (ZtZi %*% crossprod(Z, xv)))
  ZtY  <- as.matrix(crossprod(Z, Y))
  num  <- as.numeric(crossprod(Y, rx))
  ssY  <- Matrix::colSums(Y)
  ssR  <- ssY - colSums(ZtY * (ZtZi %*% ZtY))
  out  <- num / (sqrt(sum(rx^2)) * sqrt(pmax(ssR, 0)))
  out[!is.finite(out) | ssR <= 1e-10] <- NA_real_
  out
}

# =============================================================================
# Checkpoints
# =============================================================================
# A checkpoint stores the build parameters with the result. It is reused only
# when those parameters are unchanged, so a stale cache is never read silently.

read_checkpoint <- function(path, params, force_rebuild, replot_only) {
  if (force_rebuild && !replot_only) return(NULL)
  if (!file.exists(path)) {
    if (replot_only) stop("REPLOT_ONLY is TRUE but no checkpoint exists: ", path)
    return(NULL)
  }
  ck <- readRDS(path)
  if (!identical(ck$params, params)) {
    if (replot_only) stop("REPLOT_ONLY is TRUE but the checkpoint was built with ",
                          "different parameters: ", path)
    msg("[cache] parameters changed since the checkpoint; rebuilding")
    return(NULL)
  }
  msg("[cache] reading %s (set FORCE_REBUILD <- TRUE to recompute)", path)
  ck$data
}

write_checkpoint <- function(path, params, data) {
  saveRDS(list(params = params, data = data), path)
  invisible(path)
}

# =============================================================================
# Output
# =============================================================================

# Unformatted workbook, one sheet per list element.
write_xlsx <- function(sheets, path) {
  wb <- openxlsx::createWorkbook()
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    d <- sheets[[nm]]
    if (is.null(d) || (is.data.frame(d) && nrow(d) == 0))
      d <- data.frame(note = "no rows")
    openxlsx::writeData(wb, nm, d)
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  msg("    written: %s", path)
  invisible(path)
}

theme_niche <- function(base_size) {
  theme_classic(base_size = base_size) +
    theme(axis.line  = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(linewidth = 0.4, colour = "black"),
          axis.text  = element_text(colour = "black"),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = base_size))
}

# Opens a PDF device
open_pdf <- function(path, width, height) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
}

close_pdf <- function(path) {
  grDevices::dev.off()
  msg("    written: %s", path)
}

write_session_info <- function(out_dir) {
  f <- file.path(out_dir, "session_info.txt")
  writeLines(utils::capture.output(utils::sessionInfo()), f)
  invisible(f)
}

HEAT_LOW  <- "#1965B0"
HEAT_MID  <- "white"
HEAT_HIGH <- "#DC050C"
