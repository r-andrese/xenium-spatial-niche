#!/usr/bin/env Rscript
# =============================================================================
# 05_neighbourhood_genes.R -- genes tied to the local neighbourhood
#
# Question
#   Which genes are detected more often in index cells (INDEX_STATES) that
#   have more neighbour cells (NEIGHBOUR_STATES) within RADIUS_UM? This is an
#   unsupervised screen for neighbourhood-associated genes. The result is an
#   association with the neighbourhood, not evidence of signalling. The
#   reverse direction is a second run with the two sets swapped.
#
# Per-sample statistic (the sample is the replicate unit)
#   For every index cell: detection (count > 0) of each gene; exposure =
#   log1p(number of neighbour-set cells within RADIUS_UM); covariates =
#   log(n_count) of the index cell (depth) and log1p(number of all in-tissue
#   cells within RADIUS_UM, excluding itself; local density). The statistic 
#   is the partial Pearson correlation of detection with exposure given the 
#   covariates, for all genes at once by residualisation. Detection is used 
#   rather than counts because per-cell depth is low relative to panel size. 
#   A gene contributes in a sample only if detected in >= MIN_DETECT_FRAC of 
#   its index cells.
#
# Cohort summary and selection
#   Per gene: median per-sample correlation, fraction of samples positive,
#   Wilcoxon signed-rank against 0 (if tested in >= MIN_SAMPLES_GENE
#   samples), BH across genes. Selected = tested in >= MIN_SAMPLES_GENE
#   samples, BH q < FDR_MAX, median correlation > 0, positive in >=
#   MIN_FRAC_POS of samples, and passes the spillover filter; ranked by
#   median correlation, top N_TOP kept.
#
# Spillover filter
#   Imperfect segmentation passes neighbour transcripts into index cells, so
#   neighbour-compartment genes would otherwise be selected.
#   WITH A REFERENCE (SPILL_REFERENCE: gene, then one column of mean
#     expression per compartment; a matched scRNA-seq reference is
#     preferred): index-compartment mean >= SPILL_RATIO_MIN x the highest
#     other compartment. Genes absent from the reference fail.
#   WITHOUT (in situ): median over samples of the detection rate in index
#     cells with no neighbour-set cell within RADIUS_UM >= SPILL_RATIO_MIN x
#     the median detection rate in neighbour-compartment cells. 
#   Genes removed by the filter are reported.
#
# Inclusion
#   Samples: median n_count of labelled cells >= MIN_MEDIAN_DEPTH. Debris
#   excluded (01). Per sample: >= MIN_INDEX_CELLS index cells (labelled, in
#   tissue, n_count > 0) and >= MIN_NEIGHBOUR_CELLS neighbour-set cells.
#
# Reads   INPUT_DIR (cells.rds, counts.rds), RESULTS_DIR/01_geometry_qc/,
#         SPILL_REFERENCE (optional)
# Writes  RESULTS_DIR/05_neighbourhood_genes/
#           05_neighbourhood_genes.xlsx  sample_inclusion, per_sample_cor,
#                                        gene_summary, selected_genes,
#                                        spillover_removed, selection_funnel
#           05_gene_volcano.pdf          median correlation vs -log10 q
#           05_selected_heatmap.pdf      per-sample correlations, selected genes
#           cache/05_checkpoint.rds      per-sample correlations: thresholds
#                                        can be changed without re-reading data
#
# Usage   Rscript 05_neighbourhood_genes.R [config.R]
# Flags   FORCE_REBUILD, REPLOT_ONLY (config). No downsampling. ASCII only.
# =============================================================================

source("helpers.R")
load_config()
require_param("INDEX_STATES", "NEIGHBOUR_STATES")
OUT <- module_dir(RESULTS_DIR, "05_neighbourhood_genes")

samples   <- read_samples(INPUT_DIR, SAMPLES)
label_map <- read_label_map(INPUT_DIR)
miss <- setdiff(c(INDEX_STATES, NEIGHBOUR_STATES), label_map$label)
if (length(miss)) stop("labels not in label_map.csv: ", paste(miss, collapse = ", "))
INDEX_COMP <- unique(label_map$compartment[match(INDEX_STATES, label_map$label)])
NB_COMP    <- unique(label_map$compartment[match(NEIGHBOUR_STATES, label_map$label)])
if (length(INDEX_COMP) != 1) stop("INDEX_STATES must all belong to one compartment")
if (any(INDEX_STATES %in% NEIGHBOUR_STATES)) stop("a label is both index and neighbour")
msg("[1] %d samples; index = %s (%s); neighbours = %s (%s)", nrow(samples),
    paste(INDEX_STATES, collapse = ","), INDEX_COMP,
    paste(NEIGHBOUR_STATES, collapse = ","), paste(NB_COMP, collapse = ","))

params <- list(samples = samples$sample_id, label_map = label_map,
               INDEX_STATES = INDEX_STATES, NEIGHBOUR_STATES = NEIGHBOUR_STATES,
               RADIUS_UM = RADIUS_UM, MIN_MEDIAN_DEPTH = MIN_MEDIAN_DEPTH,
               MIN_INDEX_CELLS = MIN_INDEX_CELLS,
               MIN_NEIGHBOUR_CELLS = MIN_NEIGHBOUR_CELLS,
               components = file.info(file.path(RESULTS_DIR, "01_geometry_qc",
                                                "components.rds"))$mtime)
CK <- file.path(OUT, "cache", "05_checkpoint.rds")
built <- read_checkpoint(CK, params, FORCE_REBUILD, REPLOT_ONLY)

# =============================================================================
# Build: per-sample partial correlations and in-situ detection rates
# =============================================================================
if (is.null(built)) {
  rho <- list(); det <- list(); det_iso <- list(); det_nb <- list(); incl <- list()
  for (sid in samples$sample_id) {
    s  <- load_sample(sid, INPUT_DIR, label_map, counts = TRUE)
    cl <- s$cells
    med <- stats::median(cl$n_count[!is.na(cl$label)], na.rm = TRUE)
    row <- data.frame(sample_id = sid, median_n_count_labelled = med,
                      passes_depth_floor = is.finite(med) && med >= MIN_MEDIAN_DEPTH,
                      n_index = NA_integer_, n_neighbour_cells = NA_integer_,
                      included = FALSE, mean_neighbours = NA_real_,
                      frac_index_with_neighbour = NA_real_,
                      mean_all_cells_in_radius = NA_real_, stringsAsFactors = FALSE)
    if (!row$passes_depth_floor) {
      incl[[sid]] <- row; msg("    %s excluded (median depth %.0f)", sid, med); next
    }
    comp   <- read_components(RESULTS_DIR, sid, cl$cell_id)
    tissue <- comp > 0
    idx <- which(tissue & cl$label %in% INDEX_STATES & is.finite(cl$n_count) & cl$n_count > 0)
    nbr <- which(tissue & cl$label %in% NEIGHBOUR_STATES)
    all_i <- which(tissue)
    row$n_index <- length(idx); row$n_neighbour_cells <- length(nbr)
    if (length(idx) < MIN_INDEX_CELLS || length(nbr) < MIN_NEIGHBOUR_CELLS) {
      incl[[sid]] <- row; msg("    %s excluded (too few index or neighbour cells)", sid); next
    }
    x <- cl$x_um; y <- cl$y_um
    n_nb  <- count_within(x[idx], y[idx], x[nbr], y[nbr], RADIUS_UM)
    n_all <- count_within(x[idx], y[idx], x[all_i], y[all_i], RADIUS_UM) - 1L
    row$included <- TRUE
    row$mean_neighbours <- mean(n_nb); row$frac_index_with_neighbour <- mean(n_nb > 0)
    row$mean_all_cells_in_radius <- mean(n_all)
    incl[[sid]] <- row

    Y <- Matrix::t(s$counts[, idx, drop = FALSE])
    Y <- methods::as(Y, "CsparseMatrix"); Y@x <- rep(1, length(Y@x))
    Z <- cbind(1, log(cl$n_count[idx]), log1p(n_all))
    r <- partial_cor(Y, log1p(n_nb), Z)
    d <- Matrix::colMeans(Y)
    r[d < MIN_DETECT_FRAC] <- NA_real_
    rho[[sid]] <- stats::setNames(r, rownames(s$counts))
    det[[sid]] <- stats::setNames(d, rownames(s$counts))

    # in-situ spillover reference: isolated index cells vs neighbour compartment
    iso <- idx[n_nb == 0]
    nbc <- which(tissue & cl$compartment %in% NB_COMP)
    det_rate <- function(cells) {
      if (!length(cells)) return(stats::setNames(rep(NA_real_, nrow(s$counts)),
                                                 rownames(s$counts)))
      m <- s$counts[, cells, drop = FALSE]
      stats::setNames(Matrix::rowMeans(m > 0), rownames(s$counts))
    }
    det_iso[[sid]] <- det_rate(iso)
    det_nb[[sid]]  <- det_rate(nbc)
    msg("    %s  index=%d  neighbours=%d  isolated index=%d", sid, length(idx),
        length(nbr), length(iso))
  }
  to_mat <- function(lst) {
    if (!length(lst)) return(NULL)
    g <- unique(unlist(lapply(lst, names)))
    m <- sapply(lst, function(v) v[g]); rownames(m) <- g; m
  }
  built <- list(rho = to_mat(rho), det = to_mat(det), det_iso = to_mat(det_iso),
                det_nb = to_mat(det_nb), incl = do.call(rbind, incl))
  write_checkpoint(CK, params, built)
}
if (is.null(built$rho)) stop("no sample passed the inclusion rules")

# =============================================================================
# Cohort summary, spillover filter, selection
# =============================================================================
R <- built$rho
res <- data.frame(
  gene = rownames(R),
  n_samples = rowSums(is.finite(R)),
  median_cor = apply(R, 1, function(v) if (any(is.finite(v))) stats::median(v, na.rm = TRUE) else NA_real_),
  frac_positive = apply(R, 1, function(v) { v <- v[is.finite(v)]; if (length(v)) mean(v > 0) else NA_real_ }),
  median_detect_index = apply(built$det, 1, stats::median, na.rm = TRUE),
  p_wilcoxon = apply(R, 1, function(v) {
    v <- v[is.finite(v)]
    if (length(v) < MIN_SAMPLES_GENE) return(NA_real_)
    suppressWarnings(stats::wilcox.test(v, mu = 0, exact = FALSE)$p.value) }),
  stringsAsFactors = FALSE)
res$q_BH <- stats::p.adjust(res$p_wilcoxon, method = "BH")

if (!is.null(SPILL_REFERENCE)) {
  ref <- utils::read.csv(SPILL_REFERENCE, stringsAsFactors = FALSE, check.names = FALSE)
  names(ref)[1] <- "gene"
  if (!INDEX_COMP %in% names(ref))
    stop("SPILL_REFERENCE has no column for the index compartment: ", INDEX_COMP)
  others <- setdiff(names(ref)[-1], INDEX_COMP)
  rr <- ref[match(res$gene, ref$gene), , drop = FALSE]
  om <- as.matrix(rr[, others, drop = FALSE])
  res$spill_filter <- "reference"
  res$spill_ratio  <- rr[[INDEX_COMP]] / apply(om, 1, max)
  res$spill_top_other <- others[apply(om, 1, function(v) if (all(is.na(v))) NA else which.max(v))]
  res$in_reference <- !is.na(rr$gene)
} else {
  iso <- apply(built$det_iso, 1, stats::median, na.rm = TRUE)[res$gene]
  nbd <- apply(built$det_nb,  1, stats::median, na.rm = TRUE)[res$gene]
  res$spill_filter <- "in situ"
  res$detect_isolated_index <- unname(iso)
  res$detect_neighbour_compartment <- unname(nbd)
  res$spill_ratio <- unname(ifelse(nbd > 0, iso / nbd, ifelse(iso > 0, Inf, NA_real_)))
}
res$pass_spill    <- !is.na(res$spill_ratio) & res$spill_ratio >= SPILL_RATIO_MIN
res$pass_samples  <- res$n_samples >= MIN_SAMPLES_GENE
res$pass_assoc    <- res$pass_samples & !is.na(res$q_BH) & res$q_BH < FDR_MAX &
                     res$median_cor > 0
res$pass_consist  <- res$pass_assoc & res$frac_positive >= MIN_FRAC_POS
res$passes_all    <- res$pass_consist & res$pass_spill
res$removed_by_spillover <- res$pass_consist & !res$pass_spill
res <- res[order(!res$passes_all, -res$median_cor), ]
res$rank <- NA_integer_
res$rank[res$passes_all] <- seq_len(sum(res$passes_all))
res$selected <- !is.na(res$rank) & res$rank <= N_TOP
rownames(res) <- NULL

funnel <- data.frame(
  step = c("genes", sprintf("tested in >= %d samples", MIN_SAMPLES_GENE),
           sprintf("BH q < %g and median correlation > 0", FDR_MAX),
           sprintf("and positive in >= %.0f%% of samples", 100 * MIN_FRAC_POS),
           sprintf("and spillover ratio >= %g (%s)", SPILL_RATIO_MIN, res$spill_filter[1]),
           sprintf("selected (top %d)", N_TOP), "removed by the spillover filter"),
  n_genes = c(nrow(res), sum(res$pass_samples), sum(res$pass_assoc),
              sum(res$pass_consist), sum(res$passes_all), sum(res$selected),
              sum(res$removed_by_spillover)))
msg("[2] Selection funnel:"); print(funnel, row.names = FALSE)

per_sample <- data.frame(gene = rownames(R), R, check.names = FALSE)
write_xlsx(list(sample_inclusion = built$incl, per_sample_cor = per_sample,
                gene_summary = res, selected_genes = res[res$selected, ],
                spillover_removed = res[res$removed_by_spillover, ],
                selection_funnel = funnel),
           file.path(OUT, "05_neighbourhood_genes.xlsx"))

# =============================================================================
# Figures
# =============================================================================
V <- res[is.finite(res$q_BH) & is.finite(res$median_cor), ]
V$class <- factor(ifelse(V$selected, "selected",
                  ifelse(V$removed_by_spillover, "removed by spillover filter", "other")),
                  levels = c("selected", "removed by spillover filter", "other"))
V <- V[order(V$class, decreasing = TRUE), ]
lbl <- V[V$class != "other", ]
p_vol <- ggplot(V, aes(median_cor, -log10(q_BH), colour = class)) +
  geom_hline(yintercept = -log10(FDR_MAX), linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_point(size = 1.6, alpha = 0.8) +
  geom_text(data = lbl, aes(label = gene), size = BASE_SIZE / 5, vjust = -0.7,
            show.legend = FALSE, check_overlap = TRUE) +
  scale_colour_manual(values = c(selected = "#DC050C",
                                 `removed by spillover filter` = "#1965B0",
                                 other = "grey70"), name = NULL) +
  labs(x = "median per-sample partial correlation",
       y = expression(-log[10] ~ "BH q"),
       caption = sprintf("Index: %s. Neighbours: %s within %d um. Spillover filter: %s.",
                         paste(INDEX_STATES, collapse = ", "),
                         paste(NEIGHBOUR_STATES, collapse = ", "), RADIUS_UM,
                         res$spill_filter[1])) +
  theme_niche(BASE_SIZE) + theme(legend.position = "bottom")
f <- file.path(OUT, "05_gene_volcano.pdf")
open_pdf(f, width = 7.5, height = 6.5); print(p_vol); close_pdf(f)

sel <- res$gene[res$selected]
if (length(sel)) {
  Hm <- R[sel, , drop = FALSE]
  H <- data.frame(gene = rep(sel, ncol(Hm)), sample_id = rep(colnames(Hm), each = length(sel)),
                  cor = as.vector(Hm), stringsAsFactors = FALSE)
  H$gene <- factor(H$gene, levels = rev(sel))
  lim <- max(abs(H$cor), na.rm = TRUE)
  p_hm <- ggplot(H, aes(sample_id, gene, fill = cor)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = HEAT_LOW, mid = HEAT_MID, high = HEAT_HIGH, midpoint = 0,
                         limits = c(-lim, lim), na.value = "grey85",
                         name = "partial\ncorrelation") +
    labs(x = "sample", y = NULL) +
    theme_niche(BASE_SIZE) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  f <- file.path(OUT, "05_selected_heatmap.pdf")
  open_pdf(f, width = 0.35 * ncol(Hm) + 4, height = 0.28 * length(sel) + 2.2)
  print(p_hm); close_pdf(f)
} else msg("    no gene selected; heatmap not drawn")
write_session_info(OUT)
msg("[3] Done.")
