#!/usr/bin/env Rscript
# =============================================================================
# 04_exposure_contrast.R -- are B cells closer to A cells that express gene g?
#
# Question
#   Module 02 asks which B state sits near which A state. This module asks
#   whether that proximity is greater around the A cells that express
#   EXPOSURE_GENE (for example a ligand) than around those that do not.
#
# Design (A-centred)
#   Focal cells are labelled in-tissue cells of an A state, split into gene+
#   (count > 0) and gene- arms. 
#
# Readouts, per focal A state and B state
#   PRIMARY    distance from each focal cell to its nearest cell of the B
#              state (radius-free, censored at MAX_D_UM). Reported as
#              gene- minus gene+, so POSITIVE = B closer to gene+ A cells.
#   COMPANION  number of cells of the B state within RADIUS_UM. Reported as
#              gene+ minus gene-, so POSITIVE = more B around gene+ A cells.
#   Sign convention: across both readouts, positive = the gene+ effect.
#
# Depth 
#   Gene detection rises with a cell's own sequencing depth, so "gene+"
#   partly means "deeply sequenced", and depth is often spatially clustered.
#   To overcoem this limitation, focal cells are split
#   into N_DEPTH_STRATA within-sample quantile strata of n_count, and the
#   arm difference is the size-weighted mean of within-stratum differences.
#   The unstratified difference and the median depth of each arm are reported
#   alongside, so the reader can see whether depth is a cofounder.
#
# Cohort inference (the sample is the replicate unit)
#   Wilcoxon signed-rank of per-sample differences against 0, HL estimate and
#   95% CI. BH within each readout (stratified and unstratified separately),
#   for pre-specified pairs (PRIMARY_PAIRS) and the exploratory family
#   separately. Stars use the depth-stratified q.
#
# Inclusion
#   Samples: median n_count of labelled cells >= MIN_MEDIAN_DEPTH (the
#   exposure is expression-based). Cells: labelled, in tissue, n_count > 0.
#   A state: >= MIN_STATE_CELLS cells and >= MIN_ARM cells in each arm.
#   B state: >= MIN_STATE_CELLS cells. Distance readout: >= MIN_ARM uncensored
#   focal cells per arm. Cohort test: >= MIN_SAMPLES samples.
#
# Reads   INPUT_DIR (cells.rds, counts.rds), RESULTS_DIR/01_geometry_qc/
# Writes  RESULTS_DIR/04_exposure_contrast/
#           04_exposure_contrast.xlsx  sample_inclusion, distance_per_sample,
#                                      distance_summary, count_per_sample,
#                                      count_summary, state_diagnostics
#           04_effect_distance.pdf   per-sample points + HL estimate and CI
#           04_effect_count.pdf
#           04_heatmap_distance.pdf  every A x B pair
#           04_heatmap_count.pdf
#           cache/04_checkpoint.rds
#
# Usage   Rscript 04_exposure_contrast.R [config.R]
# Flags   FORCE_REBUILD, REPLOT_ONLY (config). No downsampling. ASCII only.
# =============================================================================

source("helpers.R")
load_config()
require_param("COMP_A", "COMP_B", "EXPOSURE_GENE")
OUT <- module_dir(RESULTS_DIR, "04_exposure_contrast")

samples   <- read_samples(INPUT_DIR, SAMPLES)
label_map <- read_label_map(INPUT_DIR)
if (COMP_A == COMP_B) stop("COMP_A and COMP_B must be different compartments")
ST_A <- if (is.null(STATES_A)) label_map$label[label_map$compartment == COMP_A] else STATES_A
ST_B <- if (is.null(STATES_B)) label_map$label[label_map$compartment == COMP_B] else STATES_B
msg("[1] %d samples; exposure gene %s; %d A states x %d B states",
    nrow(samples), EXPOSURE_GENE, length(ST_A), length(ST_B))

params <- list(samples = samples$sample_id, label_map = label_map,
               COMP_A = COMP_A, COMP_B = COMP_B, ST_A = ST_A, ST_B = ST_B,
               EXPOSURE_GENE = EXPOSURE_GENE, RADIUS_UM = RADIUS_UM,
               MAX_D_UM = MAX_D_UM, MIN_MEDIAN_DEPTH = MIN_MEDIAN_DEPTH,
               MIN_STATE_CELLS = MIN_STATE_CELLS, MIN_ARM = MIN_ARM,
               N_DEPTH_STRATA = N_DEPTH_STRATA,
               components = file.info(file.path(RESULTS_DIR, "01_geometry_qc",
                                                "components.rds"))$mtime)
CK <- file.path(OUT, "cache", "04_checkpoint.rds")
built <- read_checkpoint(CK, params, FORCE_REBUILD, REPLOT_ONLY)

# =============================================================================
# Build
# =============================================================================
if (is.null(built)) {
  crow <- list(); drow <- list(); incl <- list(); diag <- list()
  for (sid in samples$sample_id) {
    s  <- load_sample(sid, INPUT_DIR, label_map, counts = TRUE)
    cl <- s$cells
    med <- stats::median(cl$n_count[!is.na(cl$label)], na.rm = TRUE)
    has_gene <- EXPOSURE_GENE %in% rownames(s$counts)
    passes <- is.finite(med) && med >= MIN_MEDIAN_DEPTH && has_gene
    incl[[sid]] <- data.frame(sample_id = sid, median_n_count_labelled = med,
                              gene_on_panel = has_gene, included = passes,
                              stringsAsFactors = FALSE)
    if (!passes) {
      msg("    %s excluded (median depth %.0f%s)", sid, med,
          if (!has_gene) ", gene not on panel" else "")
      next
    }
    comp <- read_components(RESULTS_DIR, sid, cl$cell_id)
    keep <- !is.na(cl$label) & comp > 0 & is.finite(cl$n_count) & cl$n_count > 0
    g    <- as.numeric(s$counts[EXPOSURE_GENE, keep])
    x <- cl$x_um[keep]; y <- cl$y_um[keep]; lab <- cl$label[keep]
    cmp <- cl$compartment[keep]; nc <- cl$n_count[keep]
    tb <- table(lab)

    for (a in ST_A) {
      fi <- which(lab == a)
      if (!length(fi)) next
      diag[[length(diag) + 1]] <- data.frame(
        sample_id = sid, state_A = a, n_cells = length(fi),
        n_gene_pos = sum(g[fi] > 0), stringsAsFactors = FALSE)
    }

    bi <- which(cmp == COMP_B)
    b_present <- ST_B[ST_B %in% names(tb)]
    b_present <- b_present[tb[b_present] >= MIN_STATE_CELLS]
    if (!length(bi) || !length(b_present)) next

    for (a in ST_A) {
      if (is.na(tb[a]) || tb[a] < MIN_STATE_CELLS) next
      fi  <- which(lab == a)
      arm <- g[fi] > 0
      if (sum(arm) < MIN_ARM || sum(!arm) < MIN_ARM) next
      strata <- depth_strata(nc[fi], N_DEPTH_STRATA)

      pr <- bipartite_pairs(x[fi], y[fi], x[bi], y[bi], RADIUS_UM)
      for (b in b_present) {
        hit <- pr[lab[bi][pr$t] == b, , drop = FALSE]
        v <- tabulate(hit$q, nbins = length(fi))
        crow[[length(crow) + 1]] <- data.frame(
          sample_id = sid, state_A = a, state_B = b,
          n_focal_pos = sum(arm), n_focal_neg = sum(!arm),
          diff_unstratified = mean(v[arm]) - mean(v[!arm]),
          diff_depth_stratified = strat_diff(v, arm, strata),
          median_depth_pos = stats::median(nc[fi][arm]),
          median_depth_neg = stats::median(nc[fi][!arm]),
          stringsAsFactors = FALSE)

        tS <- bi[lab[bi] == b]
        dS <- nearest_dist(x[fi], y[fi], x[tS], y[tS], MAX_D_UM)
        ok <- is.finite(dS)
        if (sum(arm & ok) < MIN_ARM || sum((!arm) & ok) < MIN_ARM) next
        drow[[length(drow) + 1]] <- data.frame(
          sample_id = sid, state_A = a, state_B = b,
          n_focal_pos = sum(arm & ok), n_focal_neg = sum((!arm) & ok),
          diff_unstratified = mean(dS[!arm & ok]) - mean(dS[arm & ok]),
          diff_depth_stratified = -strat_diff(dS[ok], arm[ok], strata[ok]),
          frac_censored = mean(!ok),
          median_depth_pos = stats::median(nc[fi][arm]),
          median_depth_neg = stats::median(nc[fi][!arm]),
          stringsAsFactors = FALSE)
      }
    }
    msg("    %s  labelled=%d", sid, length(x))
  }
  built <- list(count = do.call(rbind, crow), dist = do.call(rbind, drow),
                incl = do.call(rbind, incl), diag = do.call(rbind, diag))
  write_checkpoint(CK, params, built)
}

# =============================================================================
# Cohort statistics
# =============================================================================
summarise_readout <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  k <- paste(df$state_A, df$state_B, sep = "|")
  out <- do.call(rbind, lapply(unique(k), function(kk) {
    d <- df[k == kk, ]
    rbind(
      { h <- hl_test(d$diff_depth_stratified, MIN_SAMPLES)
        if (is.null(h)) NULL else cbind(data.frame(state_A = d$state_A[1],
          state_B = d$state_B[1], adjustment = "depth_stratified",
          stringsAsFactors = FALSE), h) },
      { h <- hl_test(d$diff_unstratified, MIN_SAMPLES)
        if (is.null(h)) NULL else cbind(data.frame(state_A = d$state_A[1],
          state_B = d$state_B[1], adjustment = "unstratified",
          stringsAsFactors = FALSE), h) })
  }))
  if (is.null(out)) return(NULL)
  out$is_primary_pair <- is_primary_pair(out$state_A, out$state_B, PRIMARY_PAIRS)
  out$family <- ifelse(out$is_primary_pair, "pre-specified", "exploratory")
  out$q_BH <- bh_by_family(out$p_raw, out$is_primary_pair, out$adjustment)
  out$star <- stars(out$q_BH)
  out[order(out$adjustment != "depth_stratified", !out$is_primary_pair, -out$estimate), ]
}
SD <- summarise_readout(built$dist)
SC <- summarise_readout(built$count)

DG <- built$diag
state_diag <- if (is.null(DG)) NULL else do.call(rbind, lapply(split(DG, DG$state_A), function(z)
  data.frame(state_A = z$state_A[1], n_samples_present = nrow(z),
             n_samples_ge_min_cells = sum(z$n_cells >= MIN_STATE_CELLS),
             n_samples_testable = sum(z$n_cells >= MIN_STATE_CELLS &
                                        z$n_gene_pos >= MIN_ARM &
                                        (z$n_cells - z$n_gene_pos) >= MIN_ARM),
             median_cells = stats::median(z$n_cells),
             median_gene_pos = stats::median(z$n_gene_pos),
             stringsAsFactors = FALSE)))

if (!is.null(SD)) {
  msg("[2] Distance readout, depth-stratified (positive um = B closer to gene+ A):")
  print(SD[SD$adjustment == "depth_stratified",
           c("state_A", "state_B", "family", "n_samples", "estimate", "ci_low",
             "ci_high", "q_BH")], row.names = FALSE)
}

write_xlsx(list(sample_inclusion = built$incl, distance_per_sample = built$dist,
                distance_summary = SD, count_per_sample = built$count,
                count_summary = SC, state_diagnostics = state_diag),
           file.path(OUT, "04_exposure_contrast.xlsx"))

# =============================================================================
# Figures
# =============================================================================
effect_plot <- function(per, summ, ylab) {
  if (is.null(summ)) return(NULL)
  s <- summ[summ$adjustment == "depth_stratified", ]
  d <- per; d$value <- d$diff_depth_stratified
  d <- d[paste(d$state_A, d$state_B) %in% paste(s$state_A, s$state_B), ]
  for (z in list("d", "s")) {
    w <- get(z); w$state_A <- factor(w$state_A, levels = ST_A)
    w$state_B <- factor(w$state_B, levels = ST_B); assign(z, w)
  }
  # stars sit above every point of their column, not on top of them
  span <- diff(range(c(d$value, s$ci_high, s$ci_low), na.rm = TRUE))
  top  <- tapply(d$value, paste(d$state_A, d$state_B), max, na.rm = TRUE)
  s$y_star <- pmax(s$ci_high, top[paste(s$state_A, s$state_B)], na.rm = TRUE) + 0.06 * span
  ggplot(d, aes(state_B, value)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_jitter(width = 0.12, height = 0, size = 1.1, alpha = 0.55, colour = "grey35") +
    geom_errorbar(data = s, aes(y = estimate, ymin = ci_low, ymax = ci_high),
                  width = 0.25, linewidth = 0.6) +
    geom_point(data = s, aes(y = estimate), size = 2.6, shape = 18) +
    geom_text(data = s, aes(y = y_star, label = star), size = BASE_SIZE / 3.5) +
    facet_wrap(~ state_A, labeller = labeller(state_A = function(v) paste(COMP_A, "state", v))) +
    labs(x = paste0(COMP_B, " state"), y = ylab,
         caption = paste0("Points: samples. Diamond and bar: HL estimate and 95% CI. ",
                          "Depth-stratified; stars: BH q.")) +
    theme_niche(BASE_SIZE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
heat_plot <- function(summ, title) {
  if (is.null(summ)) return(NULL)
  h <- summ[summ$adjustment == "depth_stratified", ]
  h$state_A <- factor(h$state_A, levels = rev(ST_A))
  h$state_B <- factor(h$state_B, levels = ST_B)
  lim <- max(abs(h$estimate), na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1
  ggplot(h, aes(state_B, state_A, fill = estimate)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = star), size = BASE_SIZE / 3.5, vjust = 0.75) +
    coord_fixed() +
    scale_fill_gradient2(low = HEAT_LOW, mid = HEAT_MID, high = HEAT_HIGH,
                         midpoint = 0, limits = c(-lim, lim), name = title) +
    labs(x = paste0(COMP_B, " state"), y = paste0(COMP_A, " state")) +
    theme_niche(BASE_SIZE) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}
lab_d <- sprintf("nearest B: %s- minus %s+ (um)", EXPOSURE_GENE, EXPOSURE_GENE)
lab_c <- sprintf("B within %d um: %s+ minus %s-", RADIUS_UM, EXPOSURE_GENE, EXPOSURE_GENE)
w <- 2.2 * min(length(ST_A), 3) + 0.6 * length(ST_B) + 2
for (spec in list(list("04_effect_distance.pdf", effect_plot(built$dist, SD, lab_d)),
                  list("04_effect_count.pdf",    effect_plot(built$count, SC, lab_c)),
                  list("04_heatmap_distance.pdf", heat_plot(SD, "HL estimate\n(um)")),
                  list("04_heatmap_count.pdf",    heat_plot(SC, "HL estimate\n(cells)")))) {
  if (is.null(spec[[2]])) next
  f <- file.path(OUT, spec[[1]])
  open_pdf(f, width = if (grepl("heatmap", spec[[1]])) 1.1 * length(ST_B) + 4 else w,
           height = if (grepl("heatmap", spec[[1]])) 0.7 * length(ST_A) + 3 else 5)
  print(spec[[2]])
  close_pdf(f)
}
write_session_info(OUT)
msg("[3] Done.")
