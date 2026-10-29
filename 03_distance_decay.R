#!/usr/bin/env Rscript
# =============================================================================
# 03_distance_decay.R -- over what range does the association persist
#
# Aim
#   This analysis adds the spatial RANGE of an association (at what distance it 
#   vanishes) and its SHAPE (a graded decline, or a step at a niche boundary) 
#   to the results in module 02.
#
# Statistic
#   For every labelled in-tissue cell, the distance to the nearest cell of A
#   state a, censored at MAX_D_UM, binned at BIN_WIDTH_UM. Per bin and B
#   state: the count of B-state cells, their fraction of all labelled cells,
#   and their fraction of compartment-B cells (the primary readout, which
#   conditions out B abundance and tissue edges).
#
# Null
#   B labels shuffled among compartment-B cells within tissue component,
#   positions fixed (distances unchanged; each permutation is a relabel and
#   re-tabulate). Gives the expected count per bin and a 95% envelope.
#   Effect per sample and bin: log2((obs + 1) / (exp + 1)).
#
# Cohort inference (the sample is the replicate unit)
#   PER BIN     Wilcoxon signed-rank of per-sample log2 obs/exp against 0,
#               HL estimate and 95% CI; BH across bins within an A x B pair.
#   SLOPE       per sample, least-squares slope of log2 obs/exp on bin
#               midpoint (per 100 um), over at least MIN_BINS_FIT bins;
#               Wilcoxon of the slopes against 0, BH across pairs. 
#   CROSSING    per sample, the interpolated distance at which log2 obs/exp
#               first crosses zero from above; median and IQR across samples.
#
# Star rule, fixed before the run (STAR_DIRECTION)
#   The primary readout is a proportion of compartment-B cells, so it is
#   compositional: enrichment in some bins arithmetically forces depletion in
#   others. Bins are starred only when BH q < 0.05 AND the estimate lies in
#   the hypothesised direction: "positive" (default; attraction, e.g. a
#   niche), "negative" (avoidance, e.g. exclusion) or "both" (exploratory;
#   the figure caption then says the two directions are not independent).
#   Workbooks report every bin whatever is starred.
#
# Inclusion
#   A state and B state: at least MIN_STATE_CELLS cells in the sample. Bins
#   with fewer than MIN_BIN_B compartment-B cells are dropped for that
#   sample. At least MIN_SAMPLES samples per bin. No depth floor (labels only).
#   Caveat: samples where A is abundant have a compressed distance
#   range, so far bins are supported by fewer samples; n per bin is exported.
#
# Reads   INPUT_DIR (cells.rds), RESULTS_DIR/01_geometry_qc/components.rds
# Writes  RESULTS_DIR/03_distance_decay/
#           03_distance_decay.xlsx  sample_info, per_sample_bin, bin_tests,
#                                   slopes_per_sample, slope_tests, crossing
#           03_decay_primary_pairs.pdf   pre-specified pairs (if any)
#           03_decay_all_pairs.pdf       every A x B pair
#           cache/03_checkpoint.rds
#
# Usage   Rscript 03_distance_decay.R [config.R]
# Flags   FORCE_REBUILD, REPLOT_ONLY (config). No downsampling. ASCII only.
# =============================================================================

source("helpers.R")
load_config()
require_param("COMP_A", "COMP_B")
OUT <- module_dir(RESULTS_DIR, "03_distance_decay")

samples   <- read_samples(INPUT_DIR, SAMPLES)
label_map <- read_label_map(INPUT_DIR)
if (COMP_A == COMP_B) stop("COMP_A and COMP_B must be different compartments")
ST_A <- if (is.null(STATES_A)) label_map$label[label_map$compartment == COMP_A] else STATES_A
ST_B <- if (is.null(STATES_B)) label_map$label[label_map$compartment == COMP_B] else STATES_B
directional_stars(0, 1, STAR_DIRECTION)   # validates STAR_DIRECTION early

BREAKS  <- seq(0, MAX_D_UM, by = BIN_WIDTH_UM)
NBIN    <- length(BREAKS) - 1L
BIN_MID <- BREAKS[-length(BREAKS)] + BIN_WIDTH_UM / 2
msg("[1] %d samples; %d A states x %d B states; %d bins of %g um",
    nrow(samples), length(ST_A), length(ST_B), NBIN, BIN_WIDTH_UM)

params <- list(samples = samples$sample_id, label_map = label_map,
               COMP_A = COMP_A, COMP_B = COMP_B, ST_A = ST_A, ST_B = ST_B,
               BREAKS = BREAKS, N_PERM_DECAY = N_PERM_DECAY, SEED = SEED,
               MIN_STATE_CELLS = MIN_STATE_CELLS, MIN_BIN_B = MIN_BIN_B,
               components = file.info(file.path(RESULTS_DIR, "01_geometry_qc",
                                                "components.rds"))$mtime)
CK <- file.path(OUT, "cache", "03_checkpoint.rds")
built <- read_checkpoint(CK, params, FORCE_REBUILD, REPLOT_ONLY)

# =============================================================================
# Build
# =============================================================================
if (is.null(built)) {
  rows <- list(); info <- list()
  for (si in seq_len(nrow(samples))) {
    sid <- samples$sample_id[si]
    set.seed(SEED + si)
    cl   <- load_sample(sid, INPUT_DIR, label_map)$cells
    comp <- read_components(RESULTS_DIR, sid, cl$cell_id)
    keep <- !is.na(cl$label) & comp > 0
    x <- cl$x_um[keep]; y <- cl$y_um[keep]; lab <- cl$label[keep]
    isB <- cl$compartment[keep] == COMP_B; comp <- comp[keep]
    info[[sid]] <- data.frame(sample_id = sid, n_labelled_in_tissue = length(x),
                              n_B = sum(isB), n_components = length(unique(comp)),
                              stringsAsFactors = FALSE)
    tb <- table(lab)
    b_present <- ST_B[ST_B %in% names(tb)]
    b_present <- b_present[tb[b_present] >= MIN_STATE_CELLS]
    if (length(b_present) == 0) next

    b_idx  <- which(isB); b_lab <- lab[b_idx]
    idx_by_comp <- strata_index(comp[b_idx])

    for (a in ST_A) {
      if (is.na(tb[a]) || tb[a] < MIN_STATE_CELLS) next
      ti <- which(lab == a)
      d_all <- nearest_dist(x, y, x[ti], y[ti], MAX_D_UM)
      bin <- findInterval(d_all, BREAKS, rightmost.closed = TRUE)
      bin[is.na(d_all) | bin < 1 | bin > NBIN] <- NA_integer_
      n_all_bin <- tabulate(bin[!is.na(bin)], nbins = NBIN)
      bin_b <- bin[b_idx]; valid <- !is.na(bin_b)
      n_b_bin <- tabulate(bin_b[valid], nbins = NBIN)

      tab_states <- function(labs) {
        m <- matrix(0L, NBIN, length(b_present), dimnames = list(NULL, b_present))
        for (j in seq_along(b_present)) {
          sel <- valid & labs == b_present[j]
          if (any(sel)) m[, j] <- tabulate(bin_b[sel], nbins = NBIN)
        }
        m
      }
      obs  <- tab_states(b_lab)
      perm <- array(0, dim = c(N_PERM_DECAY, NBIN, length(b_present)))
      for (k in seq_len(N_PERM_DECAY))
        perm[k, , ] <- tab_states(shuffle_within(b_lab, idx_by_comp))
      ex <- apply(perm, c(2, 3), mean)
      lo <- apply(perm, c(2, 3), stats::quantile, probs = 0.025)
      hi <- apply(perm, c(2, 3), stats::quantile, probs = 0.975)
      ex <- matrix(ex, NBIN); lo <- matrix(lo, NBIN); hi <- matrix(hi, NBIN)

      ok_bin <- n_b_bin >= MIN_BIN_B
      if (!any(ok_bin)) next
      for (j in seq_along(b_present)) {
        rows[[length(rows) + 1]] <- data.frame(
          sample_id = sid, state_A = a, state_B = b_present[j],
          bin_mid_um = BIN_MID[ok_bin],
          n_target = obs[ok_bin, j], n_B = n_b_bin[ok_bin],
          n_all_cells = n_all_bin[ok_bin],
          prop_of_B = obs[ok_bin, j] / n_b_bin[ok_bin],
          prop_of_all = obs[ok_bin, j] / pmax(n_all_bin[ok_bin], 1),
          exp_target = ex[ok_bin, j], exp_low = lo[ok_bin, j], exp_high = hi[ok_bin, j],
          exp_prop_of_B = ex[ok_bin, j] / n_b_bin[ok_bin],
          log2_obs_exp = log2((obs[ok_bin, j] + 1) / (ex[ok_bin, j] + 1)),
          stringsAsFactors = FALSE)
      }
    }
    msg("    %s  labelled=%d  B=%d", sid, length(x), sum(isB))
  }
  built <- list(D = do.call(rbind, rows), info = do.call(rbind, info))
  write_checkpoint(CK, params, built)
}
D  <- built$D
PI <- built$info
if (is.null(D) || nrow(D) == 0) stop("no A x B pair passed the inclusion rules")
D$is_primary_pair <- is_primary_pair(D$state_A, D$state_B, PRIMARY_PAIRS)

# =============================================================================
# Cohort statistics
# =============================================================================
kb <- paste(D$state_A, D$state_B, D$bin_mid_um, sep = "|")
BIN <- do.call(rbind, lapply(unique(kb), function(k) {
  d <- D[kb == k, ]
  h <- hl_test(d$log2_obs_exp, MIN_SAMPLES); if (is.null(h)) return(NULL)
  cbind(data.frame(state_A = d$state_A[1], state_B = d$state_B[1],
                   bin_mid_um = d$bin_mid_um[1],
                   is_primary_pair = d$is_primary_pair[1],
                   median_prop_of_B = stats::median(d$prop_of_B),
                   median_prop_of_all = stats::median(d$prop_of_all),
                   median_n_target = stats::median(d$n_target),
                   stringsAsFactors = FALSE), h)
}))
if (is.null(BIN)) stop("no bin reached MIN_SAMPLES samples")
BIN$q_BH <- NA_real_
cmb <- paste(BIN$state_A, BIN$state_B)
for (c1 in unique(cmb)) BIN$q_BH[cmb == c1] <- stats::p.adjust(BIN$p_raw[cmb == c1], "BH")
BIN$star <- directional_stars(BIN$estimate, BIN$q_BH, STAR_DIRECTION)
BIN$star_direction <- STAR_DIRECTION
BIN <- BIN[order(BIN$state_A, BIN$state_B, BIN$bin_mid_um), ]

ks <- paste(D$state_A, D$state_B, D$sample_id, sep = "|")
SLO <- do.call(rbind, lapply(unique(ks), function(k) {
  d <- D[ks == k, ]; d <- d[is.finite(d$log2_obs_exp), ]
  if (nrow(d) < MIN_BINS_FIT) return(NULL)
  fit <- stats::lm(log2_obs_exp ~ bin_mid_um, data = d)
  data.frame(state_A = d$state_A[1], state_B = d$state_B[1], sample_id = d$sample_id[1],
             slope_per_100um = 100 * unname(stats::coef(fit)[2]), n_bins = nrow(d),
             crossing_um = crossing(d$bin_mid_um, d$log2_obs_exp),
             stringsAsFactors = FALSE)
}))

kc <- paste(SLO$state_A, SLO$state_B, sep = "|")
SLOPE <- do.call(rbind, lapply(unique(kc), function(k) {
  d <- SLO[kc == k, ]
  h <- hl_test(d$slope_per_100um, MIN_SAMPLES); if (is.null(h)) return(NULL)
  cbind(data.frame(state_A = d$state_A[1], state_B = d$state_B[1],
                   is_primary_pair = is_primary_pair(d$state_A[1], d$state_B[1],
                                                     PRIMARY_PAIRS),
                   stringsAsFactors = FALSE), h)
}))
if (!is.null(SLOPE)) {
  SLOPE$q_BH <- bh_by_family(SLOPE$p_raw, SLOPE$is_primary_pair)
  SLOPE$family <- ifelse(SLOPE$is_primary_pair, "pre-specified", "exploratory")
  SLOPE$trend <- ifelse(is.na(SLOPE$q_BH) | SLOPE$q_BH >= 0.05, "n.s.",
                        ifelse(SLOPE$estimate < 0, "decay", "increase"))
}

CROSS <- do.call(rbind, lapply(unique(kc), function(k) {
  d <- SLO[kc == k, ]; v <- d$crossing_um[is.finite(d$crossing_um)]
  if (length(v) < MIN_SAMPLES) return(NULL)
  data.frame(state_A = d$state_A[1], state_B = d$state_B[1],
             n_samples_crossing = length(v), median_crossing_um = stats::median(v),
             iqr_low = unname(stats::quantile(v, 0.25)),
             iqr_high = unname(stats::quantile(v, 0.75)), stringsAsFactors = FALSE)
}))

if (!is.null(SLOPE)) {
  msg("[2] Decay slopes (log2 obs/exp per 100 um):")
  print(SLOPE[order(SLOPE$estimate), c("state_A", "state_B", "n_samples", "estimate",
                                       "ci_low", "ci_high", "q_BH", "trend")],
        row.names = FALSE)
}

write_xlsx(list(sample_info = PI, per_sample_bin = D, bin_tests = BIN,
                slopes_per_sample = SLO, slope_tests = SLOPE, crossing = CROSS),
           file.path(OUT, "03_distance_decay.xlsx"))

# =============================================================================
# Figures: one line per facet (A state rows x B state columns), CI band,
# stars on bins that pass the star rule, colour = trend class
# =============================================================================
TREND_COL <- c(decay = "#DC050C", increase = "#1965B0", `n.s.` = "grey45")
decay_plot <- function(B, caption) {
  B$trend <- if (is.null(SLOPE)) "n.s." else
    SLOPE$trend[match(paste(B$state_A, B$state_B), paste(SLOPE$state_A, SLOPE$state_B))]
  B$trend[is.na(B$trend)] <- "n.s."
  B$trend <- factor(B$trend, levels = names(TREND_COL))
  B$state_A <- factor(B$state_A, levels = ST_A)
  B$state_B <- factor(B$state_B, levels = ST_B)
  yl <- range(c(B$ci_low, B$ci_high), na.rm = TRUE)
  mk <- B[B$star != "", ]
  if (nrow(mk)) mk$y_star <- mk$ci_high + 0.04 * diff(yl)
  p <- ggplot(B, aes(bin_mid_um, estimate)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey30") +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = trend), alpha = 0.2,
                colour = NA) +
    geom_line(aes(colour = trend), linewidth = 0.6) +
    geom_point(aes(colour = trend), size = 1.4)
  if (nrow(mk)) p <- p + geom_text(data = mk, aes(y = y_star, label = star),
                                   vjust = 0, size = BASE_SIZE / 3.5) +
    expand_limits(y = max(mk$y_star) + 0.08 * diff(yl))
  p + facet_grid(state_A ~ state_B) +
    scale_colour_manual(values = TREND_COL, drop = FALSE, name = "slope") +
    scale_fill_manual(values = TREND_COL, drop = FALSE, name = "slope") +
    labs(x = expression(paste("distance to nearest A-state cell (", mu, "m)")),
         y = expression(log[2] ~ "observed / expected"), caption = caption) +
    theme_niche(BASE_SIZE) +
    theme(panel.border = element_rect(fill = NA, linewidth = 0.4))
}
cap <- sprintf(paste0("HL estimate and 95%% CI across samples. Stars: BH q < 0.05 ",
                      "across bins, %s direction only%s."),
               STAR_DIRECTION,
               if (STAR_DIRECTION == "both")
                 "; enrichment and depletion are not independent (compositional)" else "")

f <- file.path(OUT, "03_decay_all_pairs.pdf")
open_pdf(f, width = 2.6 * length(ST_B) + 2.5, height = 2.4 * length(ST_A) + 1.8)
print(decay_plot(BIN, cap))
close_pdf(f)

if (any(BIN$is_primary_pair)) {
  BP <- BIN[BIN$is_primary_pair, ]
  f <- file.path(OUT, "03_decay_primary_pairs.pdf")
  open_pdf(f, width = 2.8 * length(unique(BP$state_B)) + 2.5,
           height = 2.6 * length(unique(BP$state_A)) + 1.8)
  print(decay_plot(BP, cap) + facet_grid(state_A ~ state_B, drop = TRUE))
  close_pdf(f)
}
write_session_info(OUT)
msg("[3] Done.")
