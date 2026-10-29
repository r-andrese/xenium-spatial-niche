#!/usr/bin/env Rscript
# =============================================================================
# 02_colocalisation.R -- which B states sit next to which A states
#
# Question
#   Which B-cell state is closer to which A-cell state than expected given each 
#   compartment’s existing tissue location? Comparing states to states, rather 
#   than to empty space, controls for the tissue’s broad architecture.
#
# Statistic
#   Per sample, one fixed-radius neighbour graph at RADIUS_UM over labelled
#   cells in tissue (component > 0, from 01). For each A state x B state pair,
#   the observed statistic is the number of edges joining the pair. Effect
#   size: log2((obs + 1) / (exp + 1)), with exp the permutation mean; plus a
#   z-score and a two-sided empirical p. The graph is built once per sample;
#   each permutation only relabels nodes and re-tabulates edges.
#
# Nulls (labels never move between tissue components)
#   PRIMARY  within compartment. A labels are shuffled among compartment-A
#            cells and B labels among compartment-B cells, positions fixed.
#            Asks whether THIS B state is preferentially next to THIS A state,
#            given where each compartment sits.
#   GLOBAL   all labels shuffled among all labelled cells of the component.
#            This also recovers the coarse architecture, so most A x B pairs
#            read as enriched.
#
# Cohort inference (the sample is the replicate unit)
#   Per pair, Wilcoxon signed-rank of per-sample log2 enrichment against 0,
#   with Hodges-Lehmann estimate and 95% CI, median, IQR and fraction of
#   samples positive; at least MIN_SAMPLES samples. BH within each null,
#   separately for the pre-specified pairs (PRIMARY_PAIRS) and the
#   exploratory family (all other pairs).
#
# Inclusion
#   A pair contributes in a sample only if both states have at least
#   MIN_STATE_CELLS labelled in-tissue cells there. No depth floor: this
#   module uses labels only. Unlabelled cells occupy space but are never
#   shuffled or counted.
#
# Reads   INPUT_DIR (cells.rds), RESULTS_DIR/01_geometry_qc/components.rds
# Writes  RESULTS_DIR/02_colocalisation/
#           02_colocalisation.xlsx  sample_info, per_sample, cohort_summary
#           02_heatmaps.pdf         within-compartment and global null, side
#                                   by side
#           02_per_sample.pdf       per-sample enrichment, one page per A state
#           cache/02_checkpoint.rds
#
# Usage   Rscript 02_colocalisation.R [config.R]
# Flags   FORCE_REBUILD, REPLOT_ONLY (config). No downsampling. ASCII only.
# =============================================================================

source("helpers.R")
load_config()
require_param("COMP_A", "COMP_B")
OUT <- module_dir(RESULTS_DIR, "02_colocalisation")

samples   <- read_samples(INPUT_DIR, SAMPLES)
label_map <- read_label_map(INPUT_DIR)
for (cc in c(COMP_A, COMP_B))
  if (!cc %in% label_map$compartment) stop("compartment not in label_map.csv: ", cc)
if (COMP_A == COMP_B) stop("COMP_A and COMP_B must be different compartments")

ST_A <- if (is.null(STATES_A)) label_map$label[label_map$compartment == COMP_A] else STATES_A
ST_B <- if (is.null(STATES_B)) label_map$label[label_map$compartment == COMP_B] else STATES_B
# A listed state assigned to another compartment is a configuration error; a
# listed state absent from label_map.csv is dropped with a message.
check_states <- function(st, comp) {
  wrong <- st[st %in% label_map$label &
              label_map$compartment[match(st, label_map$label)] != comp]
  if (length(wrong)) stop("states not in compartment ", comp, " in label_map.csv: ",
                          paste(wrong, collapse = ", "))
  absent <- setdiff(st, label_map$label)
  if (length(absent)) msg("    states not in label_map.csv, dropped: %s", paste(absent, collapse = ", "))
  intersect(st, label_map$label)
}
ST_A <- check_states(ST_A, COMP_A)
ST_B <- check_states(ST_B, COMP_B)
nA <- length(ST_A); nB <- length(ST_B)
msg("[1] %d samples; A = %s (%d states); B = %s (%d states)",
    nrow(samples), COMP_A, nA, COMP_B, nB)

params <- list(samples = samples$sample_id, label_map = label_map,
               COMP_A = COMP_A, COMP_B = COMP_B, ST_A = ST_A, ST_B = ST_B,
               RADIUS_UM = RADIUS_UM, N_PERM_PRIMARY = N_PERM_PRIMARY,
               N_PERM_GLOBAL = N_PERM_GLOBAL, SEED = SEED,
               components = file.info(file.path(RESULTS_DIR, "01_geometry_qc",
                                                "components.rds"))$mtime)
CK <- file.path(OUT, "cache", "02_checkpoint.rds")
built <- read_checkpoint(CK, params, FORCE_REBUILD, REPLOT_ONLY)

# =============================================================================
# Build: per-sample permutation tests
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
    cmp <- cl$compartment[keep]; comp <- comp[keep]

    info[[sid]] <- data.frame(
      sample_id = sid, n_cells = nrow(cl), n_labelled_in_tissue = sum(keep),
      n_A = sum(cmp == COMP_A), n_B = sum(cmp == COMP_B),
      n_components = length(unique(comp)),
      median_n_count_labelled = stats::median(cl$n_count[keep]),
      pct_unlabelled = 100 * mean(is.na(cl$label)), stringsAsFactors = FALSE)

    code <- integer(length(lab))
    ha <- match(lab, ST_A); code[!is.na(ha)] <- ha[!is.na(ha)]
    hb <- match(lab, ST_B); code[!is.na(hb)] <- nA + hb[!is.na(hb)]
    n_state <- table(factor(lab, levels = c(ST_A, ST_B)))

    E <- radius_edges(x, y, RADIUS_UM)
    msg("    %s  labelled=%d  edges=%d", sid, length(x), nrow(E))
    if (nrow(E) == 0) next

    # A -> "A", B -> "B", every other compartment pooled as "O"
    grp <- ifelse(cmp == COMP_A, "A", ifelse(cmp == COMP_B, "B", "O"))
    r_pri <- edge_permutation(code, strata_index(paste(comp, grp, sep = ".")),
                              E[, 1], E[, 2], nA, nB, N_PERM_PRIMARY)
    r_glo <- edge_permutation(code, strata_index(as.character(comp)),
                              E[, 1], E[, 2], nA, nB, N_PERM_GLOBAL)

    for (a in seq_len(nA)) for (b in seq_len(nB)) {
      na_ <- as.integer(n_state[ST_A[a]]); nb_ <- as.integer(n_state[ST_B[b]])
      for (mdl in c("within_compartment", "global")) {
        rr <- if (mdl == "within_compartment") r_pri else r_glo
        rows[[length(rows) + 1]] <- data.frame(
          sample_id = sid, null_model = mdl, state_A = ST_A[a], state_B = ST_B[b],
          n_A = na_, n_B = nb_,
          included = na_ >= MIN_STATE_CELLS & nb_ >= MIN_STATE_CELLS,
          obs_edges = rr$obs[a, b], exp_edges = rr$exp[a, b], sd_edges = rr$sd[a, b],
          log2_enrichment = rr$log2fc[a, b], z = rr$z[a, b], p_perm = rr$p[a, b],
          mean_B_neighbours_per_A = if (na_ > 0) rr$obs[a, b] / na_ else NA_real_,
          stringsAsFactors = FALSE)
      }
    }
  }
  built <- list(R = do.call(rbind, rows), info = do.call(rbind, info))
  write_checkpoint(CK, params, built)
}

R  <- built$R
PI <- built$info
R$is_primary_pair <- is_primary_pair(R$state_A, R$state_B, PRIMARY_PAIRS)
RI <- R[R$included, ]

# =============================================================================
# Cohort summary
# =============================================================================
key <- paste(RI$null_model, RI$state_A, RI$state_B, sep = "|")
M <- do.call(rbind, lapply(unique(key), function(k) {
  d <- RI[key == k, ]
  h <- hl_test(d$log2_enrichment, MIN_SAMPLES)
  if (is.null(h)) return(NULL)
  v <- d$log2_enrichment[is.finite(d$log2_enrichment)]
  cbind(data.frame(null_model = d$null_model[1], state_A = d$state_A[1],
                   state_B = d$state_B[1], is_primary_pair = d$is_primary_pair[1],
                   stringsAsFactors = FALSE),
        h, iqr_low = unname(stats::quantile(v, 0.25)),
        iqr_high = unname(stats::quantile(v, 0.75)))
}))
if (is.null(M)) stop("no pair reached MIN_SAMPLES samples; nothing to summarise")
M$family <- ifelse(M$is_primary_pair, "pre-specified", "exploratory")
M$q_BH   <- bh_by_family(M$p_raw, M$is_primary_pair, M$null_model)
M$star   <- stars(M$q_BH)
M <- M[order(M$null_model != "within_compartment", !M$is_primary_pair,
             -M$median), ]
rownames(M) <- NULL

msg("[2] Within-compartment null, pre-specified pairs:")
pp <- M[M$null_model == "within_compartment" & M$is_primary_pair, ]
if (nrow(pp)) print(pp[, c("state_A", "state_B", "n_samples", "median",
                           "frac_positive", "p_raw", "q_BH")], row.names = FALSE)

write_xlsx(list(sample_info = PI, per_sample = R, cohort_summary = M),
           file.path(OUT, "02_colocalisation.xlsx"))

# =============================================================================
# Figures
# =============================================================================
H <- M
H$null_lab <- factor(ifelse(H$null_model == "within_compartment",
                            "within-compartment null", "global null"),
                     levels = c("within-compartment null", "global null"))
H$state_A <- factor(H$state_A, levels = rev(ST_A))
H$state_B <- factor(H$state_B, levels = ST_B)
lim <- max(abs(H$median), na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1
p_heat <- ggplot(H, aes(state_B, state_A, fill = median)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = star), size = BASE_SIZE / 3.5, vjust = 0.75) +
  facet_wrap(~ null_lab) +
  coord_fixed() +
  scale_fill_gradient2(low = HEAT_LOW, mid = HEAT_MID, high = HEAT_HIGH,
                       midpoint = 0, limits = c(-lim, lim),
                       name = "median\nlog2 enrichment") +
  labs(x = paste0(COMP_B, " state"), y = paste0(COMP_A, " state"),
       caption = sprintf("%d um neighbourhood; stars: BH q within null and family",
                         RADIUS_UM)) +
  theme_niche(BASE_SIZE) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
f <- file.path(OUT, "02_heatmaps.pdf")
open_pdf(f, width = max(8, 1.2 * nB * 2 + 4), height = max(4, 0.6 * nA + 3))
print(p_heat)
close_pdf(f)

f <- file.path(OUT, "02_per_sample.pdf")
open_pdf(f, width = max(7, 0.9 * nB + 4), height = 5.5)
for (a in ST_A) {
  d <- RI[RI$state_A == a, ]
  if (nrow(d) == 0) next
  d$state_B <- factor(d$state_B, levels = ST_B)
  d$null_lab <- ifelse(d$null_model == "within_compartment",
                       "within-compartment null", "global null")
  print(ggplot(d, aes(state_B, log2_enrichment)) +
          geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
          geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.4) +
          geom_jitter(width = 0.12, height = 0, size = 1.2, alpha = 0.7) +
          facet_wrap(~ null_lab) +
          labs(x = paste0(COMP_B, " state"), y = "log2 enrichment vs null",
               title = paste0(COMP_A, " state ", a), subtitle = "one point per sample") +
          theme_niche(BASE_SIZE) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}
close_pdf(f)
write_session_info(OUT)
msg("[3] Done.")
