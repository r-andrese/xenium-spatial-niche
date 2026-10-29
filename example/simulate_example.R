#!/usr/bin/env Rscript
# =============================================================================
# example/simulate_example.R -- simulated cohorts with known truth
#
# Writes two cohorts of 12 samples in the standard input format, plus the
# truth that example/validate_example.R checks the toolkit against.
#
#   example/input/       PLANTED cohort
#   example/input_null/  NULL cohort: same tissue architecture, no planted
#                        state-level structure
#   example/truth/       truth.rds, spill_reference.csv
#
# Tissue (both cohorts). Each sample is a 2 x 2 mm section at ~3,500 cells per
# mm2, with six tumour nests in a stroma. Compartment A (states A1, A2) and
# compartment B (states B1, B2, B3) are both enriched in stroma, so the
# architecture alone makes every A state look close to every B state. Sample
# S08 is two pieces 300 um apart plus a clump of debris and scattered single
# cells. Sequencing depth (n_count, median ~150) follows a smooth spatial
# field, so deeply sequenced cells are spatially clustered. 10% of cells are
# unlabelled.
#
# Planted cohort only:
#   - Niche. Extra B1 cells are placed around A1 cells (2D Gaussian offset,
#     sd 30 um), so B1 is enriched near A1 and the enrichment fades within
#     about 50-75 um. A1 is 12% of compartment A, so A1 cells are sparse
#     enough for the decay to be visible.
#   - Exposure. 35% of A1 cells truly express gene GENEX and receive more B1
#     around them (mean 2.0 versus 0.5 extra B1 cells).
#   - Depth confound. GENEX is detected in every other cell at a rate that
#     rises with depth, and B3 is more frequent where depth is high. A2 cells
#     detected GENEX+ therefore sit in deeper regions with more B3, with no
#     biological link: a stratified analysis should remove this.
#   - Genes. TRUE01-10: detection in A cells rises with the number of B1
#     cells within 50 um. BLEED01-10: B-specific genes that leak into A cells
#     in proportion to the B cells within 15 um (segmentation spillover).
#
# Null cohort: B1 placed at random among B cells (no extra B1 around A1), B
# states independent of depth, GENEX detection depth-only in all cells, and
# no TRUE-gene response. BLEED genes leak as in the planted cohort.
#
# Run from the repository root:  Rscript example/simulate_example.R
# ASCII only.
# =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(dbscan) })

OUT        <- "example"
N_SAMPLES  <- 12
TWO_PIECES <- "S08"
DENSITY    <- 3500e-6          # cells per um2
SEED       <- 20261027

GENE_SETS <- list(
  exposure = "GENEX",
  true     = sprintf("TRUE%02d", 1:10),
  bleed    = sprintf("BLEED%02d", 1:10),
  bg_A     = sprintf("BGA%02d", 1:30),
  bg_B     = sprintf("BGB%02d", 1:30),
  bg_T     = sprintf("BGT%02d", 1:20),
  bg_U     = sprintf("BGU%02d", 1:20)
)
GENES <- unlist(GENE_SETS, use.names = FALSE)

LABEL_MAP <- data.frame(
  label       = c("A1", "A2", "B1", "B2", "B3", "Tumour", "Tcell", "Endo"),
  compartment = c("A",  "A",  "B",  "B",  "B",  "tumour", "other", "other"),
  stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
in_rects <- function(x, y, rects) {
  ok <- rep(FALSE, length(x))
  for (r in rects) ok <- ok | (x >= r[1] & x <= r[2] & y >= r[3] & y <= r[4])
  ok
}

depth_field <- function(x, y, centres, amp) {
  f <- numeric(length(x))
  for (k in seq_len(nrow(centres)))
    f <- f + amp[k] * exp(-((x - centres[k, 1])^2 + (y - centres[k, 2])^2) / (2 * 300^2))
  f
}

simulate_sample <- function(sid, planted) {
  two  <- sid == TWO_PIECES
  rects <- if (two) list(c(0, 1000, 0, 2000), c(1300, 2300, 0, 2000))
           else     list(c(0, 2000, 0, 2000))
  xr <- range(unlist(lapply(rects, `[`, 1:2)))

  # base cells
  xs <- ys <- numeric(0)
  for (r in rects) {
    n <- rpois(1, (r[2] - r[1]) * (r[4] - r[3]) * DENSITY)
    xs <- c(xs, runif(n, r[1], r[2])); ys <- c(ys, runif(n, r[3], r[4]))
  }
  n <- length(xs)

  # tumour nests
  cen <- cbind(runif(6, xr[1], xr[2]), runif(6, 0, 2000)); rad <- runif(6, 150, 280)
  nest <- rep(FALSE, n)
  for (k in 1:6) nest <- nest | ((xs - cen[k, 1])^2 + (ys - cen[k, 2])^2 <= rad[k]^2)

  # compartments by region
  u <- runif(n)
  comp <- ifelse(nest,
                 ifelse(u < 0.85, "tumour", ifelse(u < 0.88, "A", ifelse(u < 0.95, "B", "other"))),
                 ifelse(u < 0.05, "tumour", ifelse(u < 0.45, "A", ifelse(u < 0.65, "B", "other"))))

  # depth field (standardised over the section)
  dc  <- cbind(runif(5, xr[1], xr[2]), runif(5, 0, 2000)); amp <- sample(c(-1, 1), 5, TRUE)
  fld <- function(x, y) depth_field(x, y, dc, amp)
  f0  <- fld(xs, ys); fm <- mean(f0); fs <- sd(f0)
  zf  <- function(x, y) (fld(x, y) - fm) / fs
  f   <- zf(xs, ys)

  # labels
  lab <- character(n)
  lab[comp == "tumour"] <- "Tumour"
  lab[comp == "other"]  <- ifelse(runif(sum(comp == "other")) < 0.5, "Tcell", "Endo")
  ia <- which(comp == "A")
  lab[ia] <- ifelse(runif(length(ia)) < 0.12, "A1", "A2")
  ib <- which(comp == "B")
  if (planted) {
    p3 <- plogis(-0.3 + 1.2 * f[ib])
    v  <- runif(length(ib))
    lab[ib] <- ifelse(v < 0.2, "B1", ifelse(runif(length(ib)) < p3, "B3", "B2"))
  } else {
    v <- runif(length(ib))
    lab[ib] <- ifelse(v < 0.2, "B1", ifelse(v < 0.6, "B2", "B3"))
  }

  # true GENEX status of A1 cells
  gtrue <- rep(FALSE, n)
  if (planted) {
    a1 <- which(lab == "A1")
    gtrue[a1] <- runif(length(a1)) < 0.35
  }

  # planted niche: extra B1 around A1
  if (planted) {
    a1  <- which(lab == "A1")
    noff <- rpois(length(a1), ifelse(gtrue[a1], 2.0, 0.5))
    par <- rep(a1, noff)
    ox <- xs[par] + rnorm(length(par), 0, 30); oy <- ys[par] + rnorm(length(par), 0, 30)
    ok <- in_rects(ox, oy, rects)
    ox <- ox[ok]; oy <- oy[ok]
    m <- length(ox)
    xs <- c(xs, ox); ys <- c(ys, oy)
    comp <- c(comp, rep("B", m)); lab <- c(lab, rep("B1", m))
    gtrue <- c(gtrue, rep(FALSE, m)); f <- c(f, zf(ox, oy))
    nest_new <- rep(FALSE, m)
    for (k in 1:6) nest_new <- nest_new | ((ox - cen[k, 1])^2 + (oy - cen[k, 2])^2 <= rad[k]^2)
    nest <- c(nest, nest_new)
  }

  # debris (S08 only)
  debris <- rep(FALSE, length(xs))
  if (two) {
    cx <- 3000 + rnorm(40, 0, 25); cy <- 1000 + rnorm(40, 0, 25)
    sx <- runif(10, 2600, 4000);   sy <- runif(10, 2200, 3000)
    dx <- c(cx, sx); dy <- c(cy, sy); m <- length(dx)
    xs <- c(xs, dx); ys <- c(ys, dy)
    comp <- c(comp, rep("tumour", m)); lab <- c(lab, rep("Tumour", m))
    gtrue <- c(gtrue, rep(FALSE, m)); f <- c(f, rep(0, m)); nest <- c(nest, rep(FALSE, m))
    debris <- c(debris, rep(TRUE, m))
  }
  n <- length(xs)

  # depth
  nc <- pmax(5L, rnbinom(n, mu = 150 * exp(0.5 * f), size = 8))
  sc <- nc / 150

  # neighbourhood quantities used by the gene model (true labels)
  isA <- comp == "A"; isB <- comp == "B"
  nB1_50 <- count_n(xs, ys, which(lab == "B1"), 50)
  nB_15  <- count_n(xs, ys, which(isB), 15)

  # expression rates (expected counts per cell)
  rate <- matrix(0, nrow = length(GENES), ncol = n, dimnames = list(GENES, NULL))
  # GENEX
  rate["GENEX", ] <- 0.12 * sc^1.5
  if (planted) rate["GENEX", gtrue] <- 1.0 * sc[gtrue]
  # TRUE genes
  for (g in GENE_SETS$true) {
    r <- 0.005 * sc
    r[isA] <- (0.04 + if (planted) 0.05 * nB1_50[isA] else 0) * sc[isA]
    rate[g, ] <- r
  }
  # BLEED genes
  for (g in GENE_SETS$bleed) {
    r <- 0.004 * sc
    r[isB] <- 0.5 * sc[isB]
    r[isA] <- (0.004 + 0.04 * nB_15[isA]) * sc[isA]
    rate[g, ] <- r
  }
  # background
  for (g in GENE_SETS$bg_A) rate[g, ] <- ifelse(isA, 0.3, 0.01) * sc
  for (g in GENE_SETS$bg_B) rate[g, ] <- ifelse(isB, 0.3, 0.01) * sc
  for (g in GENE_SETS$bg_T) rate[g, ] <- ifelse(comp == "tumour", 0.3, 0.01) * sc
  for (g in GENE_SETS$bg_U) rate[g, ] <- 0.15 * sc

  # counts, sparse
  ii <- jj <- vv <- vector("list", length(GENES))
  for (k in seq_along(GENES)) {
    v <- rpois(n, rate[k, ]); nz <- which(v > 0)
    ii[[k]] <- rep(k, length(nz)); jj[[k]] <- nz; vv[[k]] <- v[nz]
  }
  cell_id <- sprintf("%s_c%06d", sid, seq_len(n))
  counts <- sparseMatrix(i = unlist(ii), j = unlist(jj), x = as.numeric(unlist(vv)),
                         dims = c(length(GENES), n), dimnames = list(GENES, cell_id))

  # 10% unlabelled
  true_label <- lab
  lab[runif(n) < 0.10] <- NA_character_

  list(cells = data.frame(cell_id = cell_id, x_um = xs, y_um = ys, label = lab,
                          n_count = as.integer(nc), stringsAsFactors = FALSE),
       counts = counts,
       truth = data.frame(cell_id = cell_id, true_label = true_label, debris = debris,
                          genex_true = gtrue, in_nest = nest, stringsAsFactors = FALSE))
}

count_n <- function(x, y, targets, r) {
  if (length(targets) == 0) return(integer(length(x)))
  nn <- frNN(cbind(x[targets], y[targets]), eps = r, query = cbind(x, y), sort = FALSE)
  cnt <- lengths(nn$id)
  # a target is not its own neighbour
  self <- seq_along(x) %in% targets
  as.integer(cnt - self)
}

write_cohort <- function(dir, planted, seed) {
  set.seed(seed)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  sids <- sprintf("S%02d", seq_len(N_SAMPLES))
  truth <- list()
  for (sid in sids) {
    s <- simulate_sample(sid, planted)
    dir.create(file.path(dir, sid), showWarnings = FALSE)
    saveRDS(s$cells,  file.path(dir, sid, "cells.rds"))
    saveRDS(s$counts, file.path(dir, sid, "counts.rds"))
    truth[[sid]] <- s$truth
    message(sprintf("  %s %s  cells=%d  labelled=%d", basename(dir), sid,
                    nrow(s$cells), sum(!is.na(s$cells$label))))
  }
  utils::write.csv(data.frame(sample_id = sids,
                              group = rep(c("X", "Y"), length.out = N_SAMPLES)),
                   file.path(dir, "samples.csv"), row.names = FALSE)
  utils::write.csv(LABEL_MAP, file.path(dir, "label_map.csv"), row.names = FALSE)
  truth
}

# Reference table of mean expression per compartment (what a matched scRNA-seq
# atlas would give): rates of an isolated cell at the median depth.
spill_reference <- function() {
  ref <- data.frame(gene = GENES, A = 0, B = 0, tumour = 0, other = 0)
  set_row <- function(g, A, B, Tm, O) {
    k <- ref$gene %in% g
    ref$A[k] <<- A; ref$B[k] <<- B; ref$tumour[k] <<- Tm; ref$other[k] <<- O
  }
  set_row("GENEX", 0.20, 0.12, 0.12, 0.12)
  set_row(GENE_SETS$true,  0.06, 0.005, 0.005, 0.005)
  set_row(GENE_SETS$bleed, 0.004, 0.5, 0.004, 0.004)
  set_row(GENE_SETS$bg_A, 0.3, 0.01, 0.01, 0.01)
  set_row(GENE_SETS$bg_B, 0.01, 0.3, 0.01, 0.01)
  set_row(GENE_SETS$bg_T, 0.01, 0.01, 0.3, 0.01)
  set_row(GENE_SETS$bg_U, 0.15, 0.15, 0.15, 0.15)
  ref
}

# -----------------------------------------------------------------------------
message("Simulating the planted cohort")
tp <- write_cohort(file.path(OUT, "input"), planted = TRUE, seed = SEED)
message("Simulating the null cohort")
tn <- write_cohort(file.path(OUT, "input_null"), planted = FALSE, seed = SEED + 1)

dir.create(file.path(OUT, "truth"), showWarnings = FALSE)
saveRDS(list(planted = tp, null = tn, gene_sets = GENE_SETS,
             two_pieces = TWO_PIECES),
        file.path(OUT, "truth", "truth.rds"))
utils::write.csv(spill_reference(), file.path(OUT, "truth", "spill_reference.csv"),
                 row.names = FALSE)
message("Done: example/input, example/input_null, example/truth")
