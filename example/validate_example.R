#!/usr/bin/env Rscript
# =============================================================================
# example/validate_example.R -- does the toolkit recover known truth?
#
# Runs every module on the two simulated cohorts written by
# simulate_example.R and checks each result against what was planted:
#
#   PLANTED cohort   effects that should be found, and a depth artefact that
#                    should disappear after stratification
#   NULL cohort      no state-level structure: nothing should be found
#                    except the architecture, which the global null (02)
#                    reports and the within-compartment null does not
#   ADAPTERS         the simulated samples written out as Xenium Ranger
#                    output folders and as a merged Seurat object (pixel
#                    coordinates), imported back, and compared with the
#                    original standard input
#
# Writes  example/validation_report.xlsx   one row per check: expected,
#                                          observed, pass
# Exit    status 1 if any check fails.
#
# Run from the repository root:
#   Rscript example/simulate_example.R
#   Rscript example/validate_example.R
# Takes about 10 minutes on one core (module 02 dominates).
# ASCII only.
# =============================================================================

suppressPackageStartupMessages({ library(openxlsx); library(Matrix) })
CHECKS <- list()
check <- function(cohort, module, what, expected, observed, pass) {
  CHECKS[[length(CHECKS) + 1]] <<- data.frame(
    cohort = cohort, module = module, check = what, expected = expected,
    observed = observed, pass = isTRUE(pass), stringsAsFactors = FALSE)
  message(sprintf("  [%s] %-8s %-3s %s", if (isTRUE(pass)) "PASS" else "FAIL",
                  cohort, module, what))
}
run <- function(script, cfg) {
  message(sprintf("--- Rscript %s %s", script, cfg))
  st <- system2("Rscript", c(script, cfg), stdout = FALSE, stderr = FALSE)
  if (st != 0) stop(script, " failed with ", cfg, "; rerun it by hand to see the error")
}
xl <- function(f, sheet) read.xlsx(f, sheet = sheet)
fmt <- function(v, d = 3) paste(formatC(v, digits = d, format = "g"), collapse = "; ")

truth <- readRDS("example/truth/truth.rds")
gs <- truth$gene_sets
background <- c(gs$bg_A, gs$bg_B, gs$bg_T, gs$bg_U)

# =============================================================================
# Run the modules
# =============================================================================
for (cfg in c("example/config_example.R", "example/config_null.R"))
  for (s in c("01_geometry_qc.R", "02_colocalisation.R", "03_distance_decay.R",
              "04_exposure_contrast.R", "05_neighbourhood_genes.R")) run(s, cfg)
dir.create("example/results_reference/01_geometry_qc", recursive = TRUE, showWarnings = FALSE)
file.copy("example/results/01_geometry_qc/components.rds",
          "example/results_reference/01_geometry_qc/components.rds", overwrite = TRUE)
run("05_neighbourhood_genes.R", "example/config_reference.R")

# =============================================================================
# Module checks
# =============================================================================
module_checks <- function(cohort, res) {
  planted <- cohort == "planted"
  tr <- truth[[cohort]]

  # ---- 01 ------------------------------------------------------------------
  ps <- xl(file.path(res, "01_geometry_qc/01_geometry_qc.xlsx"), "per_sample")
  cmp <- readRDS(file.path(res, "01_geometry_qc/components.rds"))
  two <- truth$two_pieces
  check(cohort, "01", "two-piece sample has 2 components; all others 1",
        "2 / 1", sprintf("%d / %s", ps$n_components[ps$sample_id == two],
                         paste(unique(ps$n_components[ps$sample_id != two]), collapse = ",")),
        ps$n_components[ps$sample_id == two] == 2 && all(ps$n_components[ps$sample_id != two] == 1))
  t2 <- tr[[two]]; c2 <- cmp[cmp$sample_id == two, ]
  zero <- c2$cell_id[c2$component == 0]
  check(cohort, "01", "excluded cells are exactly the planted debris",
        sprintf("%d debris cells", sum(t2$debris)),
        sprintf("%d excluded, %d of them debris", length(zero),
                sum(zero %in% t2$cell_id[t2$debris])),
        setequal(zero, t2$cell_id[t2$debris]))

  # ---- 02 ------------------------------------------------------------------
  m <- xl(file.path(res, "02_colocalisation/02_colocalisation.xlsx"), "cohort_summary")
  wc <- m[m$null_model == "within_compartment", ]; gl <- m[m$null_model == "global", ]
  if (planted) {
    r <- wc[wc$state_A == "A1" & wc$state_B == "B1", ]
    check(cohort, "02", "A1 x B1 enriched, within-compartment null",
          "median > 0, q < 0.05, largest median",
          sprintf("median %.3f, q %.3g, rank %d", r$median, r$q_BH,
                  rank(-wc$median)[wc$state_A == "A1" & wc$state_B == "B1"]),
          r$median > 0 && r$q_BH < 0.05 && r$median == max(wc$median))
  } else {
    check(cohort, "02", "no pair significant, within-compartment null",
          "0 pairs with q < 0.05", sprintf("%d pairs", sum(wc$q_BH < 0.05)),
          sum(wc$q_BH < 0.05) == 0)
  }
  check(cohort, "02", "architecture detected by the global null",
        "every pair median > 0 and q < 0.05",
        sprintf("%d of %d pairs", sum(gl$median > 0 & gl$q_BH < 0.05), nrow(gl)),
        all(gl$median > 0 & gl$q_BH < 0.05))

  # ---- 03 ------------------------------------------------------------------
  f3 <- file.path(res, "03_distance_decay/03_distance_decay.xlsx")
  sl <- xl(f3, "slope_tests"); bt <- xl(f3, "bin_tests"); cr <- xl(f3, "crossing")
  if (planted) {
    s <- sl[sl$state_A == "A1" & sl$state_B == "B1", ]
    b <- bt[bt$state_A == "A1" & bt$state_B == "B1", ]
    c1 <- cr[cr$state_A == "A1" & cr$state_B == "B1", ]
    check(cohort, "03", "A1 x B1 decays with distance",
          "slope < 0, q < 0.05", sprintf("slope %.3f per 100 um, q %.3g", s$estimate, s$q_BH),
          s$estimate < 0 && s$q_BH < 0.05)
    check(cohort, "03", "A1 x B1 starred in the nearest bin, not beyond 75 um",
          "star at 12.5 um; none past 75 um",
          sprintf("starred bins: %s um", fmt(b$bin_mid_um[b$star != ""])),
          b$star[b$bin_mid_um == min(b$bin_mid_um)] != "" &&
            !any(b$star[b$bin_mid_um > 75] != ""))
    check(cohort, "03", "A1 x B1 enrichment crosses zero within the planted range",
          "median crossing 25-75 um", sprintf("%.1f um", c1$median_crossing_um),
          nrow(c1) == 1 && c1$median_crossing_um >= 25 && c1$median_crossing_um <= 75)
  } else {
    check(cohort, "03", "no decay slope significant", "0 pairs with q < 0.05",
          sprintf("%d pairs", sum(sl$q_BH < 0.05)), sum(sl$q_BH < 0.05) == 0)
    check(cohort, "03", "no bin starred", "0 starred bins",
          sprintf("%d bins", sum(bt$star != "")), sum(bt$star != "") == 0)
  }

  # ---- 04 ------------------------------------------------------------------
  f4 <- file.path(res, "04_exposure_contrast/04_exposure_contrast.xlsx")
  sd <- xl(f4, "distance_summary"); sc <- xl(f4, "count_summary")
  if (planted) {
    r <- sd[sd$state_A == "A1" & sd$state_B == "B1" & sd$adjustment == "depth_stratified", ]
    check(cohort, "04", "B1 closer to GENEX+ A1 cells (depth-stratified)",
          "estimate > 0 um, q < 0.05", sprintf("%.2f um, q %.3g", r$estimate, r$q_BH),
          r$estimate > 0 && r$q_BH < 0.05)
    u <- sd[sd$state_A == "A2" & sd$state_B == "B3" & sd$adjustment == "unstratified", ]
    s <- sd[sd$state_A == "A2" & sd$state_B == "B3" & sd$adjustment == "depth_stratified", ]
    check(cohort, "04", "depth artefact A2 x B3 present without stratification",
          "unstratified estimate > 0, q < 0.05", sprintf("%.2f um, q %.3g", u$estimate, u$q_BH),
          u$estimate > 0 && u$q_BH < 0.05)
    check(cohort, "04", "depth artefact A2 x B3 removed by stratification",
          "depth-stratified q >= 0.05", sprintf("%.2f um, q %.3g", s$estimate, s$q_BH),
          s$q_BH >= 0.05)
  } else {
    n_sig <- sum(sd$q_BH[sd$adjustment == "depth_stratified"] < 0.05) +
             sum(sc$q_BH[sc$adjustment == "depth_stratified"] < 0.05)
    check(cohort, "04", "no depth-stratified effect significant",
          "0 tests with q < 0.05", sprintf("%d tests", n_sig), n_sig == 0)
  }

  # ---- 05 ------------------------------------------------------------------
  g05 <- function(dir, label) {
    gsum <- xl(file.path(dir, "05_neighbourhood_genes/05_neighbourhood_genes.xlsx"), "gene_summary")
    sel <- gsum$gene[gsum$selected]; rem <- gsum$gene[gsum$removed_by_spillover]
    if (planted) {
      check(cohort, "05", paste0("true neighbourhood genes selected (", label, ")"),
            "10 of 10", sprintf("%d of 10", sum(gs$true %in% sel)), all(gs$true %in% sel))
    } else {
      check(cohort, "05", paste0("no true gene selected without a response (", label, ")"),
            "0", sprintf("%d", sum(gs$true %in% sel)), !any(gs$true %in% sel))
    }
    check(cohort, "05", paste0("bleed genes not selected (", label, ")"), "0 selected",
          sprintf("%d selected; %d removed by the spillover filter",
                  sum(gs$bleed %in% sel), sum(gs$bleed %in% rem)),
          !any(gs$bleed %in% sel))
    nb <- sum(background %in% sel)
    check(cohort, "05", paste0("background genes selected (", label, ")"),
          "at most 1 of 100 (chance at BH 5%)", sprintf("%d (%s)", nb,
                                                        paste(background[background %in% sel], collapse = ", ")),
          nb <= 1)
  }
  g05(res, "in-situ filter")
  if (planted) g05("example/results_reference", "reference filter")
}

module_checks("planted", "example/results")
module_checks("null", "example/results_null")

# =============================================================================
# Adapter round trips (first three planted samples)
# =============================================================================
RT <- "example/roundtrip"
unlink(RT, recursive = TRUE); dir.create(RT)
sids <- c("S01", "S02", "S03")
orig <- lapply(setNames(sids, sids), function(s)
  list(cells = readRDS(file.path("example/input", s, "cells.rds")),
       counts = readRDS(file.path("example/input", s, "counts.rds"))))

same_input <- function(dir, coord_tol) {
  ok <- TRUE; worst <- 0
  for (s in sids) {
    a <- orig[[s]]; b <- list(cells = readRDS(file.path(dir, s, "cells.rds")),
                              counts = readRDS(file.path(dir, s, "counts.rds")))
    j <- match(a$cells$cell_id, sub("^S[0-9]+_(S[0-9]+_)", "\\1", b$cells$cell_id))
    if (anyNA(j)) return(list(ok = FALSE, worst = NA))
    bc <- b$cells[j, ]
    worst <- max(worst, abs(a$cells$x_um - bc$x_um), abs(a$cells$y_um - bc$y_um))
    ok <- ok && identical(a$cells$label, bc$label) &&
      isTRUE(all.equal(as.numeric(a$cells$n_count), as.numeric(bc$n_count))) &&
      max(abs(a$cells$x_um - bc$x_um), abs(a$cells$y_um - bc$y_um)) <= coord_tol &&
      identical(rownames(a$counts), rownames(b$counts)) &&
      max(abs(a$counts - b$counts[, j, drop = FALSE])) == 0
  }
  list(ok = ok, worst = worst)
}

# ---- Xenium Ranger folders: cells.csv.gz + cell_feature_matrix.h5 + labels ---
labs <- list()
for (s in sids) {
  d <- file.path(RT, "xr", s); dir.create(d, recursive = TRUE)
  cl <- orig[[s]]$cells
  z <- gzfile(file.path(d, "cells.csv.gz"), "w")
  utils::write.csv(data.frame(cell_id = cl$cell_id, x_centroid = cl$x_um,
                              y_centroid = cl$y_um, transcript_counts = cl$n_count,
                              control_probe_counts = 0L), z, row.names = FALSE)
  close(z)
  m <- orig[[s]]$counts
  ctl <- Matrix(rpois(2 * ncol(m), 0.05), nrow = 2, sparse = TRUE)
  mm <- methods::as(rbind(m, ctl), "CsparseMatrix")
  h <- hdf5r::H5File$new(file.path(d, "cell_feature_matrix.h5"), mode = "w")
  g <- h$create_group("matrix")
  g[["barcodes"]] <- cl$cell_id
  g[["data"]] <- as.integer(mm@x); g[["indices"]] <- as.integer(mm@i)
  g[["indptr"]] <- as.integer(mm@p); g[["shape"]] <- as.integer(dim(mm))
  fg <- g$create_group("features")
  fg[["name"]] <- c(rownames(m), "NegControlProbe_1", "UnassignedCodeword_1")
  fg[["id"]] <- c(rownames(m), "NegControlProbe_1", "UnassignedCodeword_1")
  fg[["feature_type"]] <- c(rep("Gene Expression", nrow(m)), "Negative Control Probe",
                            "Unassigned Codeword")
  h$close_all()
  labs[[s]] <- data.frame(sample_id = s, cell_id = cl$cell_id[!is.na(cl$label)],
                          label = cl$label[!is.na(cl$label)])
}
utils::write.csv(do.call(rbind, labs), file.path(RT, "xr_labels.csv"), row.names = FALSE)
writeLines(c('source("config.R", local = TRUE)',
             sprintf('INPUT_DIR <- "%s/xr_input"', RT),
             sprintf('XR_DIRS <- c(%s)', paste(sprintf('%s = "%s/xr/%s"', sids, RT, sids),
                                               collapse = ", ")),
             sprintf('XR_LABELS <- "%s/xr_labels.csv"', RT)),
           file.path(RT, "config_xr.R"))
run("00a_import_xenium_ranger.R", file.path(RT, "config_xr.R"))
r <- same_input(file.path(RT, "xr_input"), 1e-6)
check("adapter", "00a", "Xenium Ranger round trip reproduces the standard input",
      "identical cells, labels, n_count and gene counts",
      sprintf("max coordinate difference %.2g um", r$worst), r$ok)

# ---- Seurat: one object per sample, pixel coordinates, merged --------------
suppressPackageStartupMessages(library(SeuratObject))
PX <- 0.2125   # um per pixel
objs <- lapply(sids, function(s) {
  cl <- orig[[s]]$cells; m <- orig[[s]]$counts
  o <- suppressWarnings(CreateSeuratObject(counts = m, assay = "Xenium"))
  o$sample <- s; o$cell_label <- cl$label; o$total_tx <- cl$n_count
  cen <- CreateCentroids(data.frame(x = cl$x_um / PX, y = cl$y_um / PX, cell = cl$cell_id))
  o[[paste0("fov_", s)]] <- suppressWarnings(CreateFOV(cen, type = "centroids",
                                                       assay = "Xenium", key = "fov_"))
  o
})
mrg <- suppressWarnings(merge(objs[[1]], objs[-1], add.cell.ids = sids))
saveRDS(mrg, file.path(RT, "merged_seurat.rds"))
writeLines(c('source("config.R", local = TRUE)',
             sprintf('INPUT_DIR <- "%s/seurat_input"', RT),
             sprintf('SEURAT_RDS <- "%s/merged_seurat.rds"', RT),
             'SEURAT_LABEL_COL <- "cell_label"', 'SEURAT_SAMPLE_COL <- "sample"',
             'SEURAT_COUNT_COL <- "total_tx"', sprintf('COORD_SCALE <- %s', PX)),
           file.path(RT, "config_seurat.R"))
run("00b_import_seurat.R", file.path(RT, "config_seurat.R"))
r <- same_input(file.path(RT, "seurat_input"), 1e-6)
check("adapter", "00b", "merged Seurat object (pixels x COORD_SCALE) reproduces the standard input",
      "identical cells, labels, n_count and gene counts; coordinates within 1e-6 um",
      sprintf("max coordinate difference %.2g um", r$worst), r$ok)

# =============================================================================
CK <- do.call(rbind, CHECKS)
wb <- createWorkbook(); addWorksheet(wb, "checks"); writeData(wb, "checks", CK)
saveWorkbook(wb, "example/validation_report.xlsx", overwrite = TRUE)
message(sprintf("\n%d of %d checks passed. Report: example/validation_report.xlsx",
                sum(CK$pass), nrow(CK)))
if (!all(CK$pass)) quit(status = 1)
