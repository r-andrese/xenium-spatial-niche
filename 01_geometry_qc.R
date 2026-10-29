#!/usr/bin/env Rscript
# =============================================================================
# 01_geometry_qc.R -- tissue geometry and contiguity
#
# Question
#   Is each section one contiguous piece of tissue or several fragments, 
#   and at what spatial scale are its cells packed? Every later module
#   permutes labels within a tissue component.
#
# Method
#   Cells are binned on a BIN_UM grid. A bin is tissue if it holds at least
#   MIN_CELLS_BIN cells. Occupied bins are joined by rook (CONNECTIVITY = 4)
#   or queen (8) adjacency. Components smaller than MIN_REGION_BINS bins are
#   debris and get component 0, as do cells in sparse bins; such cells are
#   excluded from every later module. Components are found on the grid rather
#   than on the cells so that a few stray segmented objects cannot bridge two
#   separate pieces of tissue.
#
# Also reported, per sample
#   bounding box; tissue area (occupied bins); cell density; median
#   nearest-neighbour distance (the scale reference for RADIUS_UM); fraction
#   of cells unlabelled; median n_count overall and per label (the reference
#   for MIN_MEDIAN_DEPTH). The number of components is recomputed with
#   MIN_CELLS_BIN in {1, 3, 5, 10} so its sensitivity to that threshold is
#   visible.
#
# Reads   INPUT_DIR/samples.csv, label_map.csv, <sample>/cells.rds
# Writes  RESULTS_DIR/01_geometry_qc/
#           01_geometry_qc.xlsx  per_sample, per_component, depth_by_label,
#                                sensitivity
#           01_tissue_maps.pdf   occupied bins coloured by component
#           components.rds       sample_id, cell_id, component (read by 02-05)
#
# Usage   Rscript 01_geometry_qc.R [config.R]
# =============================================================================

source("helpers.R")
load_config()
OUT <- module_dir(RESULTS_DIR, "01_geometry_qc")
SENS_MIN_CELLS <- c(1, 3, 5, 10)

samples   <- read_samples(INPUT_DIR, SAMPLES)
label_map <- read_label_map(INPUT_DIR)
msg("[1] %d samples", nrow(samples))

per_sample <- list(); per_comp <- list(); depth_lab <- list()
sens <- list(); comp_rows <- list(); maps <- list()

for (sid in samples$sample_id) {
  s  <- load_sample(sid, INPUT_DIR, label_map)
  cl <- s$cells
  ok <- is.finite(cl$x_um) & is.finite(cl$y_um)
  comp <- rep(0L, nrow(cl))
  comp[ok] <- tissue_components(cl$x_um[ok], cl$y_um[ok], BIN_UM, MIN_CELLS_BIN,
                                CONNECTIVITY, MIN_REGION_BINS)
  comp_rows[[sid]] <- data.frame(sample_id = sid, cell_id = cl$cell_id,
                                 component = comp, stringsAsFactors = FALSE)

  tis <- comp > 0
  bkey <- paste(floor(cl$x_um[tis] / BIN_UM), floor(cl$y_um[tis] / BIN_UM))
  area_mm2 <- length(unique(bkey)) * (BIN_UM / 1000)^2

  per_sample[[sid]] <- data.frame(
    sample_id = sid, n_cells = nrow(cl),
    n_in_tissue = sum(tis), n_debris_or_sparse = sum(ok & !tis),
    n_no_coordinates = sum(!ok),
    n_components = length(unique(comp[tis])),
    x_min_um = min(cl$x_um[ok]), x_max_um = max(cl$x_um[ok]),
    y_min_um = min(cl$y_um[ok]), y_max_um = max(cl$y_um[ok]),
    width_mm = diff(range(cl$x_um[ok])) / 1000,
    height_mm = diff(range(cl$y_um[ok])) / 1000,
    tissue_area_mm2 = area_mm2,
    cells_per_mm2 = if (area_mm2 > 0) sum(tis) / area_mm2 else NA_real_,
    median_nn_distance_um = median_nn_distance(cl$x_um[tis], cl$y_um[tis]),
    pct_unlabelled = 100 * mean(is.na(cl$label)),
    median_n_count_all = stats::median(cl$n_count, na.rm = TRUE),
    median_n_count_labelled = stats::median(cl$n_count[!is.na(cl$label)], na.rm = TRUE),
    stringsAsFactors = FALSE)

  for (k in sort(unique(comp[tis]))) {
    sel <- comp == k
    per_comp[[length(per_comp) + 1]] <- data.frame(
      sample_id = sid, component = k, n_cells = sum(sel),
      n_bins = length(unique(paste(floor(cl$x_um[sel] / BIN_UM),
                                   floor(cl$y_um[sel] / BIN_UM)))),
      x_centre_um = mean(cl$x_um[sel]), y_centre_um = mean(cl$y_um[sel]),
      stringsAsFactors = FALSE)
  }

  lab <- ifelse(is.na(cl$label), "(unlabelled)", cl$label)
  dl <- stats::aggregate(cl$n_count[tis], by = list(label = lab[tis]),
                         FUN = function(v) c(n = length(v), med = stats::median(v)))
  depth_lab[[sid]] <- data.frame(sample_id = sid, label = dl$label,
                                 compartment = label_map$compartment[match(dl$label, label_map$label)],
                                 n_cells = dl$x[, "n"], median_n_count = dl$x[, "med"],
                                 stringsAsFactors = FALSE)

  for (m in SENS_MIN_CELLS) {
    cs <- tissue_components(cl$x_um[ok], cl$y_um[ok], BIN_UM, m, CONNECTIVITY,
                            MIN_REGION_BINS)
    sens[[length(sens) + 1]] <- data.frame(
      sample_id = sid, min_cells_bin = m, n_components = length(unique(cs[cs > 0])),
      n_cells_excluded = sum(cs == 0), stringsAsFactors = FALSE)
  }

  # occupancy map for the figure: every bin with at least one cell
  bx <- floor(cl$x_um[ok] / BIN_UM); by <- floor(cl$y_um[ok] / BIN_UM)
  key <- paste(bx, by)
  first <- !duplicated(key)
  bc <- tapply(comp[ok], key, max)
  maps[[sid]] <- data.frame(sample_id = sid, x = (bx[first] + 0.5) * BIN_UM / 1000,
                            y = (by[first] + 0.5) * BIN_UM / 1000,
                            component = as.integer(bc[key[first]]))

  msg("    %s  cells=%d  components=%d  excluded=%d  median NN=%.1f um",
      sid, nrow(cl), per_sample[[sid]]$n_components,
      per_sample[[sid]]$n_debris_or_sparse, per_sample[[sid]]$median_nn_distance_um)
}

PS <- do.call(rbind, per_sample)
saveRDS(do.call(rbind, comp_rows), file.path(OUT, "components.rds"))
msg("    written: %s", file.path(OUT, "components.rds"))

write_xlsx(list(per_sample = PS, per_component = do.call(rbind, per_comp),
                depth_by_label = do.call(rbind, depth_lab),
                sensitivity = do.call(rbind, sens)),
           file.path(OUT, "01_geometry_qc.xlsx"))

# -----------------------------------------------------------------------------
# Figure: one panel per sample, bins coloured by component, debris in grey
# -----------------------------------------------------------------------------
MP <- do.call(rbind, maps)
MP$component_f <- ifelse(MP$component == 0, "excluded", paste0("c", MP$component))
lv <- c("excluded", setdiff(sort(unique(MP$component_f)), "excluded"))
MP$component_f <- factor(MP$component_f, levels = lv)
pal <- c(excluded = "grey80",
         stats::setNames(rep(c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
                               "#DDCC77", "#CC6677", "#882255", "#AA4499"),
                             length.out = length(lv) - 1), lv[-1]))

# One page per sample, true aspect ratio (sections differ in size).
f <- file.path(OUT, "01_tissue_maps.pdf")
open_pdf(f, width = 7, height = 6)
for (sid in samples$sample_id) {
  d <- MP[MP$sample_id == sid, ]
  print(ggplot(d, aes(x, y, fill = component_f)) +
          geom_tile(width = BIN_UM / 1000, height = BIN_UM / 1000) +
          coord_fixed() +
          scale_fill_manual(values = pal, name = "component", drop = TRUE) +
          labs(x = "x (mm)", y = "y (mm)", title = sid) +
          theme_niche(BASE_SIZE))
}
close_pdf(f)

write_session_info(OUT)
msg("[2] Done.")
