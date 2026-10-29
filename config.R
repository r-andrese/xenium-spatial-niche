# =============================================================================
# config.R -- the only file to edit
#
# Every script reads this file (or the file named as its first argument):
#   Rscript 02_colocalisation.R               # uses config.R
#   Rscript 02_colocalisation.R my_config.R   # uses my_config.R
#
# DEFAULT THRESHOLDS. Every threshold below is the value used in the analysis
# this toolkit was built from: one breast cancer cohort, Xenium 5,000-gene
# panel, median ~180 transcripts per cell, ~3,000-5,000 cells per mm2. Check 
# each against your own data (01_geometry_qc.R reports the density and depth 
# figures needed).
# =============================================================================

# -----------------------------------------------------------------------------
# Folders
# -----------------------------------------------------------------------------
INPUT_DIR   <- "input"      # standard input, written by 00a or 00b
RESULTS_DIR <- "results"    # one sub-folder per module
SAMPLES     <- NULL         # NULL = every sample in input/samples.csv

# -----------------------------------------------------------------------------
# 1. What is compared (modules 02-04)
# -----------------------------------------------------------------------------
# input/label_map.csv assigns every label to a broad compartment, e.g.
#   label,compartment
#   fibro-prematrix,fibroblast
#   CAFs,fibroblast
#   mac-M2-CXCL,myeloid
#   mac-angio,myeloid
#   cancer-basal,tumour
#
# COMP_A: reference compartment, the cells whose neighbourhood is examined.
# COMP_B: target compartment, the cells counted around them.
# Labels are shuffled only within their own compartment (and tissue
# component), so the tests ask which B state sits nearest which A state,
# given where each compartment already is.
COMP_A <- "fibroblast"
COMP_B <- "myeloid"

# Fine labels to test. NULL = every label of the compartment. Each
# A state x B state combination is one test.
STATES_A <- NULL
STATES_B <- NULL

# A x B pairs for the primary analysis. They are corrected for
# multiple testing as their own family; every other pair is exploratory.
# list() makes the whole run exploratory.
PRIMARY_PAIRS <- list()
# e.g. PRIMARY_PAIRS <- list(c("fibro-prematrix", "mac-M2-CXCL"))

# -----------------------------------------------------------------------------
# 2. Tissue components (all modules)
# -----------------------------------------------------------------------------
BIN_UM          <- 100   # occupancy grid bin (um)
MIN_CELLS_BIN   <- 5     # cells for a bin to count as tissue
CONNECTIVITY    <- 4     # 4 = rook (splits diagonally touching pieces,
                         # conservative); 8 = queen (joins them, permissive)
MIN_REGION_BINS <- 10    # components smaller than this are debris
                         # (10 bins of 100 um = 0.1 mm2)

# -----------------------------------------------------------------------------
# 3. Shared thresholds
# -----------------------------------------------------------------------------
MIN_STATE_CELLS <- 50    # cells of a state in a sample for it to contribute
MIN_SAMPLES     <- 5     # samples contributing for a cohort test (02-04)
RADIUS_UM       <- 50    # neighbourhood radius (02, 04, 05)
SEED            <- 20260819
BASE_SIZE       <- 14    # figure base font size

# -----------------------------------------------------------------------------
# Module 02 - co-localisation
# -----------------------------------------------------------------------------
N_PERM_PRIMARY <- 1000   # within-compartment null
N_PERM_GLOBAL  <- 200    # global null (all labels shuffled within component)

# -----------------------------------------------------------------------------
# Module 03 - distance decay
# -----------------------------------------------------------------------------
BIN_WIDTH_UM   <- 25     # distance bin width
MAX_D_UM       <- 200    # cells further than this are censored (also 04)
N_PERM_DECAY   <- 200    # permutations for the expected curve and envelope
MIN_BIN_B      <- 20     # compartment-B cells in a bin for a proportion
MIN_BINS_FIT   <- 3      # bins needed for a per-sample decay slope
STAR_DIRECTION <- "positive"   # "positive" (attraction, e.g. a niche),
                               # "negative" (avoidance, e.g. exclusion) or
                               # "both" (exploratory). Set before the run.

# -----------------------------------------------------------------------------
# Module 04 - exposure contrast
# -----------------------------------------------------------------------------
EXPOSURE_GENE    <- NULL # gene whose detection (count > 0) splits A cells
MIN_MEDIAN_DEPTH <- 100  # sample excluded if the median n_count of labelled
                         # cells is below this (04, 05)
MIN_ARM          <- 20   # A cells per arm (gene+ and gene-) per sample
N_DEPTH_STRATA   <- 10   # within-sample depth quantile strata

# -----------------------------------------------------------------------------
# Module 05 - neighbourhood-associated genes
# -----------------------------------------------------------------------------
INDEX_STATES        <- NULL   # labels pooled as index cells (one compartment)
NEIGHBOUR_STATES    <- NULL   # labels pooled as neighbour cells
MIN_INDEX_CELLS     <- 50
MIN_NEIGHBOUR_CELLS <- 50
MIN_DETECT_FRAC     <- 0.02   # index-cell detection rate for a gene to be
                              # tested in a sample
MIN_SAMPLES_GENE    <- 10     # samples in which a gene must be tested
FDR_MAX             <- 0.05
MIN_FRAC_POS        <- 2 / 3  # fraction of samples with a positive correlation
N_TOP               <- 25     # selected genes kept, by median correlation
SPILL_REFERENCE     <- NULL   # optional CSV: gene, then one column per
                              # compartment of mean expression (a matched
                              # scRNA-seq reference is preferred). NULL = the
                              # in-situ filter.
SPILL_RATIO_MIN     <- 2

# -----------------------------------------------------------------------------
# Re-running
# -----------------------------------------------------------------------------
FORCE_REBUILD <- FALSE   # TRUE recomputes even if a valid checkpoint exists
REPLOT_ONLY   <- FALSE   # TRUE redraws figures and tables from the checkpoint

# -----------------------------------------------------------------------------
# Import: 00a_import_xenium_ranger.R
# -----------------------------------------------------------------------------
# One Xenium Ranger / Xenium Onboard Analysis output folder per sample,
# named by sample_id.
XR_DIRS <- c()
# e.g. XR_DIRS <- c(S01 = "/data/run1/output-XETG..._S01", S02 = "...")

# Cell labels: a CSV with cell_id and label (one file for all samples with a
# sample_id column, or a named vector of one file per sample).
XR_LABELS <- NULL
XR_IMPORT_COUNTS <- TRUE  # FALSE skips the counts matrix (modules 04, 05 need it)

# -----------------------------------------------------------------------------
# Import: 00b_import_seurat.R
# -----------------------------------------------------------------------------
SEURAT_RDS        <- NULL   # path to a saved Seurat v5 object
SEURAT_LABEL_COL  <- NULL   # meta.data column with the cell labels
SEURAT_SAMPLE_COL <- NULL   # meta.data column with the sample id
SEURAT_COUNT_COL  <- NULL   # meta.data column with total counts; NULL =
                            # nCount_<assay>
SEURAT_ASSAY      <- NULL   # NULL = the default assay
SEURAT_IMPORT_COUNTS <- TRUE
COORD_SCALE       <- NULL   # REQUIRED: micrometres per coordinate unit.
                            # 1 for objects built with LoadXenium(). See the
                            # README for how to check it.
