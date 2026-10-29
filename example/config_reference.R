# =============================================================================
# example/config_reference.R -- module 05 on the planted cohort with a
# reference table for the spillover filter instead of the in-situ filter.
# validate_example.R copies the 01 component assignments into RESULTS_DIR
# before running 05 with this file. ASCII only.
# =============================================================================
source("example/config_example.R", local = TRUE)

RESULTS_DIR     <- "example/results_reference"
SPILL_REFERENCE <- "example/truth/spill_reference.csv"
