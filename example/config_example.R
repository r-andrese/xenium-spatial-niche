# =============================================================================
# example/config_example.R -- configuration for the simulated planted cohort
#
# Starts from the defaults in config.R and sets only what this dataset needs.
# Run from the repository root, e.g.
#   Rscript 02_colocalisation.R example/config_example.R
# ASCII only.
# =============================================================================
source("config.R", local = TRUE)

INPUT_DIR   <- "example/input"
RESULTS_DIR <- "example/results"

# What is compared
COMP_A   <- "A"
COMP_B   <- "B"
STATES_A <- c("A1", "A2")
STATES_B <- c("B1", "B2", "B3")
PRIMARY_PAIRS <- list(c("A1", "B1"))

# Module 04
EXPOSURE_GENE <- "GENEX"

# Module 05
INDEX_STATES     <- c("A1", "A2")
NEIGHBOUR_STATES <- "B1"
SPILL_REFERENCE  <- NULL     # in-situ filter; validate_example.R also runs
                             # it with example/truth/spill_reference.csv
