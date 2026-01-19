# =============================================================================
# GULF OF LION BIRD MIGRATION ANALYSIS
# =============================================================================
#
# Main analysis script for GPS tracking data processing and visualization.
#
# Pipeline:
#   1. Load and process all data sources (Movebank + non-Movebank)
#   2. Generate summary statistics
#   3. Create visualizations
#   4. (Optional) Run Bayesian modeling
#
# =============================================================================


# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

# Paths
movebank_dir     <- "D:/migralion_mb_data/individuals"
non_movebank_dir <- "data/none_MB"
output_dir       <- "output"
plots_dir        <- "output/plots"

# Processing parameters
DOWNSAMPLE_MIN   <- 10   # Minutes between GPS fixes
MAX_GAP_DAYS     <- 5    # Days - larger gaps split trajectories
MIN_TRAJ_POINTS  <- 5    # Minimum points to keep a trajectory
DBSCAN_EPS       <- 0.1  # DBSCAN epsilon (degrees)
DBSCAN_MINPTS    <- 10   # DBSCAN minimum cluster points

# Parallel processing
N_WORKERS        <- 3 # Number of workers (NULL = auto-detect: cores - 1)
BATCH_SIZE       <- 10  # Individuals per batch for memory management

# Species to exclude
EXCLUDED_SPECIES <- c(
  "Puffinus yelkouan",
  "Calonectris diomedea",
  "Thalasseus sandvicensis",
  "Larus michahellis",
  "Aquila clanga",
  "Aquila chrysaetos",
  "Ciconia ciconia"
)

# Run Bayesian modeling? (computationally intensive)
RUN_MODELING <- TRUE


# -----------------------------------------------------------------------------
# SETUP
# -----------------------------------------------------------------------------

message("\n", strrep("=", 70))
message("GULF OF LION BIRD MIGRATION ANALYSIS")
message(strrep("=", 70), "\n")

# Load packages and functions
source("R/00_setup.R")
source("R/01_formatting_functions.R")
source("R/02_spatial_functions.R")
source("R/03_data_readers.R")
source("R/04_visualization_functions.R")
source("R/05_processing_pipeline.R")

# Create output directories
dir.create(output_dir, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Load spatial reference data
message("[SETUP] Loading spatial data...")
spatial <- load_spatial_data()


# -----------------------------------------------------------------------------
# PHASE 1: DATA PROCESSING
# -----------------------------------------------------------------------------

message("\n", strrep("=", 70))
message("PHASE 1: Data Processing")
message(strrep("=", 70), "\n")

output_csv <- file.path(output_dir, "gps_migration.csv")
if (file.exists(output_csv)) file.remove(output_csv)

# Common processing arguments
proc_args <- list(
  output_csv        = output_csv,
  excluded_species  = EXCLUDED_SPECIES,
  area_res          = spatial$restricted_area,
  area_ext          = spatial$extended_area,
  coastline         = spatial$coastline,
  downsample_interval = DOWNSAMPLE_MIN,
  max_gap_days      = MAX_GAP_DAYS,
  min_points        = MIN_TRAJ_POINTS,
  dbscan_eps        = DBSCAN_EPS,
  dbscan_minpts     = DBSCAN_MINPTS,
  # Parallel processing (only used for Movebank)
  n_workers         = N_WORKERS,
  batch_size        = BATCH_SIZE
)

stats <- list(total = 0, kept = 0)


# 1.1 Movebank data
message("[1/5] Movebank data...")
if (dir.exists(movebank_dir)) {
  result <- do.call(process_movebank_individuals, c(list(base_dir = movebank_dir), proc_args))
  stats$total <- stats$total + result$total
  stats$kept  <- stats$kept + result$kept
} else {
  message("  Directory not found, skipping.")
}


# 1.2 Turtle Dove & Stone Curlew
message("\n[2/5] Turtle Dove & Stone Curlew...")
td_dir <- file.path(non_movebank_dir, "Tourterelle&Oedicnemes")
if (dir.exists(td_dir)) {
  result <- do.call(process_non_movebank_source, 
                    c(list(reader_function = read_turtle_dove_stone_curlew, data_dir = td_dir), proc_args))
  stats$total <- stats$total + result$total
  stats$kept  <- stats$kept + result$kept
} else {
  message("  Directory not found, skipping.")
}


# 1.3 Honey Buzzard (Patrick Byholm)
message("\n[3/5] Honey Buzzard (Byholm)...")
hb_dir <- file.path(non_movebank_dir, "Bondrée_apivore_patrick_byholm")
if (dir.exists(hb_dir)) {
  result <- do.call(process_non_movebank_source,
                    c(list(reader_function = read_honey_buzzard_byholm, data_dir = hb_dir), proc_args))
  stats$total <- stats$total + result$total
  stats$kept  <- stats$kept + result$kept
} else {
  message("  Directory not found, skipping.")
}


# 1.4 BIRDMAP
message("\n[4/5] BIRDMAP...")
bm_dir <- file.path(non_movebank_dir, "BIRDMAP")
if (dir.exists(bm_dir)) {
  result <- do.call(process_non_movebank_source,
                    c(list(reader_function = read_birdmap, data_dir = bm_dir), proc_args))
  stats$total <- stats$total + result$total
  stats$kept  <- stats$kept + result$kept
} else {
  message("  Directory not found, skipping.")
}


# 1.5 Montagu's Harrier
message("\n[5/5] Montagu's Harrier...")
mh_dir <- file.path(non_movebank_dir, "schlaich_et_al.2016")
if (dir.exists(mh_dir)) {
  result <- do.call(process_non_movebank_source,
                    c(list(reader_function = read_montagu_harrier, data_dir = mh_dir), proc_args))
  stats$total <- stats$total + result$total
  stats$kept  <- stats$kept + result$kept
} else {
  message("  Directory not found, skipping.")
}


message("\n", strrep("-", 50))
message(sprintf("PHASE 1 COMPLETE: %d/%d individuals retained", stats$kept, stats$total))
message(strrep("-", 50))


# -----------------------------------------------------------------------------
# PHASE 2: SUMMARY STATISTICS
# -----------------------------------------------------------------------------

message("\n", strrep("=", 70))
message("PHASE 2: Summary Statistics")
message(strrep("=", 70), "\n")

summary_stats <- get_csv_summary(output_csv)

message(sprintf("  Total GPS records: %s", format(summary_stats$n_records, big.mark = ",")))
message(sprintf("  Unique individuals: %d", summary_stats$n_individuals))
message(sprintf("  Unique species: %d\n", summary_stats$n_species))

print(summary_stats$species_counts, n = 25)


# -----------------------------------------------------------------------------
# PHASE 3: VISUALIZATION
# -----------------------------------------------------------------------------

message("\n", strrep("=", 70))
message("PHASE 3: Visualization")
message(strrep("=", 70), "\n")

# Load processed data
gps <- load_processed_data(output_csv)

# Add ecological groups if available
ec_file <- "data/supp/EC.xlsx"
if (file.exists(ec_file)) {
  ec_mapping <- load_ecological_groups(ec_file)
  gps <- add_ecological_groups(gps, ec_mapping)
  message(sprintf("  Added ecological groups for %d species\n", n_distinct(ec_mapping$species)))
}

# Generate plots
message("  Generating visualizations...")
create_all_visualizations(gps, spatial$coastline, plots_dir)
message(sprintf("  Plots saved to: %s", plots_dir))


# -----------------------------------------------------------------------------
# PHASE 4: BAYESIAN MODELING (Optional)
# -----------------------------------------------------------------------------

if (RUN_MODELING) {
  message("\n", strrep("=", 70))
  message("PHASE 4: Bayesian Modeling")
  message(strrep("=", 70), "\n")
  
  source("R/06_bayesian_modeling.R")
  results <- run_flight_height_model(gps_data = gps, output_dir = output_dir, plots_dir = plots_dir, n_iter = 1000000, burnin_prop = 0.9, thin = 20, chains = 10)
} else {
  message("\n[NOTE] Bayesian modeling skipped (set RUN_MODELING = TRUE to enable)")
}


# -----------------------------------------------------------------------------
# COMPLETE
# -----------------------------------------------------------------------------

message("\n", strrep("=", 70))
message("ANALYSIS COMPLETE")
message(strrep("=", 70))
message(sprintf("\nOutputs: %s, %s/", output_csv, plots_dir))