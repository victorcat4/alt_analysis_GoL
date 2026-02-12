# =============================================================================
# GULF OF LION BIRD MIGRATION ANALYSIS
# =============================================================================
#
# Main analysis script for GPS tracking data processing and visualization.
#
# Pipeline:
#   1. (Optional) Process all data sources (Movebank + non-Movebank)
#   2. Summary statistics
#   3. Load processed data & add ecological groups
#   4. Visualization
#   5a. (Optional) Bayesian model fitting
#   5b. (Optional) Generate plots from model results
#
# =============================================================================


# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

# Paths
movebank_dir     <- "D:/migralion_mb_data/individuals"
non_movebank_dir <- "D:/migralion_none_mb_data"
output_dir       <- "output"
plots_dir        <- "output/plots"

# Processing parameters
DOWNSAMPLE_MIN   <- 10   # Minutes between GPS fixes
MAX_GAP_DAYS     <- 5    # Days - larger gaps split trajectories
MIN_TRAJ_POINTS  <- 5    # Minimum points to keep a trajectory
DBSCAN_EPS       <- 0.1  # DBSCAN epsilon (degrees)
DBSCAN_MINPTS    <- 10   # DBSCAN minimum cluster points

# Parallel processing
N_WORKERS        <- 10   # Number of workers (NULL = auto-detect: cores - 1)
BATCH_SIZE       <- 10   # Individuals per batch for memory management

# MCMC
n_chains <- 5
n_it <- 500000
p_burnin <- 0.9
thin <- 20

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

# Run flags (set FALSE to skip already-completed steps)
RUN_PROCESSING   <- FALSE  # Phase 1: data processing (slow, only needed once)
RUN_MODELING     <- TRUE  # Phase 5a: Bayesian model fitting (computationally intensive)
RUN_MODEL_PLOTS  <- TRUE   # Phase 5b: generate plots from saved model results


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
#
# Strategy: Each data source processes individuals and saves to temp files,
# then merges all temp files into the final gps_migration.csv.
# This avoids file access conflicts during parallel processing.
#
# -----------------------------------------------------------------------------

output_csv <- file.path(output_dir, "gps_migration.csv")

if (RUN_PROCESSING) {

message("\n", strrep("=", 70))
message("PHASE 1: Data Processing")
message(strrep("=", 70), "\n")

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
hb_dir <- file.path(non_movebank_dir, "BondrÃƒÂ©e_apivore_patrick_byholm")
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

} else {
  message("\n[NOTE] Data processing skipped (set RUN_PROCESSING = TRUE to enable)")
  message(sprintf("  Using existing output: %s", output_csv))
}


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
# PHASE 3: LOAD DATA & ADD ECOLOGICAL GROUPS
# -----------------------------------------------------------------------------

message("\n", strrep("=", 70))
message("PHASE 3: Load Processed Data")
message(strrep("=", 70), "\n")

# Load processed GPS data
message("  Loading processed GPS data...")
gps <- read_csv(
  output_csv,
  col_types = cols(.default = col_character()),
  progress = FALSE,
  na = c("", "NA", "NULL", "na", "null")
)

message(sprintf("  Raw records loaded: %s", format(nrow(gps), big.mark = ",")))

# Convert columns to appropriate types
gps <- gps %>%
  mutate(
    location_long = as.numeric(location_long),
    location_lat = as.numeric(location_lat),
    height_above_msl = as.numeric(height_above_msl),
    distance_to_shore = as.numeric(distance_to_shore),
    in_GOL = as.logical(in_GOL),
    db_cluster = as.integer(db_cluster),
    gps_hdop = as.numeric(gps_hdop),
    gps_vdop = as.numeric(gps_vdop),
    gps_satellite_count = as.numeric(gps_satellite_count)
  )

# Validate coordinates (should be clean after pipeline fix)
n_na_coords <- sum(is.na(gps$location_long) | is.na(gps$location_lat))
if (n_na_coords > 0) {
  warning(sprintf("Found %d rows with NA coordinates - filtering out (investigate pipeline bug!)", n_na_coords))
  gps <- gps %>% filter(!is.na(location_long), !is.na(location_lat))
}

# Remove duplicated records (same location + timestamp, keep first occurrence)
n_before_dedup <- nrow(gps)
gps <- gps %>%
  distinct(location_long, location_lat, timestamp, .keep_all = TRUE)
n_dupes <- n_before_dedup - nrow(gps)
if (n_dupes > 0) {
  message(sprintf("  Removed %s duplicate records (same location & timestamp)",
                  format(n_dupes, big.mark = ",")))
}

# Add ecological groups
message("  Adding ecological groups...")
ec_file <- "data/supp/EC.xlsx"
if (file.exists(ec_file)) {
  
  # Load ecological group mapping
  ec_mapping <- read_excel(ec_file) %>%
    select(species = specie_scientific, ecological_group = ecological_group)
  
  # Join
  gps <- gps %>%
    left_join(ec_mapping, by = c("individual_taxon_canonical_name" = "species"))
  
  n_with_group <- sum(!is.na(gps$ecological_group))
  n_species_matched <- n_distinct(gps$individual_taxon_canonical_name[!is.na(gps$ecological_group)])
  message(sprintf("  Matched %d species (%s records with ecological group)",
                  n_species_matched, format(n_with_group, big.mark = ",")))
  
} else {
  warning("Ecological groups file not found: ", ec_file)
  gps$ecological_group <- NA_character_
}

# Summary statistics
message("\n  Data Summary:")
message(sprintf("    Total records: %s", format(nrow(gps), big.mark = ",")))
message(sprintf("    Unique individuals: %d", n_distinct(gps$individual_id)))
message(sprintf("    Unique species: %d", n_distinct(gps$individual_taxon_canonical_name)))
message(sprintf("    Date range: %s to %s", min(gps$timestamp), max(gps$timestamp)))

# Records by season and state
message("\n  Records by Season & State:")
gps %>%
  count(season, state) %>%
  pivot_wider(names_from = state, values_from = n, values_fill = 0) %>%
  print()

# Records by ecological group
if (any(!is.na(gps$ecological_group))) {
  message("\n  Records by Ecological Group:")
  gps %>%
    filter(!is.na(ecological_group)) %>%
    count(ecological_group, name = "n_records") %>%
    arrange(desc(n_records)) %>%
    print()
}


# -----------------------------------------------------------------------------
# PHASE 4: VISUALIZATION
# -----------------------------------------------------------------------------

message("\n", strrep("=", 70))
message("PHASE 4: Visualization")
message(strrep("=", 70), "\n")

# Generate plots
message("  Generating visualizations...")
create_all_visualizations(gps, spatial$coastline, plots_dir)
message(sprintf("  Plots saved to: %s", plots_dir))


# -----------------------------------------------------------------------------
# PHASE 5: BAYESIAN MODELING (Optional)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# PHASE 5a: BAYESIAN MODEL FITTING (Optional)
# -----------------------------------------------------------------------------

if (RUN_MODELING) {
  message("\n", strrep("=", 70))
  message("PHASE 5a: Bayesian Model Fitting")
  message(strrep("=", 70), "\n")
  
  source("R/06_bayesian_modeling.R")
  
  # Pre-breeding migration
  results_pre <- run_flight_height_model(
    gps_data = gps, 
    output_dir = output_dir, 
    season = "pre-breeding",
    n_iter = n_it, 
    burnin_prop = p_burnin, 
    thin = thin, 
    chains = n_chains
  )
  
  # Post-breeding migration
  results_post <- run_flight_height_model(
    gps_data = gps, 
    output_dir = output_dir, 
    season = "post-breeding",
    n_iter = n_it, 
    burnin_prop = p_burnin, 
    thin = thin, 
    chains = n_chains
  )
  
} else {
  message("\n[NOTE] Model fitting skipped (set RUN_MODELING = TRUE to enable)")
}


# -----------------------------------------------------------------------------
# PHASE 5b: MODEL PLOTS (from saved or freshly fitted models)
# -----------------------------------------------------------------------------

if (RUN_MODEL_PLOTS) {
  message("\n", strrep("=", 70))
  message("PHASE 5b: Model Plots")
  message(strrep("=", 70), "\n")
  
  source("R/07_plot_modeling.R")
  
  if (RUN_MODELING) {
    # Use freshly fitted results
    if (!is.null(results_pre)) {
      generate_model_plots(results_pre, plots_dir, season_filter = "pre-breeding")
    }
    if (!is.null(results_post)) {
      generate_model_plots(results_post, plots_dir, season_filter = "post-breeding")
    }
  } else {
    # Load from saved model files
    model_pre_path  <- file.path(output_dir, "model", "model_pre-breeding.rds")
    model_post_path <- file.path(output_dir, "model", "model_post-breeding.rds")
    
    if (file.exists(model_pre_path)) {
      message("  Loading saved pre-breeding model...")
      plot_from_saved_model(model_pre_path, gps, plots_dir, season_filter = "pre-breeding")
    } else {
      message("  Pre-breeding model not found: ", model_pre_path)
    }
    
    if (file.exists(model_post_path)) {
      message("  Loading saved post-breeding model...")
      plot_from_saved_model(model_post_path, gps, plots_dir, season_filter = "post-breeding")
    } else {
      message("  Post-breeding model not found: ", model_post_path)
    }
  }
  
} else {
  message("\n[NOTE] Model plots skipped (set RUN_MODEL_PLOTS = TRUE to enable)")
}


# -----------------------------------------------------------------------------
# COMPLETE
# -----------------------------------------------------------------------------

message("\n", strrep("=", 70))
message("ANALYSIS COMPLETE")
message(strrep("=", 70))
message(sprintf("\nOutputs: %s, %s/", output_csv, plots_dir))