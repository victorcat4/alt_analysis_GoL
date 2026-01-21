# =============================================================================
# Processing Pipeline Functions
# =============================================================================
#
# High-level functions that orchestrate the data processing workflow.
#
# STRATEGY: Save processed individuals to temporary files, then merge at the end.
# This avoids issues with concurrent file writes during parallel processing.
#
# =============================================================================


#' Apply Full Processing Pipeline to Individual Data
#'
#' Pipeline steps:
#' 1. Clean rows (remove invalid coordinates)
#' 2. Add equipment location
#' 3. Filter to study area
#' 4. Downsample
#' 5. Assign trajectories
#' 6. Remove short trajectories
#' 7. Add season
#' 8. Check restricted area
#' 9. Classify behavior
#' 10. Calculate distance to shore
#' 11. VALIDATE (final check for NAs)
#'
#' @param df Dataframe with GPS data for one individual
#' @param area_res Restricted study area polygon
#' @param area_ext Extended study area polygon
#' @param coastline Coastline polygon
#' @param downsample_interval Minutes between fixes
#' @param max_gap_days Days to split trajectories
#' @param min_points Minimum points per trajectory
#' @param dbscan_eps DBSCAN epsilon
#' @param dbscan_minpts DBSCAN minimum points
#' @return Processed dataframe (character columns for CSV storage)
process_individual <- function(df,
                                area_res,
                                area_ext,
                                coastline,
                                downsample_interval = 10,
                                max_gap_days = 5,
                                min_points = 5,
                                dbscan_eps = 0.1,
                                dbscan_minpts = 10) {
  
  # Step 1-2: Clean and add equipment location
  df <- df %>%
    clean_row() %>%
    first_loc()
  
  if (nrow(df) == 0) return(df)
  
  # Step 3: Spatial filter - keep only if crosses Gulf of Lion
  df <- select_pts(df, area_res = area_res, area_ext = area_ext)
  
  if (nrow(df) == 0) return(df)
  
  # Step 4: Downsample
  df <- downsample_tracks(df, time_interval = downsample_interval)
  
  if (nrow(df) == 0) return(df)
  
  # Step 5-6: Trajectories
  df <- df %>%
    assign_trajectories(max_gap_days = max_gap_days) %>%
    remove_short_trajectories(min_points = min_points)
  
  if (nrow(df) == 0) return(df)
  
  # Step 7-10: Classifications
  df <- df %>%
    add_season() %>%
    check_restricted_area(area_res = area_res) %>%
    classify_staging_migration(eps = dbscan_eps, minPts = dbscan_minpts) %>%
    calculate_distance_to_shore(coastline = coastline)
  
  # Step 11: FINAL VALIDATION - Filter any remaining NAs
  # This should never remove rows if clean_row() worked correctly,
  # but serves as a safety net
  n_before <- nrow(df)
  df <- df %>%
    filter(
      !is.na(location_long),
      !is.na(location_lat),
      is.finite(location_long),
      is.finite(location_lat)
    )
  n_removed <- n_before - nrow(df)
  
  if (n_removed > 0) {
    warning(sprintf("BUG: Final validation removed %d rows with invalid coordinates!", n_removed))
  }
  
  # Convert to character for CSV storage
  class_col(df)
}


#' Format Time Duration
#'
#' @param seconds Duration in seconds
#' @return Formatted string (e.g., "2h 15m 30s")
format_duration <- function(seconds) {
  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  secs <- round(seconds %% 60)
  
  if (hours > 0) {
    sprintf("%dh %dm %ds", hours, minutes, secs)
  } else if (minutes > 0) {
    sprintf("%dm %ds", minutes, secs)
  } else {
    sprintf("%ds", secs)
  }
}


#' Display Progress Bar
#'
#' @param current Current iteration
#' @param total Total iterations
#' @param kept Number kept so far
#' @param start_time Start time (POSIXct)
#' @param bar_width Width of progress bar
display_progress <- function(current, total, kept, start_time, bar_width = 40) {
  
  pct <- current / total
  filled <- round(pct * bar_width)
  empty <- bar_width - filled
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  eta_str <- if (current > 0) format_duration((elapsed / current) * (total - current)) else "..."
  
  bar <- paste0(
    "\r  [",
    paste(rep("=", filled), collapse = ""),
    if (filled < bar_width) ">" else "",
    paste(rep(" ", max(0, empty - 1)), collapse = ""),
    "] ",
    sprintf("%3.0f%% (%d/%d) | Kept: %d | ETA: %s    ", pct * 100, current, total, kept, eta_str)
  )
  
  cat(bar)
  flush.console()
}


# -----------------------------------------------------------------------------
# TEMPORARY FILE MANAGEMENT
# -----------------------------------------------------------------------------

#' Create Temporary Directory for Processing
#'
#' @param base_dir Base output directory
#' @param prefix Prefix for temp directory name
#' @return Path to created temp directory
create_temp_dir <- function(base_dir, prefix = "temp_processing") {
  # Handle empty or NULL base_dir
  if (is.null(base_dir) || base_dir == "" || base_dir == ".") {
    base_dir <- "output"
  }
  
  # Ensure base_dir exists
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Create unique temp directory name
  temp_dir <- file.path(base_dir, paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  temp_dir
}


#' Save Individual to Temporary File
#'
#' Saves processed individual data to a unique temp file.
#'
#' @param df Processed dataframe
#' @param temp_dir Temporary directory path
#' @param individual_id Individual identifier (used for filename)
#' @return Path to saved temp file, or NULL if df is empty
save_to_temp <- function(df, temp_dir, individual_id) {
  if (nrow(df) == 0) return(NULL)
  
  # Standardize columns before saving
  df <- standardize_output_columns(df)
  
  # Create safe filename
  safe_id <- gsub("[^a-zA-Z0-9_-]", "_", as.character(individual_id))
  temp_file <- file.path(temp_dir, paste0(safe_id, ".csv"))
  
  # Write without header (we'll add it when merging)
  write_csv(df, temp_file, col_names = TRUE, na = "")
  
  temp_file
}


#' Merge Temporary Files into Final CSV
#'
#' Reads all temp files and combines them into the final output CSV.
#'
#' @param temp_dir Directory containing temp files
#' @param output_csv Final output CSV path
#' @param cleanup Whether to remove temp directory after merging (default: TRUE)
#' @return Number of records written
merge_temp_files <- function(temp_dir, output_csv, cleanup = TRUE) {
  
  temp_files <- list.files(temp_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(temp_files) == 0) {
    message("  No temp files to merge.")
    return(0)
  }
  
  message(sprintf("  Merging %d temp files...", length(temp_files)))
  
  # Read and combine all temp files
  all_data <- lapply(temp_files, function(f) {
    tryCatch(
      read_csv(f, col_types = cols(.default = col_character()), progress = FALSE),
      error = function(e) {
        warning(sprintf("Failed to read temp file: %s", f))
        NULL
      }
    )
  })
  
  # Remove NULLs and combine
  all_data <- all_data[!sapply(all_data, is.null)]
  
  if (length(all_data) == 0) {
    message("  No valid data in temp files.")
    return(0)
  }
  
  combined_df <- bind_rows(all_data)
  
  # Ensure consistent column order
  combined_df <- standardize_output_columns(combined_df)
  
  # Append to output CSV (preserving any existing data from non-Movebank sources)
  append_to_csv(combined_df, output_csv)
  
  n_records <- nrow(combined_df)
  message(sprintf("  Merged %s records to output CSV.", format(n_records, big.mark = ",")))
  
  # Cleanup temp directory
  if (cleanup) {
    unlink(temp_dir, recursive = TRUE)
    message("  Cleaned up temp directory.")
  }
  
  n_records
}


# -----------------------------------------------------------------------------
# MOVEBANK PROCESSING (PARALLEL)
# -----------------------------------------------------------------------------

#' Process Single Movebank Individual (for parallel processing)
#'
#' Reads and processes a single individual file. Saves to temp file.
#'
#' @param file_path Path to individual CSV file
#' @param species Species name (for exclusion check)
#' @param excluded_species Species to exclude
#' @param area_res Restricted study area polygon
#' @param area_ext Extended study area polygon
#' @param coastline Coastline polygon
#' @param temp_dir Temporary directory for saving results
#' @param ... Additional parameters passed to process_individual
#' @return List with status and temp_file path (or NULL)
process_single_movebank <- function(file_path,
                                     species,
                                     excluded_species,
                                     area_res,
                                     area_ext,
                                     coastline,
                                     temp_dir,
                                     ...) {
  
  # Skip excluded species
  if (species %in% excluded_species) {
    return(list(status = "excluded_species", temp_file = NULL))
  }
  
  # Read individual file
  df <- tryCatch(
    read_movebank_individual(file_path),
    error = function(e) NULL
  )
  
  if (is.null(df) || nrow(df) == 0) {
    return(list(status = "read_error", temp_file = NULL))
  }
  
  # Skip if no altitude data
  if (!"height_above_msl" %in% names(df) || all(is.na(df$height_above_msl))) {
    return(list(status = "no_altitude", temp_file = NULL))
  }
  
  # Get individual ID for temp file naming
  individual_id <- unique(df$individual_id)[1]
  
  # Process individual
  df <- process_individual(
    df,
    area_res = area_res,
    area_ext = area_ext,
    coastline = coastline,
    ...
  )
  
  if (nrow(df) > 0) {
    # Save to temp file
    temp_file <- save_to_temp(df, temp_dir, individual_id)
    return(list(status = "kept", temp_file = temp_file))
  } else {
    return(list(status = "filtered", temp_file = NULL))
  }
}


#' Process All Movebank Individuals (Parallel with Temp Files)
#'
#' Iterates through the Movebank directory structure and processes each individual
#' using parallel processing. Saves to temp files, then merges at the end.
#'
#' @param base_dir Base directory (species/individual.csv structure)
#' @param output_csv Output CSV file path
#' @param excluded_species Species to exclude
#' @param area_res Restricted study area polygon
#' @param area_ext Extended study area polygon
#' @param coastline Coastline polygon
#' @param n_workers Number of parallel workers (default: availableCores() - 1)
#' @param batch_size Number of individuals per batch (default: 100)
#' @param ... Additional parameters passed to process_individual
#' @return List with total and kept counts
process_movebank_individuals <- function(base_dir,
                                          output_csv,
                                          excluded_species,
                                          area_res,
                                          area_ext,
                                          coastline,
                                          n_workers = NULL,
                                          batch_size = 100,
                                          ...) {
  
  # Get list of all individual files
  file_info <- list_movebank_files(base_dir)
  
  if (nrow(file_info) == 0) {
    message("  No individual files found in: ", base_dir)
    return(list(total = 0, kept = 0))
  }
  
  n_total <- nrow(file_info)
  n_species <- n_distinct(file_info$species)
  
  # Set up parallel workers
  if (is.null(n_workers)) {
    n_workers <- max(1, future::availableCores() - 1)
  }
  n_workers <- min(n_workers, n_total)  # Don't use more workers than files
  
  message(sprintf("  Found %d individuals across %d species", n_total, n_species))
  message(sprintf("  Using %d parallel workers\n", n_workers))
  
  # Create temp directory for this processing run
  temp_dir <- create_temp_dir(dirname(output_csv), "movebank_temp")
  message(sprintf("  Temp directory: %s\n", temp_dir))
  
  # Set up parallel backend
  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  
  # Set up progress reporting
  progressr::handlers(progressr::handler_progress(
    format = "  [:bar] :percent | :current/:total | ETA: :eta",
    width = 60,
    clear = FALSE
  ))
  
  # Track statistics
  stats <- list(
    n_kept = 0,
    n_excluded_species = 0,
    n_no_altitude = 0,
    n_filtered = 0,
    n_errors = 0
  )
  
  # Collect all temp file paths
  all_temp_files <- character(0)
  
  start_time <- Sys.time()
  
  # Process in batches
  n_batches <- ceiling(n_total / batch_size)
  
  progressr::with_progress({
    p <- progressr::progressor(steps = n_total)
    
    for (batch_idx in seq_len(n_batches)) {
      # Get indices for this batch
      start_idx <- (batch_idx - 1) * batch_size + 1
      end_idx <- min(batch_idx * batch_size, n_total)
      batch_info <- file_info[start_idx:end_idx, ]
      
      # Process batch in parallel - each worker saves to its own temp file
      results <- furrr::future_map(
        seq_len(nrow(batch_info)),
        function(i) {
          result <- process_single_movebank(
            file_path = batch_info$file_path[i],
            species = batch_info$species[i],
            excluded_species = excluded_species,
            area_res = area_res,
            area_ext = area_ext,
            coastline = coastline,
            temp_dir = temp_dir,
            ...
          )
          p()  # Update progress
          result
        },
        .options = furrr::furrr_options(seed = TRUE)
      )
      
      # Collect results and update stats
      for (result in results) {
        if (is.null(result)) {
          stats$n_errors <- stats$n_errors + 1
        } else if (result$status == "kept") {
          stats$n_kept <- stats$n_kept + 1
          if (!is.null(result$temp_file)) {
            all_temp_files <- c(all_temp_files, result$temp_file)
          }
        } else if (result$status == "excluded_species") {
          stats$n_excluded_species <- stats$n_excluded_species + 1
        } else if (result$status == "no_altitude") {
          stats$n_no_altitude <- stats$n_no_altitude + 1
        } else if (result$status == "filtered") {
          stats$n_filtered <- stats$n_filtered + 1
        } else {
          stats$n_errors <- stats$n_errors + 1
        }
      }
    }
  })
  
  cat("\n")
  
  # Merge all temp files into final CSV
  message("\n  Finalizing output...")
  merge_temp_files(temp_dir, output_csv, cleanup = TRUE)
  
  # Print summary
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  message("\n  Processing Summary:")
  message(sprintf("    Total: %d | Retained: %d (%.1f%%)", n_total, stats$n_kept, 100 * stats$n_kept / n_total))
  message(sprintf("    Excluded species: %d | No altitude: %d | Filtered: %d | Errors: %d",
                  stats$n_excluded_species, stats$n_no_altitude, stats$n_filtered, stats$n_errors))
  message(sprintf("    Time: %s (%.1f individuals/min)", format_duration(elapsed), 
                  n_total / (elapsed / 60)))
  
  list(total = n_total, kept = stats$n_kept)
}


# -----------------------------------------------------------------------------
# NON-MOVEBANK PROCESSING (SEQUENTIAL WITH TEMP FILES)
# -----------------------------------------------------------------------------

#' Process Non-Movebank Data Source
#'
#' Processes non-Movebank data sources sequentially, saving to temp files
#' before merging into the final output.
#'
#' @param reader_function Function to read the data source
#' @param data_dir Directory containing data
#' @param output_csv Output CSV file path
#' @param excluded_species Species to exclude
#' @param area_res Restricted study area polygon
#' @param area_ext Extended study area polygon
#' @param coastline Coastline polygon
#' @param downsample_interval Minutes between fixes (default: 10)
#' @param max_gap_days Days to split trajectories (default: 5)
#' @param min_points Minimum points per trajectory (default: 5)
#' @param dbscan_eps DBSCAN epsilon (default: 0.1)
#' @param dbscan_minpts DBSCAN minimum points (default: 10)
#' @param ... Ignored (for compatibility with shared proc_args)
#' @return List with total and kept counts
process_non_movebank_source <- function(reader_function,
                                         data_dir,
                                         output_csv,
                                         excluded_species,
                                         area_res,
                                         area_ext,
                                         coastline,
                                         downsample_interval = 10,
                                         max_gap_days = 5,
                                         min_points = 5,
                                         dbscan_eps = 0.1,
                                         dbscan_minpts = 10,
                                         ...) {
  
  df_all <- tryCatch(
    reader_function(data_dir),
    error = function(e) {
      message(sprintf("  Error reading data: %s", e$message))
      NULL
    }
  )
  
  if (is.null(df_all) || nrow(df_all) == 0) {
    message("  No data loaded.")
    return(list(total = 0, kept = 0))
  }
  
  message(sprintf("  Loaded %s records", format(nrow(df_all), big.mark = ",")))
  
  individuals <- unique(df_all$individual_id)
  n_total <- length(individuals)
  
  message(sprintf("  Processing %d individuals...\n", n_total))
  
  # Create temp directory for this source
  # Use data_dir basename for identification, sanitized for path safety
  source_name <- basename(data_dir)
  source_name <- gsub("[^a-zA-Z0-9_-]", "_", source_name)
  source_name <- substr(source_name, 1, 30)
  temp_dir <- create_temp_dir(dirname(output_csv), paste0("nonmb_temp_", source_name))
  
  n_kept <- 0
  n_excluded_species <- 0
  n_no_altitude <- 0
  n_filtered <- 0
  
  start_time <- Sys.time()
  
  for (i in seq_along(individuals)) {
    display_progress(i, n_total, n_kept, start_time)
    
    indiv <- individuals[i]
    df <- df_all %>% filter(individual_id == indiv)
    
    species <- unique(df$individual_taxon_canonical_name)
    if (length(species) > 0 && species %in% excluded_species) {
      n_excluded_species <- n_excluded_species + 1
      next
    }
    
    if (!"height_above_msl" %in% names(df) || all(is.na(df$height_above_msl))) {
      n_no_altitude <- n_no_altitude + 1
      next
    }
    
    df <- process_individual(
      df,
      area_res = area_res,
      area_ext = area_ext,
      coastline = coastline,
      downsample_interval = downsample_interval,
      max_gap_days = max_gap_days,
      min_points = min_points,
      dbscan_eps = dbscan_eps,
      dbscan_minpts = dbscan_minpts
    )
    
    if (nrow(df) > 0) {
      # Save to temp file instead of directly to output
      save_to_temp(df, temp_dir, indiv)
      n_kept <- n_kept + 1
    } else {
      n_filtered <- n_filtered + 1
    }
  }
  
  display_progress(n_total, n_total, n_kept, start_time)
  cat("\n\n")
  
  # Merge temp files into final CSV
  message("  Finalizing output...")
  merge_temp_files(temp_dir, output_csv, cleanup = TRUE)
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  message("  Processing Summary:")
  message(sprintf("    Total: %d | Retained: %d (%.1f%%)", n_total, n_kept, 100 * n_kept / n_total))
  message(sprintf("    Excluded species: %d | No altitude: %d | Filtered: %d",
                  n_excluded_species, n_no_altitude, n_filtered))
  message(sprintf("    Time: %s", format_duration(elapsed)))
  
  list(total = n_total, kept = n_kept)
}