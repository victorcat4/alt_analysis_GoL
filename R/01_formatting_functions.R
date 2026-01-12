# =============================================================================
# Formatting Functions
# =============================================================================
#
# Functions for cleaning and standardizing GPS tracking data.
#
# CRITICAL: clean_row() is the SINGLE point of coordinate validation.
# All downstream functions assume valid coordinates. NAs appearing later
# indicate a bug that must be fixed at the source.
#
# =============================================================================


# -----------------------------------------------------------------------------
# CSV UTILITIES
# -----------------------------------------------------------------------------

#' Append Dataframe to CSV
#'
#' @param df Dataframe to append
#' @param csv_path Path to CSV file
append_to_csv <- function(df, csv_path) {
  if (nrow(df) == 0) return(invisible(NULL))
  
  write_header <- !file.exists(csv_path)
  write_csv(df, csv_path, append = !write_header, col_names = write_header, na = "")
  invisible(NULL)
}


#' Get CSV Summary Statistics
#'
#' @param csv_path Path to CSV file
#' @return List with n_records, n_individuals, n_species, species_counts
get_csv_summary <- function(csv_path) {
  if (!file.exists(csv_path)) {
    stop("CSV file does not exist: ", csv_path)
  }
  
  summary_df <- read_csv(
    csv_path,
    col_select = c(individual_id, individual_taxon_canonical_name),
    col_types = cols(.default = col_character()),
    progress = FALSE
  )
  
  if (nrow(summary_df) == 0) {
    return(list(n_records = 0, n_individuals = 0, n_species = 0, species_counts = tibble()))
  }
  
  list(
    n_records = nrow(summary_df),
    n_individuals = n_distinct(summary_df$individual_id),
    n_species = n_distinct(summary_df$individual_taxon_canonical_name),
    species_counts = summary_df %>%
      count(individual_taxon_canonical_name, name = "n_records") %>%
      arrange(desc(n_records))
  )
}


# -----------------------------------------------------------------------------
# DATA CLEANING
# -----------------------------------------------------------------------------

#' Convert All Columns to Character for CSV Storage
#'
#' @param df A dataframe
#' @return Dataframe with all columns as character
class_col <- function(df) {
  df %>% mutate(across(everything(), as.character))
}


#' Clean Coordinate String
#'
#' Robust conversion of coordinate strings to numeric.
#' Handles whitespace, empty strings, NA strings, and comma decimals.
#'
#' @param x Character vector of coordinates
#' @return Numeric vector (invalid values become NA)
clean_coordinate <- function(x) {
  x <- as.character(x)
  
  # Trim whitespace
  x <- trimws(x)
  
  # Convert empty strings and common NA representations to NA
  x[x == ""] <- NA
  x[tolower(x) == "na"] <- NA
  x[tolower(x) == "null"] <- NA
  x[x == "-"] <- NA
  x[x == "."] <- NA
  
  # Handle comma decimal separators (European format)
  x <- gsub(",", ".", x)
  
  # Convert to numeric (invalid values become NA)
  suppressWarnings(as.numeric(x))
}


#' Remove Invalid Rows - SINGLE POINT OF COORDINATE VALIDATION
#'
#' Removes rows with missing or invalid coordinates/timestamps, duplicates,
#' and sorts by timestamp.
#'
#' THIS IS THE ONLY FUNCTION THAT REMOVES NA COORDINATES.
#' All downstream functions assume valid coordinates.
#'
#' @param df Dataframe with GPS data
#' @return Cleaned dataframe with numeric coordinates
clean_row <- function(df) {
  n_before <- nrow(df)
  
  # Convert coordinates using robust cleaning
  df <- df %>%
    mutate(
      location_long = clean_coordinate(location_long),
      location_lat = clean_coordinate(location_lat)
    )
  
  # Filter out invalid coordinates and timestamps
  df <- df %>%
    filter(
      # Remove NA values
      !is.na(location_long),
      !is.na(location_lat),
      !is.na(timestamp),
      # Remove infinite values
      is.finite(location_long),
      is.finite(location_lat),
      # Basic coordinate validity
      location_long >= -180,
      location_long <= 180,
      location_lat >= -90,
      location_lat <= 90
    ) %>%
    distinct() %>%
    arrange(timestamp)
  
  n_removed <- n_before - nrow(df)
  if (n_removed > 0) {
    pct_removed <- 100 * n_removed / n_before
    if (pct_removed > 1) {
      message(sprintf("      clean_row: removed %d invalid rows (%.1f%%)", n_removed, pct_removed))
    }
  }
  
  df
}


#' Validate Coordinates Are Clean
#'
#' Final validation check - returns TRUE if data has no NA coordinates.
#' Used as a sanity check before writing to CSV.
#'
#' @param df Dataframe with location_long and location_lat
#' @return TRUE if valid, FALSE otherwise (with warning)
validate_coordinates <- function(df) {
  if (nrow(df) == 0) return(TRUE)
  
  n_na <- sum(is.na(df$location_long) | is.na(df$location_lat))
  
  if (n_na > 0) {
    warning(sprintf("BUG: Found %d rows with NA coordinates after processing!", n_na))
    return(FALSE)
  }
  
  TRUE
}


#' Add Equipment Location
#'
#' Adds first recorded location as tagging site.
#'
#' @param df Dataframe with GPS data
#' @return Dataframe with "equip_loc" column
first_loc <- function(df) {
  df %>%
    group_by(individual_id) %>%
    mutate(
      equip_loc = paste(
        location_long[which.min(as.POSIXct(timestamp))],
        location_lat[which.min(as.POSIXct(timestamp))],
        sep = "_"
      )
    ) %>%
    ungroup()
}


# -----------------------------------------------------------------------------
# TEMPORAL PROCESSING
# -----------------------------------------------------------------------------

#' Downsample GPS Tracks
#'
#' Reduces temporal resolution by keeping middle point per interval.
#'
#' @param df Dataframe with GPS data for one individual
#' @param time_interval Interval in minutes (default: 30)
#' @return Downsampled dataframe
downsample_tracks <- function(df, time_interval = 30) {
  df %>%
    mutate(
      timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      interval = floor_date(timestamp, paste(time_interval, "minutes"))
    ) %>%
    group_by(individual_id, interval) %>%
    slice(ceiling(n() / 2)) %>%
    select(-interval) %>%
    ungroup()
}


#' Add Season Classification
#'
#' Classifies points as "pre-breeding" (Jan-Jun) or "post-breeding" (Jul-Dec).
#'
#' @param df Dataframe with timestamp column
#' @return Dataframe with "season" column
add_season <- function(df) {
  df %>%
    mutate(
      timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      season = if_else(month(timestamp) <= 6, "pre-breeding", "post-breeding")
    )
}


# -----------------------------------------------------------------------------
# TRAJECTORY ASSIGNMENT
# -----------------------------------------------------------------------------

#' Assign Trajectory IDs
#'
#' Splits trajectories at July 1st and gaps > max_gap_days.
#'
#' @param df Dataframe with individual_id and timestamp
#' @param max_gap_days Maximum gap in days (default: 5)
#' @return Dataframe with "trajectory" column
assign_trajectories <- function(df, max_gap_days = 5) {
  df %>%
    mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
    arrange(individual_id, timestamp) %>%
    group_by(individual_id) %>%
    mutate(
      year = year(timestamp),
      season_period = if_else(
        month(timestamp) < 7 | (month(timestamp) == 7 & day(timestamp) == 1),
        paste0(year, "_H1"),
        paste0(year, "_H2")
      ),
      time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "days")),
      new_trajectory = row_number() == 1 |
        season_period != lag(season_period) |
        (time_diff > max_gap_days),
      trajectory_num = cumsum(new_trajectory),
      trajectory = paste0(individual_id, "_traj_", trajectory_num)
    ) %>%
    select(-year, -season_period, -time_diff, -new_trajectory, -trajectory_num) %>%
    ungroup()
}


#' Remove Short Trajectories
#'
#' @param df Dataframe with trajectory column
#' @param min_points Minimum points per trajectory (default: 5)
#' @return Filtered dataframe
remove_short_trajectories <- function(df, min_points = 5) {
  df %>%
    group_by(trajectory) %>%
    filter(n() >= min_points) %>%
    ungroup()
}


# -----------------------------------------------------------------------------
# BEHAVIOR CLASSIFICATION
# -----------------------------------------------------------------------------

#' Classify Staging vs Migration
#'
#' Uses DBSCAN to classify points as "Staging" (clustered) or "Migration".
#'
#' @param df Dataframe with coordinates and trajectory
#' @param eps DBSCAN epsilon in degrees (default: 0.1)
#' @param minPts DBSCAN minimum points (default: 10)
#' @return Dataframe with "state" and "db_cluster" columns
classify_staging_migration <- function(df, eps = 0.1, minPts = 10) {
  df %>%
    group_by(trajectory) %>%
    group_modify(~ {
      if (nrow(.x) < minPts) {
        .x$db_cluster <- 0
        .x$state <- "Migration"
        return(.x)
      }
      
      coords <- .x %>% select(location_long, location_lat)
      clustering <- dbscan::dbscan(coords, eps = eps, minPts = minPts)
      
      .x$db_cluster <- clustering$cluster
      .x$state <- ifelse(clustering$cluster == 0, "Migration", "Staging")
      .x
    }) %>%
    ungroup()
}