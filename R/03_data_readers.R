# =============================================================================
# Data Reader Functions
# =============================================================================
#
# Functions to read GPS tracking data from various sources.
# These functions ONLY read and rename columns.
#
# Movebank Structure: base_dir/{species}/{individual}.csv
#
# =============================================================================


# Standard columns to retain
KEEP_COLS <- c(
  "individual_id",
  "study_id",
  "individual_taxon_canonical_name",
  "height_above_msl",
  "location_lat",
  "location_long",
  "timestamp",
  # GPS quality metrics (NAs allowed - not used for filtering)
  "gps_hdop",
  "gps_vdop",
  "gps_satellite_count"
)


# -----------------------------------------------------------------------------
# MOVEBANK DATA
# -----------------------------------------------------------------------------

#' List All Movebank Individual Files
#'
#' @param base_dir Base directory containing species folders
#' @return Dataframe with file_path, species, individual_id
list_movebank_files <- function(base_dir) {
  if (!dir.exists(base_dir)) {
    stop("Movebank directory not found: ", base_dir)
  }
  
  species_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  
  file_list <- lapply(species_dirs, function(sp_dir) {
    # Convert underscores to spaces to match taxonomic canonical names
    species_name <- gsub("_", " ", basename(sp_dir))
    csv_files <- list.files(sp_dir, pattern = "\\.csv$", full.names = TRUE)
    
    if (length(csv_files) == 0) return(NULL)
    
    tibble(
      file_path = csv_files,
      species = species_name,
      individual_id = tools::file_path_sans_ext(basename(csv_files))
    )
  })
  
  bind_rows(file_list)
}


#' Read Single Movebank Individual File
#'
#' @param file_path Path to individual CSV file
#' @param keep_cols Columns to retain
#' @return Dataframe with standardized columns
read_movebank_individual <- function(file_path, keep_cols = KEEP_COLS) {
  read_csv(
    file_path,
    col_types = cols(.default = col_character()),
    progress = FALSE
  ) %>%
    select(any_of(keep_cols))
}


# -----------------------------------------------------------------------------
# NON-MOVEBANK SOURCES
# -----------------------------------------------------------------------------

#' Read Turtle Dove & Stone Curlew Data
#'
#' @param data_dir Directory containing data files
#' @return Dataframe with standardized columns
read_turtle_dove_stone_curlew <- function(data_dir) {
  if (!dir.exists(data_dir)) stop("Directory not found: ", data_dir)
  
  files <- c(
    "CentTourterelles_location_data_all.csv",
    "Oedicnemes10_location_data_all.csv",
    "Oedicnemes30_location_data_all.csv",
    "Oedicnemes68_location_data_all.csv",
    "Oedicnemes69_location_data_all.csv",
    "Oedicnemes71_location_data_all.csv"
  )
  
  file_paths <- file.path(data_dir, files)
  file_paths <- file_paths[file.exists(file_paths)]
  
  df_list <- lapply(file_paths, function(f) {
    read.csv(f, sep = ";", stringsAsFactors = FALSE) %>%
      tibble() %>%
      mutate(
        indiv_long = gsub(",", ".", as.character(indiv_long)),
        indiv_lati = gsub(",", ".", as.character(indiv_lati)),
        indiv_elvt = gsub(",", ".", as.character(indiv_elvt))
      ) %>%
      transmute(
        individual_id = indiv_name,
        location_long = indiv_long,
        location_lat = indiv_lati,
        height_above_msl = indiv_elvt,
        timestamp = indiv_time,
        study_id = "operation100tourterelles_oedicnemes",
        individual_taxon_canonical_name = if_else(
          grepl("Tourterelles", f),
          "Streptopelia turtur",
          "Burhinus oedicnemus"
        )
      )
  })
  
  bind_rows(df_list)
}


#' Read Honey Buzzard Data (Patrick Byholm)
#'
#' @param data_dir Directory containing data files
#' @return Dataframe with standardized columns
read_honey_buzzard_byholm <- function(data_dir) {
  if (!dir.exists(data_dir)) stop("Directory not found: ", data_dir)
  
  df_all <- tibble()
  
  # File 1: 6056.csv
  file1 <- file.path(data_dir, "6056.csv")
  if (file.exists(file1)) {
    df <- read_csv(file1, show_col_types = FALSE) %>%
      transmute(
        individual_id = as.character(device_info_serial),
        timestamp = date_time,
        location_lat = latitude,
        location_long = longitude,
        height_above_msl = altitude,
        individual_taxon_canonical_name = "Pernis apivorus",
        study_id = "patrick_byholm_data"
      )
    df_all <- bind_rows(df_all, df)
  }
  
  # File 2: YoungBirds.csv
  file2 <- file.path(data_dir, "YoungBirds.csv")
  if (file.exists(file2)) {
    df <- read_csv(file2, show_col_types = FALSE) %>%
      transmute(
        individual_id = `individual-local-identifier`,
        individual_taxon_canonical_name = `individual-taxon-canonical-name`,
        location_lat = `location-lat`,
        location_long = `location-long`,
        timestamp = timestamp,
        study_id = "patrick_byholm_data"
      )
    df_all <- bind_rows(df_all, df)
  }
  
  # File 3: AdultBird(6056).csv
  file3 <- file.path(data_dir, "AdultBird(6056).csv")
  if (file.exists(file3)) {
    df <- read_csv(file3, show_col_types = FALSE) %>%
      transmute(
        individual_id = `individual-local-identifier`,
        individual_taxon_canonical_name = `individual-taxon-canonical-name`,
        location_lat = `location-lat`,
        location_long = `location-long`,
        timestamp = timestamp,
        study_id = "patrick_byholm_data"
      )
    df_all <- bind_rows(df_all, df)
  }
  
  df_all
}


#' Read BIRDMAP Data
#'
#' @param data_dir Directory containing data file
#' @param file_name File name (default: "occurrence.txt")
#' @return Dataframe with standardized columns
read_birdmap <- function(data_dir, file_name = "occurrence.txt") {
  file_path <- file.path(data_dir, file_name)
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  
  read_delim(file_path, delim = "\t", show_col_types = FALSE) %>%
    select(where(~ any(!is.na(.)))) %>%
    transmute(
      individual_id = acceptedNameUsageID,
      study_id = datasetID,
      individual_taxon_canonical_name = species,
      location_lat = decimalLatitude,
      location_long = decimalLongitude,
      timestamp = eventDate
    )
}


#' Read Montagu's Harrier Data
#'
#' @param data_dir Directory containing data file
#' @param file_name File name
#' @return Dataframe with standardized columns
read_montagu_harrier <- function(data_dir,
                                  file_name = "Schlaich_etal_2016_JAE_gps_tracking_data_20160802.csv") {
  file_path <- file.path(data_dir, file_name)
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  
  read_csv(file_path, show_col_types = FALSE) %>%
    transmute(
      location_lat = latitude,
      location_long = longitude,
      timestamp = date_time,
      individual_id = "Schlaich_etal_2016_id1",
      study_id = "Schlaich_etal_2016",
      individual_taxon_canonical_name = "Circus pygargus"
    )
}