# =============================================================================
# Visualization Functions
# =============================================================================
#
# Functions for creating summary plots and visualizations.
#
# =============================================================================


#' Load Ecological Groups Mapping
#'
#' @param ec_file Path to EC.xlsx file
#' @return Dataframe with species and ecological_group columns
load_ecological_groups <- function(ec_file = "data/supp/EC.xlsx") {
  if (!file.exists(ec_file)) stop("Ecological groups file not found: ", ec_file)
  
  read_excel(ec_file) %>%
    select(species = specie_scientific, ecological_group = ecological_group)
}


#' Add Ecological Groups to Data
#'
#' @param df Dataframe with individual_taxon_canonical_name column
#' @param ec_mapping Ecological groups dataframe
#' @return Dataframe with ecological_group column
add_ecological_groups <- function(df, ec_mapping) {
  group_translation <- c(
    "ArdÃ©idÃ©" = "Ardeid",
    "Grand Ã©chassier" = "Large wader",
    "LaridÃ©" = "Larid",
    "Limicole" = "Shorebird",
    "Terrestre" = "Terrestrial"
  )
  
  df %>%
    left_join(ec_mapping, by = c("individual_taxon_canonical_name" = "species")) %>%
    mutate(ecological_group = if_else(
      ecological_group %in% names(group_translation),
      group_translation[ecological_group],
      ecological_group
    ))
}


#' Load Processed GPS Data
#'
#' Loads processed data and filters any remaining NA coordinates.
#'
#' @param csv_path Path to processed CSV file
#' @return Dataframe with proper column types and valid coordinates only
load_processed_data <- function(csv_path) {
  df <- read_csv(
    csv_path,
    col_types = cols(.default = col_character()),
    progress = FALSE,
    na = c("", "NA", "NULL", "na", "null")
  ) %>%
    select(-starts_with("..."))
  
  # Convert to appropriate types
  df <- df %>%
    mutate(
      location_long = as.numeric(location_long),
      location_lat = as.numeric(location_lat),
      height_above_msl = as.numeric(height_above_msl),
      distance_to_shore = as.numeric(distance_to_shore),
      in_GOL = as.logical(in_GOL)
    )
  
  # Check for NA coordinates and filter them out
  n_na <- sum(is.na(df$location_long) | is.na(df$location_lat))
  if (n_na > 0) {
    warning(sprintf(
      "Found %d rows with NA coordinates. Filtering them out. (This indicates a processing bug that should be investigated.)",
      n_na
    ))
    df <- df %>%
      filter(!is.na(location_long), !is.na(location_lat))
  }
  
  df
}


#' Plot Gulf of Lion Migration Points
#'
#' @param gps GPS dataframe
#' @param eu European coastline sf polygon
#' @param state_filter Behavior state to plot (default: "Migration")
#' @return ggplot object
plot_gol_migration <- function(gps, eu, state_filter = "Migration") {
  plot_data <- gps %>%
    filter(in_GOL == TRUE, state == state_filter)
  
  if (nrow(plot_data) == 0) {
    return(ggplot() + theme_void() + labs(title = paste("No data for:", state_filter)))
  }
  
  plot_sf <- st_as_sf(plot_data, coords = c("location_long", "location_lat"), crs = 4326)
  
  ggplot() +
    geom_sf(data = eu, fill = "grey90") +
    geom_sf(data = plot_sf, aes(color = individual_id), size = 1.2, alpha = 0.7) +
    coord_sf(xlim = c(2, 7), ylim = c(41, 44), expand = FALSE) +
    facet_wrap(~ season, ncol = 2) +
    labs(x = "Longitude", y = "Latitude", title = paste("Gulf of Lion -", state_filter)) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
}


#' Plot Migration vs Staging Behavior
#'
#' @param gps GPS dataframe
#' @param eu European coastline sf polygon
#' @return ggplot object
plot_gol_migration_vs_staging <- function(gps, eu) {
  if (nrow(gps) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data available"))
  }
  
  plot_sf <- st_as_sf(gps, coords = c("location_long", "location_lat"), crs = 4326)
  
  ggplot() +
    geom_sf(data = eu, fill = "grey90") +
    geom_sf(data = plot_sf, aes(color = state), size = 1.2, alpha = 0.6) +
    coord_sf(xlim = c(2, 7), ylim = c(41, 44), expand = FALSE) +
    facet_wrap(~ state, ncol = 2) +
    scale_color_manual(values = c("Migration" = "#E69F00", "Staging" = "#56B4E9"), name = "Behavior") +
    labs(x = "Longitude", y = "Latitude", title = "Gulf of Lion - Migration vs Staging") +
    theme_classic(base_size = 12)
}


#' Plot by Ecological Group
#'
#' @param gps GPS dataframe with ecological_group column
#' @param eu European coastline sf polygon
#' @param season_filter Season to plot
#' @param state_filter Behavior state to plot
#' @return ggplot object
plot_gol_by_ecological_group <- function(gps, eu, season_filter = "pre-breeding", state_filter = "Migration") {
  plot_data <- gps %>%
    filter(in_GOL == TRUE, season == season_filter, state == state_filter) %>%
    filter(!is.na(ecological_group))
  
  if (nrow(plot_data) == 0) {
    return(ggplot() + theme_void() + labs(title = paste("No data for:", season_filter)))
  }
  
  plot_sf <- st_as_sf(plot_data, coords = c("location_long", "location_lat"), crs = 4326)
  
  ggplot() +
    geom_sf(data = eu, fill = "grey90") +
    geom_sf(data = plot_sf, aes(color = individual_id), size = 1.2, alpha = 0.7) +
    coord_sf(xlim = c(2, 7), ylim = c(41, 44), expand = FALSE) +
    facet_wrap(~ ecological_group, ncol = 3) +
    labs(x = "Longitude", y = "Latitude",
         title = paste("Gulf of Lion by Ecological Group -", season_filter)) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
}


#' Plot Flight Height Distribution
#'
#' @param gps GPS dataframe
#' @param state_filter Behavior state to plot
#' @return ggplot object
plot_flight_height_distribution <- function(gps, state_filter = "Migration") {
  plot_data <- gps %>%
    filter(in_GOL == TRUE, state == state_filter) %>%
    filter(!is.na(height_above_msl))
  
  if (nrow(plot_data) == 0) {
    return(ggplot() + theme_void() + labs(title = "No altitude data available"))
  }
  
  ggplot(plot_data, aes(x = height_above_msl)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(x = "Flight Height (m above MSL)", y = "Count",
         title = paste("Flight Height Distribution -", state_filter)) +
    theme_classic(base_size = 12)
}


#' Plot Flight Height by Ecological Group
#'
#' @param gps GPS dataframe
#' @param state_filter Behavior state to plot
#' @return ggplot object
plot_flight_height_by_ecological_group <- function(gps, state_filter = "Migration") {
  plot_data <- gps %>%
    filter(in_GOL == TRUE, state == state_filter) %>%
    filter(!is.na(ecological_group), !is.na(height_above_msl))
  
  if (nrow(plot_data) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data available"))
  }
  
  ggplot(plot_data, aes(x = height_above_msl)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    facet_wrap(~ ecological_group, ncol = 2, scales = "free_y") +
    labs(x = "Flight Height (m above MSL)", y = "Count",
         title = paste("Flight Height by Ecological Group -", state_filter)) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
}


#' Create All Visualizations
#'
#' @param gps GPS dataframe
#' @param coastline European coastline sf polygon
#' @param output_dir Directory to save plots
create_all_visualizations <- function(gps, coastline, output_dir) {
  
  message("    - Migration map (both seasons)...")
  p1 <- plot_gol_migration(gps, coastline)
  ggsave(file.path(output_dir, "01_gol_migration_both_seasons.png"), p1, width = 14, height = 7, dpi = 300)
  
  message("    - Migration vs Staging...")
  p2 <- plot_gol_migration_vs_staging(gps, coastline)
  ggsave(file.path(output_dir, "02_gol_migration_vs_staging.png"), p2, width = 10, height = 8, dpi = 300)
  
  if ("ecological_group" %in% names(gps) && any(!is.na(gps$ecological_group))) {
    message("    - By ecological group (pre-breeding)...")
    p3 <- plot_gol_by_ecological_group(gps, coastline, "pre-breeding")
    ggsave(file.path(output_dir, "03_gol_ecological_groups_pre_breeding.png"), p3, width = 14, height = 10, dpi = 300)
    
    message("    - By ecological group (post-breeding)...")
    p4 <- plot_gol_by_ecological_group(gps, coastline, "post-breeding")
    ggsave(file.path(output_dir, "04_gol_ecological_groups_post_breeding.png"), p4, width = 14, height = 10, dpi = 300)
  }
  
  message("    - Flight height distribution...")
  p5 <- plot_flight_height_distribution(gps)
  ggsave(file.path(output_dir, "05_flight_height_distribution.png"), p5, width = 10, height = 6, dpi = 300)
  
  if ("ecological_group" %in% names(gps) && any(!is.na(gps$ecological_group))) {
    message("    - Flight height by ecological group...")
    p6 <- plot_flight_height_by_ecological_group(gps)
    ggsave(file.path(output_dir, "06_flight_height_by_ecological_group.png"), p6, width = 12, height = 10, dpi = 300)
  }
  
  invisible(NULL)
}