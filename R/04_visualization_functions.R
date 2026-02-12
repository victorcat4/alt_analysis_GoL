# =============================================================================
# Visualization Functions
# =============================================================================
#
# Functions for creating summary plots and visualizations.
#
# =============================================================================


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
         title = paste("Gulf of Lion by Group of Species -", season_filter)) +
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
         title = paste("Flight Height by Group of Species -", state_filter)) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
}


# -----------------------------------------------------------------------------
# DAY/NIGHT CLASSIFICATION
# -----------------------------------------------------------------------------

#' Classify GPS Points as Day or Night (French Time)
#'
#' Converts timestamps to French time (Europe/Paris) and uses the suncalc
#' package to determine sunrise/sunset at each point's location.
#' Points between sunrise and sunset are "Day", otherwise "Night".
#'
#' @param gps GPS dataframe with timestamp, location_lat, location_long
#' @return GPS dataframe with added "day_night" column
classify_day_night <- function(gps) {
  if (!require("suncalc", quietly = TRUE)) {
    install.packages("suncalc", dependencies = TRUE)
    library(suncalc)
  }
  
  # Convert UTC timestamps to French time (Europe/Paris)
  gps <- gps %>%
    mutate(
      timestamp_posix = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      timestamp_fr = with_tz(timestamp_posix, tzone = "Europe/Paris"),
      date = as.Date(timestamp_fr)
    )
  
  sun_times <- suncalc::getSunlightTimes(
    data = data.frame(
      date = gps$date,
      lat = as.numeric(gps$location_lat),
      lon = as.numeric(gps$location_long)
    ),
    keep = c("sunrise", "sunset"),
    tz = "Europe/Paris"
  )
  
  gps$day_night <- ifelse(
    gps$timestamp_fr >= sun_times$sunrise & gps$timestamp_fr <= sun_times$sunset,
    "Day", "Night"
  )
  
  gps <- gps %>% select(-timestamp_posix, -timestamp_fr, -date)
  
  gps
}


#' Plot Proportion of Points at Day vs Night
#'
#' Creates a bar chart showing the proportion of GPS points recorded during
#' day vs night, faceted by ecological group and season.
#'
#' @param gps GPS dataframe with in_GOL, state, timestamp, location columns
#' @param state_filter Behavior state to plot (default: "Migration")
#' @return ggplot object
plot_day_night_proportion <- function(gps, state_filter = "Migration") {
  
  plot_data <- gps %>%
    filter(in_GOL == TRUE, state == state_filter)
  
  if (nrow(plot_data) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data available"))
  }
  
  # Classify day/night
  message("      Computing sunrise/sunset times...")
  plot_data <- classify_day_night(plot_data)
  
  # Compute proportions
  has_eco <- "ecological_group" %in% names(plot_data) && any(!is.na(plot_data$ecological_group))
  
  if (has_eco) {
    plot_data <- plot_data %>% filter(!is.na(ecological_group))
    
    prop_data <- plot_data %>%
      count(ecological_group, season, day_night) %>%
      group_by(ecological_group, season) %>%
      mutate(
        total = sum(n),
        proportion = n / total * 100
      ) %>%
      ungroup()
    
    p <- ggplot(prop_data, aes(x = ecological_group, y = proportion, fill = day_night)) +
      geom_col(position = "stack", alpha = 0.85) +
      geom_text(
        aes(label = sprintf("%.0f%%", proportion)),
        position = position_stack(vjust = 0.5),
        size = 3, color = "white", fontface = "bold"
      ) +
      facet_wrap(~ season, ncol = 2) +
      scale_fill_manual(
        values = c("Day" = "#F4B942", "Night" = "#2C3E50"),
        name = "Period"
      ) +
      labs(
        x = "Ecological Group",
        y = "Proportion (%)",
        title = paste("Day vs Night Migration Points -", state_filter),
        subtitle = "Gulf of Lion region"
      ) +
      theme_classic(base_size = 12) +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1)
      )
    
  } else {
    prop_data <- plot_data %>%
      count(season, day_night) %>%
      group_by(season) %>%
      mutate(
        total = sum(n),
        proportion = n / total * 100
      ) %>%
      ungroup()
    
    p <- ggplot(prop_data, aes(x = season, y = proportion, fill = day_night)) +
      geom_col(position = "stack", alpha = 0.85) +
      geom_text(
        aes(label = sprintf("%.0f%%", proportion)),
        position = position_stack(vjust = 0.5),
        size = 4, color = "white", fontface = "bold"
      ) +
      scale_fill_manual(
        values = c("Day" = "#F4B942", "Night" = "#2C3E50"),
        name = "Period"
      ) +
      labs(
        x = "Season",
        y = "Proportion (%)",
        title = paste("Day vs Night Migration Points -", state_filter),
        subtitle = "Gulf of Lion region"
      ) +
      theme_classic(base_size = 12)
  }
  
  p
}


#' Create All Visualizations
#'
#' @param gps GPS dataframe
#' @param coastline European coastline sf polygon
#' @param output_dir Directory to save plots
create_all_visualizations <- function(gps, coastline, output_dir) {
  
  message("    - Migration map (both seasons)...")
  p1 <- plot_gol_migration(gps, coastline)
  ggsave(file.path(output_dir, "01_gol_migration_both_seasons.pdf"), p1, width = 14, height = 7, dpi = 300)
  
  if ("ecological_group" %in% names(gps) && any(!is.na(gps$ecological_group))) {
    message("    - By group of species (pre-breeding)...")
    p3 <- plot_gol_by_ecological_group(gps, coastline, "pre-breeding")
    ggsave(file.path(output_dir, "03_gol_ecological_groups_pre_breeding.pdf"), p3, width = 14, height = 10, dpi = 300)
    
    message("    - By group of species (post-breeding)...")
    p4 <- plot_gol_by_ecological_group(gps, coastline, "post-breeding")
    ggsave(file.path(output_dir, "04_gol_ecological_groups_post_breeding.pdf"), p4, width = 14, height = 10, dpi = 300)
  }
  
  message("    - Flight height distribution...")
  p5 <- plot_flight_height_distribution(gps)
  ggsave(file.path(output_dir, "05_flight_height_distribution.pdf"), p5, width = 10, height = 6, dpi = 300)
  
  if ("ecological_group" %in% names(gps) && any(!is.na(gps$ecological_group))) {
    message("    - Flight height by group of species...")
    p6 <- plot_flight_height_by_ecological_group(gps)
    ggsave(file.path(output_dir, "06_flight_height_by_ecological_group.pdf"), p6, width = 12, height = 10, dpi = 300)
  }
  
  message("    - Day vs Night proportions...")
  p7 <- plot_day_night_proportion(gps)
  ggsave(file.path(output_dir, "07_day_night_proportion.pdf"), p7, width = 12, height = 7, dpi = 300)
  
  message("    - Species summary table...")
  create_species_summary_table(gps, output_dir)
  
  message("    - Migration vs Staging...")
  p2 <- plot_gol_migration_vs_staging(gps, coastline)
  ggsave(file.path(output_dir, "02_gol_migration_vs_staging.pdf"), p2, width = 10, height = 8, dpi = 300)

  invisible(NULL)
}


#' Create Species Summary Table
#'
#' Generates a summary table with one row per species showing:
#' study IDs, total records, records in GOL, and records in GOL during migration.
#'
#' @param gps GPS dataframe
#' @param output_dir Directory to save the table
#' @return gt table object (invisible)
create_species_summary_table <- function(gps, output_dir) {
  
  if (!require("gt", quietly = TRUE)) {
    warning("Package 'gt' not available. Skipping species summary table.")
    return(invisible(NULL))
  }
  
  summary_df <- gps %>%
    group_by(individual_taxon_canonical_name) %>%
    summarise(
      study_ids = paste(sort(unique(study_id)), collapse = ", "),
      total_records = n(),
      records_in_GOL = sum(in_GOL == TRUE, na.rm = TRUE),
      records_in_GOL_migration = sum(in_GOL == TRUE & state == "Migration", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(total_records))
  
  summary_table <- summary_df %>%
    gt() %>%
    cols_label(
      individual_taxon_canonical_name = "Species",
      study_ids = "Study IDs",
      total_records = "Total Records",
      records_in_GOL = "In GOL",
      records_in_GOL_migration = "In GOL & Migration"
    ) %>%
    fmt_number(
      columns = c(total_records, records_in_GOL, records_in_GOL_migration),
      use_seps = TRUE,
      decimals = 0
    ) %>%
    tab_header(
      title = "Summary of GPS records by species"
    )
  
  gtsave(summary_table, file.path(output_dir, "08_species_summary_table.html"))
  message(sprintf("    Species summary table saved: %s",
                  file.path(output_dir, "08_species_summary_table.html")))
  
  invisible(summary_table)
}