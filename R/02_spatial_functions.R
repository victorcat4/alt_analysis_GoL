# =============================================================================
# Spatial Functions
# =============================================================================
#
# Functions for spatial filtering and analysis of GPS tracking data.
#
# =============================================================================


#' Load Spatial Reference Data
#'
#' @param data_dir Directory containing spatial data
#' @return List with coastline, extended_area, restricted_area
load_spatial_data <- function(data_dir = "data/supp") {
  list(
    coastline = readRDS(file.path(data_dir, "europe.rds")),
    extended_area = readRDS(file.path(data_dir, "GOL_extended.rds")),
    restricted_area = readRDS(file.path(data_dir, "GOL_restricted.rds"))
  )
}


#' Filter Points to Gulf of Lion Region
#'
#' Keeps only individuals that cross the restricted area,
#' retaining their points within the extended study area.
#'
#' @param df Dataframe with location_long, location_lat
#' @param area_res sf polygon - restricted area (must cross)
#' @param area_ext sf polygon - extended area (points retained)
#' @return Filtered dataframe (empty if no crossing)
select_pts <- function(df, area_res, area_ext) {
  if (nrow(df) == 0) return(df)
  
  df_sf <- st_as_sf(
    df,
    coords = c("location_long", "location_lat"),
    crs = 4326
  )
  
  crossing_bool <- any(st_intersects(df_sf, area_res, sparse = FALSE))
  
  if (crossing_bool) {
    area_bool <- st_intersects(df_sf, area_ext, sparse = FALSE)
    df[area_bool, ]
  } else {
    df[0, ]
  }
}


#' Check if Points are in Restricted Area
#'
#' @param df Dataframe with coordinates
#' @param area_res sf polygon - restricted area
#' @return Dataframe with "in_GOL" logical column
check_restricted_area <- function(df, area_res) {
  if (nrow(df) == 0) return(mutate(df, in_GOL = logical()))
  
  df_sf <- st_as_sf(
    df,
    coords = c("location_long", "location_lat"),
    crs = 4326
  )
  
  df$in_GOL <- as.vector(st_intersects(df_sf, area_res, sparse = FALSE))
  df
}


#' Calculate Distance to Shoreline
#'
#' @param df Dataframe with coordinates
#' @param coastline sf polygon - European coastline
#' @return Dataframe with "distance_to_shore" column (meters)
calculate_distance_to_shore <- function(df, coastline) {
  if (nrow(df) == 0) return(mutate(df, distance_to_shore = numeric()))
  
  df_sf <- st_as_sf(
    df,
    coords = c("location_long", "location_lat"),
    crs = 4326
  )
  
  # Transform to ETRS89-LAEA Europe for accurate distance
  df_projected <- st_transform(df_sf, crs = 3035)
  coastline_projected <- st_transform(coastline, crs = 3035)
  
  distances <- st_distance(df_projected, coastline_projected)
  df$distance_to_shore <- as.numeric(apply(distances, 1, min))
  
  df
}