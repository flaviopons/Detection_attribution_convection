create_gridded_3d_array_eswd <- function(eswd_data, minlon, maxlon, minlat, maxlat, res) {
  # Load required packages
  require(dplyr)
  require(tidyr)
  
  # 1. Create grid function
  create_grid <- function(lonmin, lonmax, latmin, latmax, res) {
    expand.grid(
      lon = seq(lonmin + res/2, lonmax - res/2, by = res),
      lat = seq(latmin + res/2, latmax - res/2, by = res)
    )
  }
  
  # 2. Create the grid
  grid <- create_grid(minlon, maxlon, minlat, maxlat, res)
  
  # 3. Function to assign points to grid cells
  assign_to_grid <- function(lon, lat, grid) {
    distances <- sqrt((grid$lon - lon)^2 + (grid$lat - lat)^2)
    which.min(distances)
  }
  
  # 4. Process ESWD data into gridded format
  process_data <- function(data, grid) {
    # Assign to grid cells
    data$grid_idx <- mapply(
      function(lon, lat) assign_to_grid(lon, lat, grid),
      data$LONGITUDE,
      data$LATITUDE
    )
    
    # Add grid coordinates
    data <- data %>%
      mutate(
        grid_lon = grid$lon[grid_idx],
        grid_lat = grid$lat[grid_idx]
      )
    
    # Calculate daily counts
    gridded_counts <- data %>%
      group_by(date, grid_lon, grid_lat) %>%
      summarise(count = n(), .groups = "drop")
    
    # Create complete grid with all combinations
    all_dates_grid <- expand.grid(
      date = unique(data$date),
      grid_lon = unique(grid$lon),
      grid_lat = unique(grid$lat)
    )
    
    # Join and replace NAs with 0
    left_join(all_dates_grid, gridded_counts, by = c("date", "grid_lon", "grid_lat")) %>%
      mutate(count = ifelse(is.na(count), 0, count))
  }
  
  # 5. Convert to 3D array
  convert_to_3d <- function(gridded_data) {
    # Get unique dates for time dimension
    unique_dates <- sort(unique(gridded_data$date))
    
    # Create array
    array_3d <- gridded_data %>%
      pivot_wider(names_from = date, values_from = count) %>%
      {array(
        data = as.matrix(.[,-c(1,2)]),
        dim = c(
          length(unique(.$grid_lon)),
          length(unique(.$grid_lat)),
          length(unique_dates)
        ),
        dimnames = list(
          sort(unique(.$grid_lon)),
          sort(unique(.$grid_lat)),
          as.character(unique_dates)
        )
      )}
    
    return(list(array_3d = array_3d, dates = unique_dates))
  }
  
  # Process the data
  gridded_data <- process_data(eswd_data, grid)
  
  # Convert to 3D array and return both array and dates
  result <- convert_to_3d(gridded_data)
  return(result)
}

# Example usage:
#result <- create_gridded_3d_array(hail_data, -15, 25, 35, 60, 0.25)
#grid_array <- result$array_3d  # The 3D array [lon × lat × time]
#date_vector <- result$dates    # Vector of dates corresponding to time dimension
