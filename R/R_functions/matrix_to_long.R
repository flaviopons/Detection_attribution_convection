matrix_to_long <- function(matrix, lon, lat) {
  lon_lat <- expand.grid(lon = lon, lat = rev(lat))
  data.frame(lon_lat, value = as.vector(matrix))
}
