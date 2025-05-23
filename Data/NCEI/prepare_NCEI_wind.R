library(lubridate)
library(ncdf4)
library(raster)
library(sf)
library(dplyr)
library(ggplot2)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)

setwd('/home/fpons/Desktop/ENGIE/ENGIE_Project/DATA/NCEI')

# read wind data
NCEI_wind <- read.csv("~/Desktop/ENGIE/ENGIE_Project/DATA/NCEI/1955-2023_wind.csv")
names(NCEI_wind)
range(as.Date(NCEI_wind$date)) #" "1955-02-01" "2023-12-23"

# knots to km/h
NCEI_wind$mag_km <- NCEI_wind$mag*1.852
# 14196 of 53670 total reports have NA max diameter.
hist(NCEI_wind$mag_km)
summary(NCEI_wind$mag_km)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00   92.60   92.60   78.87  101.86  324.10 

#limit to CONUS
ind_removed <- c(which(NCEI_wind$slat<25), which(NCEI_wind$slat>50), 
                 which(NCEI_wind$slon>-60), which(NCEI_wind$slon<(-130)))
NCEI_wind <- NCEI_wind[-ind_removed,]


# obtain total number of yearly reports
NCEI_wind <- NCEI_wind %>%
  mutate(date = as.Date(date),
         ymon = format(date, "%Y%m")) 
annual_reports_wind <- NCEI_wind %>%
  group_by(yr) %>%
  summarise(total_reports = n())


# Aggregate daily data into monthly totals
monthly_reports_wind <- NCEI_wind %>%
  mutate(year_month = floor_date(date, "month")) %>% # Create a year-month column
  group_by(year_month) %>% # Group by year-month
  summarise(total_reports = n(), .groups = "drop") # Sum reports for each month

# Aggregate data into daily totals
daily_reports_wind <- NCEI_wind %>%
  group_by(date) %>%
  summarise(total_reports = n())

# annual reports by country
annual_reports_wind_by_country <- NCEI_wind %>%
  group_by(yr, st) %>%
  summarise(reports = n(), .groups = 'drop')

# SEASONALITY

# Aggregate daily data into monthly totals
seasonality_reports_wind <- NCEI_wind %>%
  group_by(mo) %>% # Group by month
  summarise(total_reports = n())
# Convert month numbers to month names
seasonality_reports_wind <- seasonality_reports_wind %>%
  mutate(month_name = factor(mo, levels = 1:12, labels = month.name))


# seasonal cycle by country
seasonality_wind_by_country <- NCEI_wind %>%
  group_by(mo, st) %>%
  summarise(reports = n(), .groups = 'drop')
# Convert month numbers to month names
seasonality_wind_by_country <- seasonality_wind_by_country %>%
  mutate(month_name = factor(mo, levels = 1:12, labels = month.name))


##### GRIDDED REPORTS
# using values from the database
#lon <- seq(from = min(floor(NCEI_wind$slon)), to = ceiling(max(NCEI_wind$slon)), by = 0.5)
#lat <- seq(from = min(floor(NCEI_wind$slat)), to = max(ceiling(NCEI_wind$slat)), by = 0.5)

# using the ERA5 grid
lon <- seq(from = -130, to = -60, by = 0.25)
lat <- seq(from = 25, to = 50, by = 0.25)
us_bbox <- c(xmin = min(lon), xmax = max(lon), ymin = min(lat), ymax = max(lat))  # Longitude and latitude

# Map wind data onto the grid
grid <- expand.grid(lon = lon, lat = lat)
wind_sf <- st_as_sf(NCEI_wind, coords = c("slon", "slat"), crs = 4326)
grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
nearest_grid_index <- st_nearest_feature(wind_sf, grid_sf)
NCEI_wind$grid_index <- nearest_grid_index

# Calculate wind statistics
wind_counts <- NCEI_wind %>%
  group_by(grid_index) %>%
  summarise(count = n())

grid_with_counts <- grid_sf %>%
  mutate(grid_index = row_number()) %>%
  left_join(wind_counts, by = "grid_index") %>%
  mutate(lon = st_coordinates(.)[, 1],  # Extract longitude
         lat = st_coordinates(.)[, 2])  # Extract latitude

#Create bins for counts
grid_with_counts <- grid_with_counts %>%
  mutate(count_bin = cut(count,
                         #                         breaks = c(0, 4, 9, 29, 49, 74, 99, 149, 299, Inf),
                         #                         labels = c("1-4", "5-9", "10-29", "30-49", "50-74", "75-99", "100-149", "150-299", "300+")))
                         breaks = c(0, 4, 9, 49, 99, 199, 349, 499, Inf),
                         labels = c("1-4", "5-9", "10-49", "50-99", "100-199", "200-349", "300-499", "500+")))



#######
# Get gridded daily reports, store in array and save to netcdf
# Get unique dates
unique_dates <- sort(unique(NCEI_wind$date))
TT <- length(unique_dates)

# Initialize 3D array (lon × lat × time)
wind_array <- array(0, dim = c(length(lon), length(lat), TT),
                    dimnames = list(lon, lat, as.character(unique_dates)))

# Fill the array
for (i in 1:TT) {
  print(i)
  current_date <- unique_dates[i]
  daily_data <- NCEI_wind %>% filter(date == current_date)
  
  # Count reports per grid cell for this date
  counts <- daily_data %>%
    group_by(grid_index) %>%
    summarise(count = n())
  
  # Convert to matrix format
  grid_counts <- grid_sf %>%
    mutate(grid_index = row_number()) %>%
    left_join(counts, by = "grid_index") %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    pull(count)
  
  # Reshape to lon×lat matrix and store in array
  wind_array[,,i] <- matrix(grid_counts, nrow = length(lon), ncol = length(lat))
}

# Save as NetCDF file
save_as_netcdf <- function(array_3d, lon, lat, dates, filename) {
  # Define dimensions
  londim <- ncdim_def("longitude", "degrees_east", lon)
  latdim <- ncdim_def("latitude", "degrees_north", lat)
  timedim <- ncdim_def("time", "days since 1970-01-01", 
                       as.numeric(dates - as.Date("1970-01-01")))
  
  # Define variable
  wind_def <- ncvar_def("wind_reports", "count", 
                        list(londim, latdim, timedim), 
                        missval = -999, 
                        longname = "Daily wind report counts")
  
  # Create NetCDF file
  ncout <- nc_create(filename, list(wind_def))
  
  # Put the array
  ncvar_put(ncout, wind_def, array_3d)
  
  # Add attributes
  ncatt_put(ncout, "longitude", "axis", "X")
  ncatt_put(ncout, "latitude", "axis", "Y")
  ncatt_put(ncout, "time", "axis", "T")
  ncatt_put(ncout, 0, "title", "Daily wind reports gridded data")
  ncatt_put(ncout, 0, "source", "ESWD database")
  ncatt_put(ncout, 0, "grid_resolution", "0.25 degrees")
  
  nc_close(ncout)
}

# Save the data
save_as_netcdf(wind_array, lon, lat, unique_dates, "wind_reports_gridded_NCEI.nc")


###############################
###### PLOTS ##################
###############################

# plot time series of yearly report
ggplot(annual_reports_wind, aes(x = yr, y = total_reports, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Yearly number of wind reports",
       x = "Year",
       y = "Report Number") +
  theme_minimal() +
  scale_x_discrete(breaks = seq(min(annual_reports_wind$yr), max(annual_reports_wind$yr), by = 5)) + # Show every 5th year
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) # Rotate and reduce font size

# plot time series of monthly reports
ggplot(monthly_reports_wind, aes(x = year_month, y = total_reports, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Monthly number of wind reports",
       x = "Year",
       y = "Report Number") +
  theme_minimal() +
  scale_x_date(breaks = seq(min(monthly_reports_wind$year_month), max(monthly_reports_wind$year_month), by = "5 years"), # Show every year
               date_labels = "%Y") + # Format labels as years
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) # Rotate and reduce font size


# plot annual reports by country
ggplot(annual_reports_wind_by_country, aes(x = yr, y = st, fill = reports)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +  # Usa una scala di colori
  labs(title = "Annual number of wind reports by country",
       x = "Year",
       y = "Country",
       fill = "Report number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot seasonal cycle
ggplot(seasonality_reports_wind, aes(x = month_name, y = total_reports)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = total_reports), vjust = -0.5, color = "black") + # Add labels
  labs(title = "Seasonal Cycle of wind Reports",
       x = "Month",
       y = "Total Reports") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# plot seasonal cycle by country
ggplot(seasonality_wind_by_country, aes(x = month_name, y = st, fill = reports)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +  # Usa una scala di colori
  labs(title = "Annual number of wind reports by state",
       x = "Month",
       y = "Country",
       fill = "Report number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot gridded total reports
world_map <- ne_countries(scale = "medium", returnclass = "sf")
region_map <- st_crop(world_map, us_bbox)

ggplot() +
  geom_tile(data = grid_with_counts, aes(x = lon, y = lat, fill = count_bin), color = NA) +
  geom_raster() +
  scale_fill_viridis_d(option = "inferno", na.value = "white", direction = -1,
                       labels = c("1-4", "5-9", "10-49", "50-99", "100-199", "200-349", "300-499", "500+"))+
  borders("world", colour = "black") +
  coord_sf(xlim = c(us_bbox["xmin"], us_bbox["xmax"]),
           ylim = c(us_bbox["ymin"], us_bbox["ymax"])) +
  labs(title = "Total number of wind Reports in the CONUS",
       x = "Longitude",
       y = "Latitude",
       fill = "wind Reports Count") +
  theme_minimal()
