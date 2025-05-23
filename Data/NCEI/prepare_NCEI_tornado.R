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

# read tornado data
NCEI_tornado <- read.csv("~/Desktop/ENGIE/ENGIE_Project/DATA/NCEI/1950-2023_torn.csv")
names(NCEI_tornado)
range(as.Date(NCEI_tornado$date)) #" "1955-01-03" "2023-12-23"

NCEI_tornado$mag[NCEI_tornado$mag<0]<-NA
hist(NCEI_tornado$mag)

#limit to CONUS
ind_removed <- c(which(NCEI_tornado$slat<25), which(NCEI_tornado$slat>50), 
                 which(NCEI_tornado$slon>-60), which(NCEI_tornado$slon<(-130)))
NCEI_tornado <- NCEI_tornado[-ind_removed,]


# obtain total number of yearly reports
NCEI_tornado <- NCEI_tornado %>%
  mutate(date = as.Date(date),
         ymon = format(date, "%Y%m")) 
annual_reports_tornado <- NCEI_tornado %>%
  group_by(yr) %>%
  summarise(total_reports = n())


# Aggregate daily data into monthly totals
monthly_reports_tornado <- NCEI_tornado %>%
  mutate(year_month = floor_date(date, "month")) %>% # Create a year-month column
  group_by(year_month) %>% # Group by year-month
  summarise(total_reports = n(), .groups = "drop") # Sum reports for each month

# Aggregate data into daily totals
daily_reports_tornado <- NCEI_tornado %>%
  group_by(date) %>%
  summarise(total_reports = n())

# annual reports by state
annual_reports_tornado_by_state <- NCEI_tornado %>%
  group_by(yr, st) %>%
  summarise(reports = n(), .groups = 'drop')

# SEASONALITY

# Aggregate daily data into monthly totals
seasonality_reports_tornado <- NCEI_tornado %>%
  group_by(mo) %>% # Group by month
  summarise(total_reports = n())
# Convert month numbers to month names
seasonality_reports_tornado <- seasonality_reports_tornado %>%
  mutate(month_name = factor(mo, levels = 1:12, labels = month.name))


# seasonal cycle by state
seasonality_tornado_by_state <- NCEI_tornado %>%
  group_by(mo, st) %>%
  summarise(reports = n(), .groups = 'drop')
# Convert month numbers to month names
seasonality_tornado_by_state <- seasonality_tornado_by_state %>%
  mutate(month_name = factor(mo, levels = 1:12, labels = month.name))


##### GRIDDED REPORTS
# using values from the database
#lon <- seq(from = min(floor(NCEI_tornado$slon)), to = ceiling(max(NCEI_tornado$slon)), by = 0.5)
#lat <- seq(from = min(floor(NCEI_tornado$slat)), to = max(ceiling(NCEI_tornado$slat)), by = 0.5)

# using the ERA5 grid
lon <- seq(from = -130, to = -60, by = 0.25)
lat <- seq(from = 25, to = 50, by = 0.25)
us_bbox <- c(xmin = min(lon), xmax = max(lon), ymin = min(lat), ymax = max(lat))  # Longitude and latitude

# Map tornado data onto the grid
grid <- expand.grid(lon = lon, lat = lat)
tornado_sf <- st_as_sf(NCEI_tornado, coords = c("slon", "slat"), crs = 4326)
grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
nearest_grid_index <- st_nearest_feature(tornado_sf, grid_sf)
NCEI_tornado$grid_index <- nearest_grid_index

# Calculate tornado statistics
tornado_counts <- NCEI_tornado %>%
  group_by(grid_index) %>%
  summarise(count = n())

grid_with_counts <- grid_sf %>%
  mutate(grid_index = row_number()) %>%
  left_join(tornado_counts, by = "grid_index") %>%
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
unique_dates <- sort(unique(NCEI_tornado$date))
TT <- length(unique_dates)

# Initialize 3D array (lon × lat × time)
tornado_array <- array(0, dim = c(length(lon), length(lat), TT),
                       dimnames = list(lon, lat, as.character(unique_dates)))

# Fill the array
for (i in 1:TT) {
  print(i)
  current_date <- unique_dates[i]
  daily_data <- NCEI_tornado %>% filter(date == current_date)
  
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
  tornado_array[,,i] <- matrix(grid_counts, nrow = length(lon), ncol = length(lat))
}

# Save as NetCDF file
save_as_netcdf <- function(array_3d, lon, lat, dates, filename) {
  # Define dimensions
  londim <- ncdim_def("longitude", "degrees_east", lon)
  latdim <- ncdim_def("latitude", "degrees_north", lat)
  timedim <- ncdim_def("time", "days since 1970-01-01", 
                       as.numeric(dates - as.Date("1970-01-01")))
  
  # Define variable
  tornado_def <- ncvar_def("tornado_reports", "count", 
                           list(londim, latdim, timedim), 
                           missval = -999, 
                           longname = "Daily tornado report counts")
  
  # Create NetCDF file
  ncout <- nc_create(filename, list(tornado_def))
  
  # Put the array
  ncvar_put(ncout, tornado_def, array_3d)
  
  # Add attributes
  ncatt_put(ncout, "longitude", "axis", "X")
  ncatt_put(ncout, "latitude", "axis", "Y")
  ncatt_put(ncout, "time", "axis", "T")
  ncatt_put(ncout, 0, "title", "Daily tornado reports gridded data")
  ncatt_put(ncout, 0, "source", "ESWD database")
  ncatt_put(ncout, 0, "grid_resolution", "0.25 degrees")
  
  nc_close(ncout)
}

# Save the data
save_as_netcdf(tornado_array, lon, lat, unique_dates, "tornado_reports_gridded_NCEI.nc")


###############################
###### PLOTS ##################
###############################

# plot time series of yearly report
ggplot(annual_reports_tornado, aes(x = yr, y = total_reports, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Yearly number of tornado reports",
       x = "Year",
       y = "Report Number") +
  theme_minimal() +
  scale_x_discrete(breaks = seq(min(annual_reports_tornado$yr), max(annual_reports_tornado$yr), by = 5)) + # Show every 5th year
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) # Rotate and reduce font size

# plot time series of monthly reports
ggplot(monthly_reports_tornado, aes(x = year_month, y = total_reports, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Monthly number of tornado reports",
       x = "Year",
       y = "Report Number") +
  theme_minimal() +
  scale_x_date(breaks = seq(min(monthly_reports_tornado$year_month), max(monthly_reports_tornado$year_month), by = "5 years"), # Show every year
               date_labels = "%Y") + # Format labels as years
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) # Rotate and reduce font size


# plot annual reports by state
ggplot(annual_reports_tornado_by_state, aes(x = yr, y = st, fill = reports)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +  # Usa una scala di colori
  labs(title = "Annual number of tornado reports by state",
       x = "Year",
       y = "State",
       fill = "Report number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot seasonal cycle
ggplot(seasonality_reports_tornado, aes(x = month_name, y = total_reports)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = total_reports), vjust = -0.5, color = "black") + # Add labels
  labs(title = "Seasonal Cycle of tornado Reports",
       x = "Month",
       y = "Total Reports") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# plot seasonal cycle by state
ggplot(seasonality_tornado_by_state, aes(x = month_name, y = st, fill = reports)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +  # Usa una scala di colori
  labs(title = "Annual number of tornado reports by state",
       x = "Month",
       y = "State",
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
  labs(title = "Total number of tornado Reports in the CONUS",
       x = "Longitude",
       y = "Latitude",
       fill = "tornado Reports Count") +
  theme_minimal()
