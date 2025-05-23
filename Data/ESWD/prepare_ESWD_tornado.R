library(lubridate)
library(ncdf4)
library(raster)
library(sf)
library(dplyr)
library(ggplot2)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)

setwd('/home/fpons/Desktop/ENGIE/ENGIE_Project/DATA/ESWD')

# read tornado data
ESWD_TO_2000_2024 <- read.csv("~/Desktop/ENGIE/ENGIE_Project/DATA/ESWD/17265135707_v1_7_ENGIE-Pons_ESWD domain_TO_20000101_20240706_QC1-2.csv")
names(ESWD_TO_2000_2024)

# check which event ID are repeated more than once
# Count occurrences of each ID
id_counts <- ESWD_TO_2000_2024 %>% count(ID)
# Filter IDs that appear more than once
duplicate_ids <- id_counts %>% filter(n > 1) %>% pull(ID)
# Extract rows with duplicate IDs
duplicate_rows <- ESWD_TO_2000_2024 %>%
  filter(ID %in% duplicate_ids)
# If you look at this sorting by ID, it is clear that events with the same ID are years apart
# ID cannot be used to characterize different events, variable to ignore.

range(as.Date(ESWD_TO_2000_2024$TIME_EVENT)) #"2000-01-13" "2024-07-06"

# Impact "I1" refers explicitly to energy infrastructure, however only 6 reports for tornadoes
# have this code.
length(which(ESWD_TO_2000_2024$IMPACTS=="I1"))

unique(ESWD_TO_2000_2024$SURFACE_CROSSED)
# We separate tornadoes from tornados that only stayed on water
#[1] "WATER,LAND"         "LAND"               "WATER"              "RURAL"              "LAKE"              
#[6] "SEA"                "RURAL,FOREST"       "FOREST"             "URBAN"              "LAND,SEA"          
#[11] "RURAL,FOREST,GRASS"

# separate tornadoes that happen on land, tornadoes and landfalling tornadoes
ind_tornado <- which(ESWD_TO_2000_2024$SURFACE_INITIAL_LOCATION == "SEA" |
                          ESWD_TO_2000_2024$SURFACE_INITIAL_LOCATION == "WATER" |
                          ESWD_TO_2000_2024$SURFACE_CROSSED == "WATER" | 
                          ESWD_TO_2000_2024$SURFACE_CROSSED == "SEA" | 
                          ESWD_TO_2000_2024$SURFACE_CROSSED == "WATER,LAND" |
                          ESWD_TO_2000_2024$SURFACE_CROSSED == "LAND,SEA")

tornado_data <- ESWD_TO_2000_2024[-ind_tornado,]
tornado_data$MONTH <- month(tornado_data$TIME_EVENT)

## there are data in all of Russia, including Siberia. Also, some reports have unrealistic latitude.
## Remove anything below 30 degrees N and further East than 50E
##ind_removed <- c(which(ESWD_LH_2000_2024$LATITUDE<30), which(ESWD_LH_2000_2024$LATITUDE>72), 
##                 which(ESWD_LH_2000_2024$LONGITUDE>50), which(ESWD_LH_2000_2024$LONGITUDE<(-25.5)))

# Remove anything outside the ERA5 domain
ind_removed <- c(which(ESWD_TO_2000_2024$LATITUDE<35), which(ESWD_TO_2000_2024$LATITUDE>60), 
                 which(ESWD_TO_2000_2024$LONGITUDE>25), which(ESWD_TO_2000_2024$LONGITUDE<(-15)))

tornado_data <- tornado_data[-ind_removed,]

## DAILY, MONTHLY AND YEARLY NUMBER OF EVENTS

# obtain total number of yearly reports
tornado_data <- tornado_data %>%
  mutate(date = as.Date(TIME_EVENT),
         year = format(date, "%Y")) 
annual_reports_tornado <- tornado_data %>%
  group_by(year) %>%
  summarise(total_reports = n())

# Aggregate daily data into monthly totals
monthly_reports_tornado <- tornado_data %>%
  mutate(year_month = floor_date(date, "month")) %>% # Create a year-month column
  group_by(year_month) %>% # Group by year-month
  summarise(total_reports = n(), .groups = "drop") # Sum reports for each month

# Aggregate data into daily totals
daily_reports_tornado <- tornado_data %>%
  group_by(date) %>%
  summarise(total_reports = n())

# annual reports by country
annual_reports_tornado_by_country <- tornado_data %>%
  group_by(year, COUNTRY) %>%
  summarise(reports = n(), .groups = 'drop')

## SEASONALITY

tornado_data$MONTH <- month(tornado_data$TIME_EVENT)

# Aggregate daily data into monthly totals
seasonality_reports_tornado <- tornado_data %>%
  group_by(MONTH) %>% # Group by year-month
  summarise(total_reports = n())
# Convert month numbers to month names
seasonality_reports_tornado <- seasonality_reports_tornado %>%
  mutate(month_name = factor(MONTH, levels = 1:12, labels = month.name))

# seasonal cycle by country
seasonality_tornado_by_country <- tornado_data %>%
  group_by(MONTH, COUNTRY) %>%
  summarise(reports = n(), .groups = 'drop')
# Convert month numbers to month names
seasonality_tornado_by_country <- seasonality_tornado_by_country %>%
  mutate(month_name = factor(MONTH, levels = 1:12, labels = month.name))



# RUSSIA HAS A HUGE NUMBER OF REPORTS AROUND 2017 WHICH HIDES SIGNAL FROM OTHER COUNTRIES
#Check what happened
tornado_data_russia <- tornado_data[which(tornado_data$COUNTRY == 'RU'),]

daily_reports_russia <- tornado_data_russia %>%
  group_by(date) %>%
  summarise(total_reports = n())


# CREATE A VERSION WITHOUT RUSSIA
tornado_data_norussia <- tornado_data[-which(tornado_data$COUNTRY == 'RU'),]

annual_reports_tornado_norussia <- tornado_data_norussia %>%
  group_by(year) %>%
  summarise(total_reports = n())

monthly_reports_tornado_norussia <- tornado_data_norussia %>%
  mutate(year_month = floor_date(date, "month")) %>% # Create a year-month column
  group_by(year_month) %>% # Group by year-month
  summarise(total_reports = n(), .groups = "drop") # Sum reports for each month

daily_reports_tornado_norussia <- tornado_data_norussia %>%
  group_by(date) %>%
  summarise(total_reports = n())


## SEASONALITY

tornado_data$MONTH <- month(tornado_data$TIME_EVENT)

# Aggregate daily data into monthly totals
seasonality_reports_tornado_norussia <- tornado_data_norussia %>%
  group_by(MONTH) %>% # Group by year-month
  summarise(total_reports = n())
# Convert month numbers to month names
seasonality_reports_tornado_norussia <- seasonality_reports_tornado_norussia %>%
  mutate(month_name = factor(MONTH, levels = 1:12, labels = month.name))

# seasonal cycle by country
seasonality_tornado_by_country_norussia <- tornado_data_norussia %>%
  group_by(MONTH, COUNTRY) %>%
  summarise(reports = n(), .groups = 'drop')
# Convert month numbers to month names
seasonality_tornado_by_country_norussia <- seasonality_tornado_by_country_norussia %>%
  mutate(month_name = factor(MONTH, levels = 1:12, labels = month.name))

# annual reports by country
annual_reports_tornado_by_country_norussia <- tornado_data_norussia %>%
  group_by(year, COUNTRY) %>%
  summarise(reports = n(), .groups = 'drop')


##### GRIDDED REPORTS
# using values from the database
#lon <- seq(from = min(floor(tornado_data$LONGITUDE)), to = ceiling(max(tornado_data$LONGITUDE)), by = 0.25)
#lat <- seq(from = min(floor(tornado_data$LATITUDE)), to = max(ceiling(tornado_data$LATITUDE)), by = 0.25)

# using the ERA5 grid
lon <- seq(from = -15, to = 25, by = 0.25)
lat <- seq(from = 35, to = 60, by = 0.25)
europe_bbox <- c(xmin = min(lon), xmax = max(lon), ymin = min(lat), ymax = max(lat))  # Longitude and latitude

# Map tornado data onto the grid
grid <- expand.grid(lon = lon, lat = lat)
tornado_sf <- st_as_sf(tornado_data_norussia, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
nearest_grid_index <- st_nearest_feature(tornado_sf, grid_sf)
tornado_data_norussia$grid_index <- nearest_grid_index

# Calculate tornado statistics
tornado_counts <- tornado_data_norussia %>%
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
                         breaks = c(0, 2, 4, 6, 8, 10, 15, 19, Inf),
                         labels = c("1-2", "3-4", "5-6", "7-8", "9-10", "11-15", "16-19", "20+")))


#######
# Get gridded daily reports, store in array and save to netcdf

# Get unique dates
unique_dates <- sort(unique(tornado_data_norussia$date))
TT <- length(unique_dates)

# Initialize 3D array (lon × lat × time)
tornado_array <- array(0, dim = c(length(lon), length(lat), TT),
                    dimnames = list(lon, lat, as.character(unique_dates)))

# Fill the array
for (i in 1:TT) {
  print(i)
  current_date <- unique_dates[i]
  daily_data <- tornado_data_norussia %>% filter(date == current_date)
  
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
save_as_netcdf(tornado_array, lon, lat, unique_dates, "tornado_reports_gridded_ESWD.nc")


###############################
###### PLOTS ##################
###############################

# plot time series of yearly reports
ggplot(annual_reports_tornado, aes(x = year, y = total_reports, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Yearly number of tornado reports",
       x = "Year",
       y = "Report Number") +
  theme_minimal()+
  scale_x_discrete(breaks = seq(min(annual_reports_tornado$year), max(annual_reports_tornado$year), by = 5)) + # Show every 5th year
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


# plot time series of daily reports
ggplot(daily_reports_tornado, aes(x = date, y = total_reports)) +
  geom_line(color = "blue") +
  geom_point(color = "red", size = 1) +
  labs(title = "Daily number of tornado reports",
       x = "Date",
       y = "Report Number") +
  theme_minimal() +
  # Set x-axis to show one label per year
  scale_x_date(
    date_breaks = "1 year",       # Major breaks every year
    date_labels = "%Y",         # only show the year
    minor_breaks = "1 month",       # Optional: show minor gridlines monthly
    expand = c(0.01, 0.01)        # Reduce padding on sides
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor.x = element_line(color = "grey90", linewidth = 0.1)
  )


# plot annual number of reports by country
ggplot(annual_reports_tornado_by_country, aes(x = year, y = COUNTRY, fill = reports)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +  # Usa una scala di colori
  labs(title = "Annual number of tornado reports by country",
       x = "Year",
       y = "Country",
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

# plot seasonal cycle by country
ggplot(seasonality_tornado_by_country, aes(x = month_name, y = COUNTRY, fill = reports)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +  # Usa una scala di colori
  labs(title = "Annual number of tornado reports by country",
       x = "Month",
       y = "Country",
       fill = "Report number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot time series of daily Russia reports
ggplot(daily_reports_russia, aes(x = date, y = total_reports, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Daily number of tornado reports",
       x = "Year",
       y = "Report Number") +
  theme_minimal()


# AS THE PLOTS ABOVE, BUT WITHOUT RUSSIA

# plot time series of yearly reports
ggplot(annual_reports_tornado_norussia, aes(x = year, y = total_reports, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Yearly number of tornado reports",
       x = "Year",
       y = "Report Number") +
  theme_minimal()+
  scale_x_discrete(breaks = seq(min(annual_reports_tornado_norussia$year), max(annual_reports_tornado_norussia$year), by = 5)) + # Show every 5th year
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) # Rotate and reduce font size


# plot time series of monthly reports
ggplot(monthly_reports_tornado_norussia, aes(x = year_month, y = total_reports, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Monthly number of tornado reports",
       x = "Year",
       y = "Report Number") +
  theme_minimal() +
  scale_x_date(breaks = seq(min(monthly_reports_tornado_norussia$year_month), max(monthly_reports_tornado_norussia$year_month), by = "5 years"), # Show every year
               date_labels = "%Y") + # Format labels as years
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) # Rotate and reduce font size


# plot time series of daily reports
ggplot(daily_reports_tornado_norussia, aes(x = date, y = total_reports)) +
  geom_line(color = "blue") +
  geom_point(color = "red", size = 1) +
  labs(title = "Daily number of tornado reports",
       x = "Date",
       y = "Report Number") +
  theme_minimal() +
  # Set x-axis to show one label per year
  scale_x_date(
    date_breaks = "1 year",       # Major breaks every year
    date_labels = "%Y",         # only show the year
    minor_breaks = "1 month",       # Optional: show minor gridlines monthly
    expand = c(0.01, 0.01)        # Reduce padding on sides
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor.x = element_line(color = "grey90", linewidth = 0.1)
  )


# plot annual number of reports by country
ggplot(annual_reports_tornado_by_country_norussia, aes(x = year, y = COUNTRY, fill = reports)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +  # Usa una scala di colori
  labs(title = "Annual number of tornado reports by country",
       x = "Year",
       y = "Country",
       fill = "Report number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# plot seasonal cycle
ggplot(seasonality_reports_tornado_norussia, aes(x = month_name, y = total_reports)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = total_reports), vjust = -0.5, color = "black") + # Add labels
  labs(title = "Seasonal Cycle of tornado Reports",
       x = "Month",
       y = "Total Reports") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot seasonal cycle by country
ggplot(seasonality_tornado_by_country_norussia, aes(x = month_name, y = COUNTRY, fill = reports)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +  # Usa una scala di colori
  labs(title = "Annual number of tornado reports by country",
       x = "Month",
       y = "Country",
       fill = "Report number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot gridded total reports
world_map <- ne_countries(scale = "medium", returnclass = "sf")
region_map <- st_crop(world_map, europe_bbox)

ggplot() +
  geom_tile(data = grid_with_counts, aes(x = lon, y = lat, fill = count_bin), color = NA) +
  geom_raster()+
  scale_fill_viridis_d(option = "plasma", na.value = "white", direction = -1,
                       labels = c("1-2", "3-4", "5-6", "7-8", "9-10", "11-15", "16-19", "20+", "")) +
  borders("world", colour = "black") +
  coord_sf(xlim = c(europe_bbox["xmin"], europe_bbox["xmax"]),
           ylim = c(europe_bbox["ymin"], europe_bbox["ymax"])) +
  labs(title = "Tornado Reports on Regular 0.25-Degree Grid (Europe)",
       x = "Longitude",
       y = "Latitude",
       fill = "Tornado Count") +
  theme_minimal()


