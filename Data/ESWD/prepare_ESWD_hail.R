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

# read hail data
ESWD_LH_2000_2024 <- read.csv("~/Desktop/ENGIE/ENGIE_Project/DATA/ESWD/1735497550_v1_7_ENGIE-Pons_ESWD domain_LH_20000101_20240706_QC1-2.csv")
names(ESWD_LH_2000_2024)

# check which event ID are repeated more than once
# Count occurrences of each ID
id_counts <- ESWD_LH_2000_2024 %>% count(ID)
# Filter IDs that appear more than once
duplicate_ids <- id_counts %>% filter(n > 1) %>% pull(ID)
# Extract rows with duplicate IDs
duplicate_rows <- ESWD_LH_2000_2024 %>%
  filter(ID %in% duplicate_ids)
# If you look at this sorting by ID, it is clear that events with the same ID are years apart
# ID cannot be used to characterize different events, variable to ignore.

range(as.Date(ESWD_LH_2000_2024$TIME_EVENT)) #"2000-01-13" "2024-07-06"

# Impact "I1" refers explicitly to energy infrastructure, however only 6 reports for hailes
# have this code.
length(which(ESWD_LH_2000_2024$IMPACTS=="I1")) #only 2

# 14196 of 53670 total reports have NA max diameter.
summary(ESWD_LH_2000_2024$MAX_HAIL_DIAMETER)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   2.000   3.000   3.457   4.000  19.000   14196 
# There are a total of 8 reports of diameters smaller than 2 cm
length(which(ESWD_LH_2000_2024$MAX_HAIL_DIAMETER>0)) 
length(which(ESWD_LH_2000_2024$MAX_HAIL_DIAMETER>=8)) 
hist(ESWD_LH_2000_2024$MAX_HAIL_DIAMETER[which(ESWD_LH_2000_2024$MAX_HAIL_DIAMETER<1)]) 
min(ESWD_LH_2000_2024$MAX_HAIL_DIAMETER[ESWD_LH_2000_2024$MAX_HAIL_DIAMETER>0], na.rm = T) # 0.8

## there are data in all of Russia, including Siberia. Also, some reports have unrealistic latitude.
## Remove anything below 30 degrees N and further East than 50E
##ind_removed <- c(which(ESWD_LH_2000_2024$LATITUDE<30), which(ESWD_LH_2000_2024$LATITUDE>72), 
##                 which(ESWD_LH_2000_2024$LONGITUDE>50), which(ESWD_LH_2000_2024$LONGITUDE<(-25.5)))

# Remove anything outside the ERA5 domain
ind_removed <- c(which(ESWD_LH_2000_2024$LATITUDE<35), which(ESWD_LH_2000_2024$LATITUDE>60), 
                 which(ESWD_LH_2000_2024$LONGITUDE>25), which(ESWD_LH_2000_2024$LONGITUDE<(-15)))
hail_data <- ESWD_LH_2000_2024[-ind_removed,]


## DAILY, MONTHLY AND YEARLY NUMBER OF EVENTS

# Aggregate data into yearly totals
hail_data <- hail_data %>%
  mutate(date = as.Date(TIME_EVENT),
         year = format(date, "%Y")) 
annual_reports_hail <- hail_data %>%
  group_by(year) %>%
  summarise(total_reports = n())

# Aggregate data into monthly totals
monthly_reports_hail <- hail_data %>%
  mutate(year_month = floor_date(date, "month")) %>% # Create a year-month column
  group_by(year_month) %>% # Group by year-month
  summarise(total_reports = n(), .groups = "drop") # Sum reports for each month

# Aggregate data into daily totals
daily_reports_hail <- hail_data %>%
  group_by(date) %>%
  summarise(total_reports = n())

# annual reports by country
annual_reports_hail_by_country <- hail_data %>%
  group_by(year, COUNTRY) %>%
  summarise(reports = n(), .groups = 'drop')


## SEASONALITY

hail_data$MONTH <- month(hail_data$TIME_EVENT)

# Aggregate daily data into monthly totals
seasonality_reports_hail <- hail_data %>%
  group_by(MONTH) %>% # Group by month
  summarise(total_reports = n())
# Convert month numbers to month names
seasonality_reports_hail <- seasonality_reports_hail %>%
  mutate(month_name = factor(MONTH, levels = 1:12, labels = month.name))

# seasonal cycle by country
seasonality_hail_by_country <- hail_data %>%
  group_by(MONTH, COUNTRY) %>%
  summarise(reports = n(), .groups = 'drop')
# Convert month numbers to month names
seasonality_hail_by_country <- seasonality_hail_by_country %>%
  mutate(month_name = factor(MONTH, levels = 1:12, labels = month.name))


##### GRIDDED REPORTS
# using values from the database
#lon <- seq(from = min(floor(hail_data$LONGITUDE)), to = ceiling(max(hail_data$LONGITUDE)), by = 0.25)
#lat <- seq(from = min(floor(hail_data$LATITUDE)), to = max(ceiling(hail_data$LATITUDE)), by = 0.25)

# using the ERA5 grid
lon <- seq(from = -15, to = 25, by = 0.25)
lat <- seq(from = 35, to = 60, by = 0.25)
europe_bbox <- c(xmin = min(lon), xmax = max(lon), ymin = min(lat), ymax = max(lat))  # Longitude and latitude

# Map hail data onto the grid
grid <- expand.grid(lon = lon, lat = lat)
hail_sf <- st_as_sf(hail_data, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
nearest_grid_index <- st_nearest_feature(hail_sf, grid_sf)
hail_data$grid_index <- nearest_grid_index

# Calculate hail statistics
total_hail_counts <- hail_data %>%
  group_by(grid_index) %>%
  summarise(count = n())

grid_with_counts <- grid_sf %>%
  mutate(grid_index = row_number()) %>%
  left_join(total_hail_counts, by = "grid_index") %>%
  mutate(lon = st_coordinates(.)[, 1],  # Extract longitude
         lat = st_coordinates(.)[, 2])  # Extract latitude

#Create bins for counts
grid_with_counts <- grid_with_counts %>%
  mutate(count_bin = cut(count,
                         breaks = c(0, 4, 9, 29, 49, 74, 99, 149, 299, Inf),
                         labels = c("1-4", "5-9", "10-29", "30-49", "50-74", "75-99", "100-149", "150-299", "300+")))

#######
# Get gridded daily reports, store in array and save to netcdf

# Get unique dates
unique_dates <- sort(unique(hail_data$date))
TT <- length(unique_dates)

# Initialize 3D array (lon × lat × time)
hail_array <- array(0, dim = c(length(lon), length(lat), TT),
                    dimnames = list(lon, lat, as.character(unique_dates)))

# Fill the array
for (i in 1:TT) {
  print(i)
  current_date <- unique_dates[i]
  daily_data <- hail_data %>% filter(date == current_date)
  
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
  hail_array[,,i] <- matrix(grid_counts, nrow = length(lon), ncol = length(lat))
}

# Save as NetCDF file
save_as_netcdf <- function(array_3d, lon, lat, dates, filename) {
  # Define dimensions
  londim <- ncdim_def("longitude", "degrees_east", lon)
  latdim <- ncdim_def("latitude", "degrees_north", lat)
  timedim <- ncdim_def("time", "days since 1970-01-01", 
                       as.numeric(dates - as.Date("1970-01-01")))
  
  # Define variable
  hail_def <- ncvar_def("hail_reports", "count", 
                        list(londim, latdim, timedim), 
                        missval = -999, 
                        longname = "Daily hail report counts")
  
  # Create NetCDF file
  ncout <- nc_create(filename, list(hail_def))
  
  # Put the array
  ncvar_put(ncout, hail_def, array_3d)
  
  # Add attributes
  ncatt_put(ncout, "longitude", "axis", "X")
  ncatt_put(ncout, "latitude", "axis", "Y")
  ncatt_put(ncout, "time", "axis", "T")
  ncatt_put(ncout, 0, "title", "Daily hail reports gridded data")
  ncatt_put(ncout, 0, "source", "ESWD database")
  ncatt_put(ncout, 0, "grid_resolution", "0.25 degrees")
  
  nc_close(ncout)
}

# Save the data
save_as_netcdf(hail_array, lon, lat, unique_dates, "hail_reports_gridded_ESWD.nc")

###############################
###### PLOTS ##################
###############################

# plot time series of yearly report
ggplot(annual_reports_hail, aes(x = year, y = total_reports, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Yearly number of hail reports",
       x = "Year",
       y = "Report Number") +
  theme_minimal() +
  scale_x_discrete(breaks = seq(min(annual_reports_hail$year), max(annual_reports_hail$year), by = 5)) + # Show every 5th year
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) # Rotate and reduce font size

# plot time series of daily reports
ggplot(daily_reports_hail, aes(x = date, y = total_reports)) +
  geom_line(color = "blue") +
  geom_point(color = "red", size = 1) +
  labs(title = "Daily number of hail reports",
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

# plot time series of monthly report
ggplot(monthly_reports_hail, aes(x = year_month, y = total_reports, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Monthly number of hail reports",
       x = "Year",
       y = "Report Number") +
  theme_minimal() +
  scale_x_date(breaks = seq(min(monthly_reports_hail$year_month), max(monthly_reports_hail$year_month), by = "5 years"), # Show every year
               date_labels = "%Y") + # Format labels as years
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) # Rotate and reduce font size

# plot annual reports by country
ggplot(annual_reports_hail_by_country, aes(x = year, y = COUNTRY, fill = reports)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +  # Usa una scala di colori
  labs(title = "Annual number of hail reports by country",
       x = "Year",
       y = "Country",
       fill = "Report number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot seasonal cycle
ggplot(seasonality_reports_hail, aes(x = month_name, y = total_reports)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = total_reports), vjust = -0.5, color = "black") + # Add labels
  labs(title = "Seasonal Cycle of Hail Reports",
       x = "Month",
       y = "Total Reports") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot seasonal cycle by country
ggplot(seasonality_hail_by_country, aes(x = month_name, y = COUNTRY, fill = reports)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1) +  # Usa una scala di colori
  labs(title = "Annual number of hail reports by country",
       x = "Month",
       y = "Country",
       fill = "Report number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Plot gridded total reports
world_map <- ne_countries(scale = "medium", returnclass = "sf")
region_map <- st_crop(world_map, europe_bbox)

ggplot() +
  geom_tile(data = grid_with_counts, aes(x = lon, y = lat, fill = count_bin), color = NA) +
  geom_raster() +
  scale_fill_viridis_d(option = "plasma", na.value = "white", direction = -1,
                       labels = c("1-4", "5-9", "10-29", "30-49", "50-74", "75-99", "100-149", "150-299", "300+"))+
  borders("world", colour = "black") +
  coord_sf(xlim = c(europe_bbox["xmin"], europe_bbox["xmax"]),
           ylim = c(europe_bbox["ymin"], europe_bbox["ymax"])) +
  labs(title = "Hail Reports on Regular 0.25-Degree Grid (Europe)",
       x = "Longitude",
       y = "Latitude",
       fill = "Hail Reports Count") +
  theme_minimal()

