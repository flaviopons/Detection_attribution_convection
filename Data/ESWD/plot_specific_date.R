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
length(which(ESWD_LH_2000_2024$MAX_HAIL_DIAMETER>0)) # 0.8
hist(ESWD_LH_2000_2024$MAX_HAIL_DIAMETER[which(ESWD_LH_2000_2024$MAX_HAIL_DIAMETER<1)]) # 0.8
min(ESWD_LH_2000_2024$MAX_HAIL_DIAMETER[ESWD_LH_2000_2024$MAX_HAIL_DIAMETER>0], na.rm = T) # 0.8

# there are data in all of Russia, including Siberia. Also, some reports have unrealistic latitude.
# Remove anything below 30 degrees N and further East than 50E
ind_removed <- c(which(ESWD_LH_2000_2024$LATITUDE<30), which(ESWD_LH_2000_2024$LATITUDE>72), 
                 which(ESWD_LH_2000_2024$LONGITUDE>50), which(ESWD_LH_2000_2024$LONGITUDE<(-25.5)))
hail_data <- ESWD_LH_2000_2024[-ind_removed,]

# Create date variable
hail_data <- hail_data %>%
  mutate(date = as.Date(TIME_EVENT),
         year = format(date, "%Y")) 

# using the ERA5 grid
lon <- seq(from = 0, to = 18, by = 0.25)
lat <- seq(from = 40, to = 50, by = 0.25)
europe_bbox <- c(xmin = min(lon), xmax = max(lon), ymin = min(lat), ymax = max(lat))  # Longitude and latitude

# Map hail data onto the grid
grid <- expand.grid(lon = lon, lat = lat)
hail_sf <- st_as_sf(hail_data, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
nearest_grid_index <- st_nearest_feature(hail_sf, grid_sf)
hail_data$grid_index <- nearest_grid_index

# Filter date
hail_24july <- hail_data %>% 
  filter(date == as.Date("2023-07-24"))

# convert to sf object
hail_sf_24july <- st_as_sf(hail_24july, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

# create grid
grid_sf <- st_as_sf(expand.grid(lon = lon, lat = lat), 
                    coords = c("lon", "lat"), 
                    crs = 4326)

highlight_region <- data.frame(
  lon = c(9, 14, 14, 9, 9),  # Longitudine (x)
  lat = c(45.2, 45.2, 46.2, 46.2, 45.2) # Latitudine (y)
)

ggplot() +
  # Disegna la griglia
  geom_sf(data = grid_sf, shape = 3, color = "gray80", size = 0.1, alpha = 0.5) +
  
  # Aggiungi confini Europa
  geom_sf(data = ne_countries(scale = "medium", returnclass = "sf") %>% 
            st_crop(europe_bbox), 
          fill = "lightgray", color = "black") +
  
  
  # Aggiungi eventi grandine come rombi rossi
  geom_sf(data = hail_sf_24july, 
          shape = 18, color = "red", size = 1.5, alpha = 0.8) +
  
  # Rettangolo evidenziato (VERDE)
  geom_polygon(data = highlight_region, 
               aes(x = lon, y = lat), 
               color = "green", 
               fill = NA, 
               linewidth = 1.2) +
  
  # Zoom sull'area di interesse
  coord_sf(xlim = c(0, 18), ylim = c(40, 50)) +
  
  # Stile e titolo
  labs(title = "ESSL hail reports, July 24, 2023",
       caption = paste("Total reports:", nrow(hail_24july))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "azure"))


ggsave("hail_events_24july2023_ERA5grid.png", 
       width = 10, height = 8, dpi = 300, bg = "white")
