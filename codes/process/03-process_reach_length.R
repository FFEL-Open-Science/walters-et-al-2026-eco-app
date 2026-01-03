# Approach ----

# - Download data
# - Read in with sf and dissolve
# - Convert to sp object
# - Convert sp object to network using line2network
# - Cleanup network using cleanup function
# - Read in reach coordinates
# - Snap coordinates to river network using xy2segvert
# - Compute distance between all pairs of points using riverdistancematbysurvey()


# Initial setup ----
library(bcdata)
library(tidyverse)
library(sf)
library(riverdist)


# Process spatial data ----

## Read in and process data ----

# Get metadata from BC Data Catalogue on the Freshwater Atlas Strea Network
stream_network_metadata <- bcdc_get_record("92344413-8035-4c08-b996-65a9b3f62fca")

# Download dataset and filter for the Stellako River
stellako <- bcdc_query_geodata(stream_network_metadata) %>%
  filter(GNIS_NAME == "Stellako River" | 
           OBJECTID %in% c(13047527, 13063397, 13074723, 13045631, 13057349)) %>% 
  collect() %>%
  # dissolves segments into one line feature (required by riverdist)
  summarize(geometry = st_union(geometry)) %>% 
  st_zm()

# Process data and compute distances with riverdist functions ----

# The network only needs to be re-created and cleaned if changes to the input 
# shapefile and cleaning parameters are needed.

#Answers to questions by cleanup(): y, y, y, 10, 5, 1476, y, n, n, n
# stellako_network <- cleanup(line2network(stellako))
# 
# saveRDS(stellako_network, "data/gis/stellako_network_riverdist.rds")

stellako_network <- readRDS("data/gis/stellako_network_riverdist.rds")

reaches <- st_read("data/gis/stellako_reaches.kml") %>%
  st_transform(crs = 3005) %>%
  mutate(rkm = NA,
         x = st_coordinates(.)[, 1],
         y = st_coordinates(.)[, 2])

reaches_snapped <- xy2segvert(x = reaches$x, y = reaches$y, 
                              rivers = stellako_network) %>%
  mutate(reach_site = reaches$Name,
         id = 1) # fake id for riverdistancematbysurvey function

# Modify code below to check snapping
plot(stellako_network)
zoomtoseg(seg=c(4, 3), rivers = stellako_network)
points(reaches$x, reaches$y, pch = 16, col = "red")
riverpoints(reaches_snapped$seg, reaches_snapped$vert, rivers = stellako_network, 
            pch = 15, col = "blue")

# Compute distances between points
reaches_distance <- riverdistancematbysurvey(indiv = 1, unique = 1, 
                                             survey = reaches_snapped$reach_site,
                                             seg = reaches_snapped$seg, 
                                             vert = reaches_snapped$vert,
                                             rivers = stellako_network)


# Get coordinates of snapped points
fnw_lines <- stellako_network$lines
lon_snapped <- rep(NA, nrow(reaches_snapped))
lat_snapped <- rep(NA, nrow(reaches_snapped))

for(i in 1:nrow(reaches_snapped)) {
  lon_snapped[i] <- fnw_lines[[reaches_snapped$seg[i]]][reaches_snapped$vert[i], 1]
  lat_snapped[i] <- fnw_lines[[reaches_snapped$seg[i]]][reaches_snapped$vert[i], 2]
}

reaches <-  st_drop_geometry(reaches) %>%
  arrange(Name) %>%
  rename(site = Name) %>%
  mutate(rm = round(reaches_distance[1, ], 2),
  lon = lon_snapped,
  lat = lat_snapped) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 3005) %>%
  st_transform(crs = 4326) %>%
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2],
         length = lead(rm) - rm,
         reach = c(1:8, NA)) %>%
  st_drop_geometry() %>%
  select(site, rm, lon, lat, reach, length)
  
saveRDS(reaches, "data/processed/reach_lengths.rds")
