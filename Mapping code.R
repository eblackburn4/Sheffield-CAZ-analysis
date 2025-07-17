
## ---------------------------
## Purpose of script: code for generation of maps in dissertation
## Author: Ned Blackburn
## Date Created: 2025-03-19

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)
library(ggthemes)
library(ggmap)
library(jsonlite)
library(sf)


# map of sensor locations -------------------------------------------------
# this code generates a map of the CAZ boundary and sensor locations in Sheffield
#make base map

sheffield <- c(-1.550367, 53.310425, -1.390238, 53.423676)
basemap <- get_stadiamap(sheffield, zoom = 15, maptype="stamen_toner_lite") |> ggmap()

#read in CAZ json

polygon_data <- fromJSON("mapfinal.geojson")

# Extract coordinates

coords <- polygon_data$features$geometry$coordinates[[1]]
polygon_coords <- coords[1, , ]

# Convert the resulting matrix into a data frame.
df_polygon <- as.data.frame(polygon_coords)
colnames(df_polygon) <- c("lon", "lat")

#load in sensor metadata (SCC etc)

SCCsensormeta <- read_csv("Data/SCC/SCC_sensor_meta.csv")|>
  select(sensor_id = Location, family. = Owner, lat = Lat, lon = Lon) %>%
  mutate(type = "Both") 

json_sensor_meta <- list.files(path = "Data/Sensormeta", pattern = "\\.json$", full.names = TRUE)
  

#Function to extract location data from JSONs
read_sensor_meta <- function(fn) {
  j <- fromJSON(fn)
  has_no2  <- any(grepl("\\.NO2$",  names(j$columns)))
  has_pm25 <- any(grepl("\\.PM25$", names(j$columns)))
  sensor_type <- if (has_no2 && has_pm25) {
    "Both"
  } else if (has_no2) {
    "NO2"
  } else if (has_pm25) {
    "PM2.5"
  } else {
    NA_character_
  }
  data.frame(
    sensor_id   = j$identity$siteID,
    family.     = j$identity$sensorFam,
    lat         = j$location$latitude_deg,
    lon         = j$location$longitude_deg,
    type        = sensor_type,
    stringsAsFactors = FALSE
  )
}

df_sensor_meta <- do.call(rbind, lapply(json_sensor_meta, read_sensor_meta))
all_sensor_meta <- rbind(SCC_sensor_meta, df_sensor_meta)

basemap + 
  geom_polygon(data = df_polygon, aes(x = lon, y = lat),
               fill = 'pink',        
               color = "red",
               alpha = 0.8,   
               size = 1) +
  geom_point(data = all_sensor_meta, aes(x = lon, y = lat, color = type), 
             size = 1) +
  scale_color_manual(values = c("Both" = "#35b779", "NO2" = "#ffd500", "PM2.5" = "#31688e")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())



# build CAZ adjacent polygon ----------------------------------------------
# This code builds a 500 m buffer around the CAZ boundary line and identifies sensors within that buffer.

#create sf objects with sensors and CAZ boundary
sensor_sf <- st_as_sf(all_sensor_meta, coords = c("lon", "lat"), crs = 4326)

caz_polygon <- st_polygon(list(as.matrix(df_polygon[, c("lon", "lat")])))
caz_sf      <- st_sfc(caz_polygon, crs = 4326) %>% st_sf()

#Project both to a metric coordinate system for metre‐based distances
sensor_m <- st_transform(sensor_sf, 27700)
caz_m    <- st_transform(caz_sf, 27700)

#Build a 500 m exterior “ring” around the CAZ
caz_buffer   <- st_buffer(caz_m, 500)
caz_outer_ring <- st_difference(caz_buffer, caz_m)

#Classify each sensor into one of three bins
sensor_m <- sensor_m %>%
  mutate(
    in_caz      = as.logical(st_within(geometry, caz_m, sparse = FALSE)),
    in_ring     = as.logical(st_within(geometry, caz_outer_ring, sparse = FALSE)),
    category    = case_when(
      in_caz                    ~ "Inside CAZ",
      !in_caz & in_ring         ~ "CAZ Adjacent",
      TRUE                      ~ "Other"
    )
  )

# 2) reproject your CAZ polygon, outer‐ring and sensors back to WGS84
caz_ll      <- st_transform(caz_m,            4326)
ring_ll     <- st_transform(caz_outer_ring,   4326)
sensor_ll   <- st_transform(sensor_m,         4326)

# 3) extract coords + grouping for plotting polygons with geom_polygon()
sf_to_df <- function(sf_poly) {
  coords <- st_coordinates(sf_poly)
  data.frame(
    lon   = coords[,"X"],
    lat   = coords[,"Y"],
    group = coords[,"L1"]
  )
}
df_caz     <- sf_to_df(caz_ll)
df_ring    <- sf_to_df(ring_ll)

# 4) pull lon/lat + category out of sensor_ll
sensor_df_ll <- sensor_ll %>%
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

# 5) plot everything on the ggmap
basemap +
  # the 500 m exterior ring
  geom_polygon(
    data    = df_ring,
    aes(x    = lon, y = lat, group = group),
    fill  = alpha("lightblue", 0.6),
    color   = "blue",
    size    = 0.8,
  ) +
  # the CAZ polygon
  geom_polygon(
    data    = df_caz,
    aes(x    = lon, y = lat, group = group),
    fill  = alpha("pink", 0.5),
    color   = "red",
    size    = 0.8
  ) +
  geom_point(
    data = sensor_df_ll,
    aes(x   = lon, y = lat, shape = category),
    size = 1
  ) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank())
