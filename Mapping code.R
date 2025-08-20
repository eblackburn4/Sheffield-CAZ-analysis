
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
library(ggrepel)
library(nngeo)
library(osmdata)


# map of sensor locations -------------------------------------------------
# this code generates a map of the CAZ boundary and sensor locations in Sheffield
#make base map

sheffield <- c(-1.524367, 53.350425, -1.404986, 53.406676)
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

SCC_sensor_meta <- read_csv("Data/Sensormeta/SCC_sensor_meta.csv")|>
  select(sensor_id = Location, family = Owner, lat = Lat, lon = Lon) |>
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
    sensor_id   = j$identity$Q_sensor_ID,
    family     = j$identity$sensorFam,
    lat         = j$location$latitude_deg,
    lon         = j$location$longitude_deg,
    type        = sensor_type,
    stringsAsFactors = FALSE
  )
}

df_sensor_meta <- do.call(rbind, lapply(json_sensor_meta, read_sensor_meta))
all_sensor_meta <- rbind(SCC_sensor_meta, df_sensor_meta)

#load in traffic sensor metadata
traffic_sensor_meta <- read_csv("Data/Sensormeta/traffic_meta.csv") |>
  select(sensorID, lon = 'long_[deg]', lat = 'lat_[deg]')

# build CAZ adjacent polygon map ----------------------------------------------
# This code builds a 1000 m buffer around the CAZ boundary line and identifies sensors within that buffer.

#create sf objects with AQ/traffic sensors and CAZ boundary
aq_sensor_sf <- st_as_sf(all_sensor_meta, coords = c("lon", "lat"), crs = 4326)
tf_sensor_sf <- st_as_sf(traffic_sensor_meta, coords = c("lon", "lat"), crs = 4326)

caz_polygon <- st_polygon(list(as.matrix(df_polygon[, c("lon", "lat")])))
caz_sf      <- st_sfc(caz_polygon, crs = 4326) |> st_sf()

#Project both to a metric coordinate system for metre‐based distances
aq_sensor_m <- st_transform(aq_sensor_sf, 27700)
tf_sensor_m <- st_transform(tf_sensor_sf, 27700)

caz_m    <- st_transform(caz_sf, 27700)

#Build an exterior “ring” around the CAZ
caz_buffer   <- st_buffer(caz_m, 1500)
caz_outer_ring <- st_difference(caz_buffer, caz_m)

#Classify each aq sensor into one of three bins
aq_sensor_m <- aq_sensor_m |>
  mutate(
    in_caz      = as.logical(st_within(geometry, caz_m, sparse = FALSE)),
    in_ring     = as.logical(st_within(geometry, caz_outer_ring, sparse = FALSE)),
    category    = case_when(
      in_caz                    ~ "Inside CAZ",
      !in_caz & in_ring         ~ "CAZ Adjacent",
      TRUE                      ~ "Other"
    )
  )

#Classify each traffic sensor into one of three bins
tf_sensor_m <- tf_sensor_m |>
  mutate(
    in_caz      = as.logical(st_within(geometry, caz_m, sparse = FALSE)),
    in_ring     = as.logical(st_within(geometry, caz_outer_ring, sparse = FALSE)),
    category    = case_when(
      in_caz                    ~ "Inside CAZ",
      !in_caz & in_ring         ~ "CAZ Adjacent",
      TRUE                      ~ "Other"
    )
  )


#reproject back to normal coords
caz_ll      <- st_transform(caz_m, 4326)
ring_ll     <- st_transform(caz_outer_ring, 4326)
aq_sensor_ll   <- st_transform(aq_sensor_m, 4326)
tf_sensor_ll   <- st_transform(tf_sensor_m, 4326)

#extract coords 
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

# 4) pull lon/lat + category out of aq_sensor_ll/tf_sensor_ll
aq_sensor_df_ll <- aq_sensor_ll %>%
  mutate(coords = st_coordinates(geometry),
         lon = coords[,1],
         lat = coords[,2]) %>%
  st_drop_geometry()

tf_sensor_df_ll <- tf_sensor_ll %>%
  mutate(coords = st_coordinates(geometry),
         lon = coords[,1],
         lat = coords[,2],
         sensorID = str_replace_all(sensorID, "[^[:alnum:]]", "")) %>%
  st_drop_geometry()

# plot aq on basemap
basemap +
  # the exterior ring
  geom_polygon(
    data    = df_ring,
    aes(x    = lon, y = lat, group = group),
    fill  = alpha("lightblue", 0.6),
    color   = "blue",
    linewidth    = 0.6,
  ) +
  # the CAZ polygon
  geom_polygon(
    data    = df_caz,
    aes(x    = lon, y = lat, group = group),
    fill  = alpha("pink", 0.5),
    color   = "red",
    linewidth   = 0.6
  ) +
  geom_point(
    data = aq_sensor_df_ll,
    aes(x   = lon, y = lat, shape = category),
    size = 3, alpha = 0.8
  ) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) 

#plot tf on basemap
basemap +
  # the exterior ring
  geom_polygon(
    data    = df_ring,
    aes(x    = lon, y = lat, group = group),
    color   = "black",
    fill = alpha("darkblue", 0.1),
    size    = 0.3,
  ) +
  # the CAZ polygon
  geom_polygon(
    data    = df_caz,
    aes(x    = lon, y = lat, group = group),
    color   = "#cc4778",
    fill   = alpha("#cc4778", 0.2),
    linewidth    = 0.3
  ) +
  geom_point(
    data = tf_sensor_df_ll,
    aes(x   = lon, y = lat, color = category),
    size = 2, alpha = 0.8
  ) +
  scale_color_manual(values = c("Inside CAZ" = "#9c335a", 
                                  "CAZ Adjacent" = "darkblue", 
                                  "Other" = "grey50")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank())


# traffic sensor road matching --------------------------------------------
#this section uses OSM data to match sensors to major roads around the CAZ
 
# Inputs assumed:
# - caz_m           : CAZ polygon in EPSG:27700
# - caz_outer_ring  : st_difference(st_buffer(caz_m, 1500), caz_m)  # 0–1500 m outside
# - sensors         : data.frame(sensorID, lon, lat, ...)

# 1) Get OSM roads around the CAZ (bbox of the outer ring) and within the CAZ, keep relevant classes, project to 27700

bb <- st_bbox(st_transform(caz_outer_ring, 4326))
od <- opq(bbox = bb) |>
  add_osm_feature(key = "highway") |>
  osmdata_sf()

roads_raw <- od$osm_lines |>
  st_transform(27700) |>
  filter(highway %in% c("motorway","trunk","primary","secondary",
                        "tertiary","unclassified","residential","living_street")) |>
  st_make_valid()

# Roads intersecting the inner CAZ
roads_inner <- st_intersection(roads_raw, caz_m) |>
  filter(!st_is_empty(geometry))

# Roads in the CAZ-adjacent ring
roads_outer <- st_intersection(roads_raw, caz_outer_ring) |>
  filter(!st_is_empty(geometry))


tf_sensor_m <- tf_sensor_m |>
  mutate(side = case_when(in_caz ~ "inner", in_ring ~ "outer", TRUE ~ NA_character_)) |>
  filter(!is.na(side))

# Helper: match sensors to roads (buffer membership first, then nearest within a threshold)
match_to_roads <- function(pts, roads, buf_m = 15, nearest_max_m = 30) {
  if (nrow(pts) == 0) {
    return(pts |> st_drop_geometry() |>
             mutate(road_id = NA, name = NA, ref = NA, highway = NA,
                    match_method = NA_character_, match_dist_m = NA_real_))
  }
  roads <- st_make_valid(roads)
  roads_buf <- st_buffer(roads, buf_m)
  
  # 1) Buffer membership join (may yield multiple candidates per sensor)
  cand <- st_join(
    pts |> mutate(.sid = row_number()),
    roads_buf[, c("osm_id","name","ref","highway")],
    left = TRUE
  )
  
  # Keep the closest centreline among any buffer matches
  has_match <- !is.na(cand$osm_id)
  if (any(has_match)) {
    idx <- match(cand$osm_id[has_match], roads$osm_id)
    d_to_line <- st_distance(st_geometry(cand)[has_match], st_geometry(roads)[idx], by_element = TRUE)
    cand$dist_m <- NA_real_
    cand$dist_m[has_match] <- as.numeric(d_to_line)
    cand <- cand |>
      group_by(.sid) |>
      slice_min(dist_m, n = 1, with_ties = FALSE) |>
      ungroup()
  } else {
    cand <- cand |> mutate(dist_m = NA_real_)
  }
  
  # 2) Fallback: nearest road within threshold for any still-unmatched
  i_na <- is.na(cand$osm_id)
  if (any(i_na)) {
    idx_nn <- st_nearest_feature(cand[i_na, ], roads)
    d_nn   <- st_distance(cand[i_na, ], roads[idx_nn, ], by_element = TRUE)
    ok     <- as.numeric(d_nn) <= nearest_max_m
    
    cand$osm_id[i_na][ok]  <- roads$osm_id[idx_nn[ok]]
    cand$name[i_na][ok]    <- roads$name[idx_nn[ok]]
    cand$ref[i_na][ok]     <- roads$ref[idx_nn[ok]]
    cand$highway[i_na][ok] <- roads$highway[idx_nn[ok]]
    cand$dist_m[i_na][ok]  <- as.numeric(d_nn[ok])
  }
  
  cand |>
    st_drop_geometry() |>
    transmute(
      sensorID, category,
      road_id      = osm_id,
      name, ref, highway,
      match_method = case_when(
        is.na(road_id)      ~ NA_character_,
        dist_m <= buf_m     ~ "buffer",
        TRUE                ~ "nearest"
      ),
      match_dist_m = dist_m
    )
}

# Apply per side and combine
matched_inner <- tf_sensor_m |> filter(side == "inner") |>
  match_to_roads(roads_inner, buf_m = 15, nearest_max_m = 30)

matched_outer <- tf_sensor_m |> filter(side == "outer") |>
  match_to_roads(roads_outer, buf_m = 15, nearest_max_m = 30)

sensor_to_road <- bind_rows(matched_inner, matched_outer)

# set all ring road sensors to be inside CAZ and as a single variable
ring_road <- c("Saint Mary's Road",
               "Saint Mary's Gate",
               "Suffolk Road",
               "Sheaf Street",
               "Sheffield Parkway",
               "Derek Dooley Way",
               "Shalesmoor",
               "Netherthorpe Road",
               "Upper Hanover Street",
               "Hanover Way",
               "Hoyle Street")
            
#wrangle data to get the categories right

sensor_to_road_RDD <- sensor_to_road |>
  mutate(category = case_when(
    name %in% ring_road ~ "Inside CAZ",
    TRUE                ~ category
  )) |>
  mutate(name = case_when(
    name %in% ring_road ~ "Ring Road",
    TRUE                ~ name
  )) |>
  mutate(ref = case_when(
    ref == "A61" & category == "Inside CAZ" ~ "A61 Ring Road",
    ref == "A61" & category == "CAZ Adjacent" ~ "A61",
    TRUE ~ ref
  )) |>
  mutate(ref = case_when(
    category == 'Inside CAZ' & ref != 'A61 Ring Road' ~ 'CAZ interior road',
    TRUE ~ ref
  )) |>
  mutate(highway = ifelse(highway == 'trunk', 'primary', highway)) |>
  filter(highway == 'primary' | highway == 'secondary') |>
  mutate(sensorID = str_replace_all(sensorID, "[^A-Za-z0-9]", "")) |>
  filter(ref !=is.na(ref))

#plot the roads on a map

roads_map <- roads_outer |>
  filter(!is.na(ref)) |>
  select(osm_id, name, ref, highway, geometry) |>
  st_transform(4326) |>
  mutate(high_cat = case_when(
    highway == "primary"   ~ "primary",
    highway == "secondary" ~ "secondary",
    highway == "trunk" ~ "primary",
    TRUE                   ~ "other"
  ))




bbox_ll <- st_as_sfc(st_bbox(c(xmin = sheffield[1], ymin = sheffield[2],
                               xmax = sheffield[3], ymax = sheffield[4]), crs = 4326))
roads_map <- st_crop(roads_map, bbox_ll)
caz_ll       <- st_intersection(caz_ll, bbox_ll)

# Plot over your ggmap basemap (basemap object already created)
basemap +
  geom_sf(data = roads_map,
          aes(color = high_cat),
          inherit.aes = FALSE, linewidth = 0.9) +
  coord_sf(xlim = c(sheffield[1], sheffield[3]),
           ylim = c(sheffield[2], sheffield[4]),
           expand = FALSE, datum = NA) +
  geom_polygon(
    data    = df_caz,
    aes(x    = lon, y = lat, group = group),
    fill  = alpha("pink", 0.5),
    color   = "red",
    linewidth    = 1
  ) +
  scale_color_manual(values = c("primary" = "darkblue", 
                                "secondary" = "limegreen")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) 
                            