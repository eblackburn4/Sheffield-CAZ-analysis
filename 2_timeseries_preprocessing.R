
## ---------------------------
## Purpose of script: data exploration and preparation code for dissertation
## Author: Ned Blackburn
## Date Created: 2025-06-27

## ---------------------------

# read in and QA air quality data from each sensor family --------------------------------------------------------

SCC_sensor_names <- list.files("Data/SCC", pattern = "\\.csv$", full.names = TRUE)
names(SCC_sensor_names) <- basename(SCC_sensor_names)

EW_sensor_names <- list.files("Data/eWatch", pattern = "\\.csv$", full.names = TRUE)
names(EW_sensor_names) <- basename(EW_sensor_names)

defra_sensor_names <- list.files("Data/defra", pattern = "\\.csv$", full.names = TRUE)
names(defra_sensor_names) <- basename(defra_sensor_names)

AM_sensor_names <- list.files("Data/AMfixed", pattern = "\\.csv$", full.names = TRUE)
names(AM_sensor_names) <- basename(AM_sensor_names)


#load in all csvs into a single dataframe for each sensor family and clean up/reformat

# SCC data
SCC_aq_raw <- map_dfr(SCC_sensor_names, read_csv, .id = "file_id") |>
  pivot_longer(
    cols         = matches("(NO2|PM25)"),
    names_to     = c("SensorID","Pollutant"),
    names_pattern= "([^,]+),\\s*(NO2|PM25),.*",
    values_to    = "Value"
  ) |>
  group_by(SensorID, DateTime, Pollutant) |>
  summarize(Value = first(na.omit(Value)), .groups = "drop") |>
  pivot_wider(
    names_from   = Pollutant,
    values_from  = Value
  ) |>
  mutate(SensorID = case_when(
    SensorID == 'Firvale' ~ 'SCC_GH1',
    SensorID == 'Lowfield' ~ 'SCC_GH3',
    SensorID == 'The Wicker' ~ 'SCC_GH4',
    SensorID == 'King Ecgbert' ~ 'SCC_GH5',
    SensorID == 'Pond Hill' ~ 'SCC_GH6',
  )) |>
  select(SensorID, DateTime, NO2, PM25)


# eWatch data
EW_aq_raw <- map_dfr(EW_sensor_names, read_csv, .id = "file_id") |>
  select(data.sensor,
         data.universalTime, 
         data.NO2, 
         'data.NO2:QC') |>
  rename(SensorID = data.sensor,
         DateTime = data.universalTime,
         NO2      = data.NO2,
         NO2_QC = `data.NO2:QC`)
         
# Defra data
defra_aq_raw <- map_dfr(defra_sensor_names, read_csv, .id = "file_id") |>
  select(data.sensor,
         data.time, 
         data.universalTime, 
         data.PM25, 
         'data.PM25:QC', 
         data.NO2, 
         'data.NO2:QC') |>
  rename(SensorID = data.sensor,
         DateTime = data.universalTime,
         NO2      = data.NO2,
         PM25     = data.PM25,
         NO2_QC = `data.NO2:QC`,
         PM25_QC = `data.PM25:QC`)


#AQMesh data
AM_aq_raw <- map_dfr(AM_sensor_names, read_csv, .id = "file_id") |>
  select(amon.sensor,
         amon.universalTime, 
         amon.PM25, 
         'amon.PM25:QC', 
         amon.NO2, 
         'amon.NO2:QC') |>
  rename(
    SensorID = amon.sensor,
    DateTime = amon.universalTime,
    NO2      = amon.NO2,
    PM25     = amon.PM25,
    NO2_QC = `amon.NO2:QC`,
    PM25_QC = `amon.PM25:QC`)
  
#combine into one master df
master_aq_raw <- bind_rows(
  SCC_aq_raw |> mutate(family = "SCC"),
  EW_aq_raw   |> mutate(family = "eWatch"),
  defra_aq_raw |> mutate(family = "Defra"),
  AM_aq_raw   |> mutate(family = "AQMesh")
)

#remove individual datasets
rm(SCC_aq_raw, EW_aq_raw, defra_aq_raw, AM_aq_raw)


#join to sensor metadata and filter for only CAZ-adjacent sensors with at least 6 months of data
#also filter out 1-5am timestamps for SCC_GH4_NO2

start <- as.POSIXct("2022-08-27 00:00:00", tz = "UTC")
end   <- as.POSIXct("2023-08-26 00:00:00", tz = "UTC")

master_aq_join <- aq_sensor_df_ll |>
  select(SensorID = sensor_id,
         type,
         category) |>
  left_join(master_aq_raw, by = "SensorID") |>
  filter(category != 'Other') |>
  group_by(SensorID) |>
  filter(max(DateTime) >= '2023-08-26 00:00:00' & 
         min(DateTime) <= '2022-09-27 00:00:00') |>
  filter(DateTime >= start, DateTime <= end) |>
  filter(mean(is.na(PM25)) | mean(is.na(NO2)) <= 0.10) |>
  filter(!(SensorID == "SCC_GH6" & !is.na(NO2) &
           hour(DateTime) >= 1 & hour(DateTime) < 5))

rm(master_aq_raw)

#summarise percentage of NA for each pollutants and sensor
na_summary <- master_aq_join |>
  ungroup() |>
  group_by(SensorID) |>
  select(SensorID, PM25, NO2) |>
  summarise(
    PM25_NA = mean(is.na(PM25)*100),
    NO2_NA = mean(is.na(NO2)*100),
    COMP_NO2 = 100 - NO2_NA,
    COMP_PM25 = 100 - PM25_NA
  )


# read in and QA traffic flow data -----------------------------------------------
traffic_sensor_names <- list.files("Data/Traffic", pattern = "\\.csv$", full.names = TRUE)
names(traffic_sensor_names) <- basename(traffic_sensor_names)

#load in all traffic csvs into a single dataframe
traffic_raw <- map_dfr(traffic_sensor_names, read_csv, .id = "file_id") |>
  select(-file_id) |>
  mutate(flows.sensor = str_sub(flows.sensor, start = 5)) |>
  rename(sensorID = flows.sensor) |>
  mutate(sensorID = str_replace_all(sensorID, "[^[:alnum:]]", ""))

#join to sensor metadata and filter for only sensors with at least 6 months of data and no more than 10% missing values

master_tf_join <- tf_sensor_df_ll |>
  select(sensorID, category) |>
  inner_join(traffic_raw, by = 'sensorID') |>
  filter(category != 'Other') |>
  rename(DateTime = 'flows.universalTime',
         flow = 'flows.flow',
         QC = 'flows.flow:QC') |>
  mutate(flow = if_else(QC == 8, NA, flow)) |>
  group_by(sensorID) |>
  filter(max(DateTime) >= end, min(DateTime) <= start) |>
  filter(DateTime >= start, DateTime <= end) |>
  filter(mean(is.na(flow)) <= 0.10) 

  
#remove unneeded datasets
rm(traffic_raw)

##check for NAs at the road level (min two reading per hour to not be NA)
tf_NA_check <- master_tf_join |>
  inner_join(sensor_to_road_RDD, by = "sensorID") |>
  ungroup() |>
  group_by(ref) |>
  complete(DateTime = seq(min(DateTime), max(DateTime), by = "30 min")) |>
  summarise(
    Is_NA = mean(is.na(flow)*100),
    complete = 100 - Is_NA
  )

#aggregate the traffic data to be per hour by grouping the timestamps into hourly bins
#so it can be weather-normalised
#then join to the sensor metadata, only including sensors on primary/secondary roads and inside the CAZ

master_tf_join_hourly <- master_tf_join |>
  mutate(hr = floor_date(DateTime, "hour")) |> 
  group_by(sensorID, hr, category) |>
  summarise(cars_per_hour = sum(5*flow, na.rm = TRUE), .groups = "drop") |>
  arrange(sensorID, hr) |>
  inner_join(sensor_to_road_RDD, by = "sensorID") |>
  select(-c('match_method', 'match_dist_m', 'category.y', 'road_id')) |>
  rename(date = hr)

rm(master_tf_join)


  
# recreate maps with only included sensors --------------------------------
aq_sensor_df_ll_RDD <- aq_sensor_df_ll |>
  filter(sensor_id %in% master_aq_join$SensorID) |>
  mutate(category = case_when(
    sensor_id == 'SCC_GH4' ~ 'Inside CAZ',
    sensor_id == 'AMF_2450229' ~ 'Inside CAZ',
    TRUE ~ category
  ))

tf_sensor_df_ll_RDD <- tf_sensor_df_ll |>
  filter(sensorID %in% master_tf_join_hourly$sensorID)

basemap +
  geom_polygon(
    data    = df_ring,
    aes(x    = lon, y = lat, group = group),
    fill  = alpha("lightblue", 0.6),
    color   = "black",
    size    = 0.3,
  ) +
  # the CAZ polygon
  geom_polygon(
    data    = df_caz,
    aes(x    = lon, y = lat, group = group),
    fill  = alpha("pink", 0.5),
    color   = "red",
    size    = 0.3
  ) +
  geom_point(
    data = aq_sensor_df_ll_RDD,
    aes(x   = lon, y = lat, color = category),
    size = 3, alpha = 0.8
  ) +
  geom_label_repel(
    data    = aq_sensor_df_ll_RDD,
    aes(x = lon, y = lat, label = sensor_id),
    point.padding = 0,
    label.padding = 0.15,
    size    = 2,             
    family = 'Roboto Condensed',
    min.segment.length = 0,
    colour  = "black"
  ) +
  scale_color_manual(labels = c('CAZ', 'Spillover'),
                     values = c("Inside CAZ" = "darkred", 
                                "CAZ Adjacent" = "darkblue")) +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_blank()) 

basemap +
  geom_polygon(
    data    = df_ring,
    aes(x    = lon, y = lat, group = group),
    color   = "black",
    fill = alpha("lightblue", 0.1),
    size    = 0.3,
  ) +
  geom_polygon(
    data    = df_caz,
    aes(x    = lon, y = lat, group = group),
    color   = "red",
    fill   = alpha("pink", 0.2),
    size    = 0.3
  ) +
  geom_point(
    data = tf_sensor_df_ll_RDD,
    aes(x   = lon, y = lat, color = category),
    size = 2, alpha = 0.8
  ) +
  scale_color_manual(labels = c('CAZ', 'Spillover'),
                     values = c("Inside CAZ" = "darkred", 
                                "CAZ Adjacent" = "darkblue")) +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_blank())


#plot only roads with sensors included in study and at least 4 sensors per ref
roads_map_RDD <- roads_map |>
  filter(ref %in% master_tf_join_hourly$ref) 
  
basemap +
  geom_sf(data = roads_map_RDD,
          aes(color = high_cat),
          inherit.aes = FALSE, linewidth = 0.9) +
  coord_sf(xlim = c(sheffield[1], sheffield[3]),
           ylim = c(sheffield[2], sheffield[4]),
           expand = FALSE, datum = NA) +
  geom_point(
     data = tf_sensor_df_ll_RDD,
     aes(x  = lon, y = lat),
     size = 2, alpha = 0.8
   ) +
  geom_polygon(
    data    = df_caz,
    aes(x    = lon, y = lat, group = group),
    fill  = alpha("pink", 0.5),
    color   = "red",
    linewidth    = 1
  ) +
  scale_color_manual(labels = c('Primary roads', 'Secondary roads'), 
                     values = c("primary" = "darkblue", 
                                "secondary" = "limegreen")) +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 8, family = 'Roboto condensed'))


