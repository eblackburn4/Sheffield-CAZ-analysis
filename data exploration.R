
## ---------------------------
## Purpose of script: data exploration and preparation code for dissertation
## Author: Ned Blackburn
## Date Created: 2025-06-27

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)
library(rmweather)
library(ggthemes)

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
  filter(mean(is.na(PM25)) | mean(is.na(NO2)) <= 0.10) 

#rm(master_aq_raw)

#summarise percentage of NA for each pollutants and sensor
master_aq_join |>
  ungroup() |>
  group_by(SensorID) |>
  select(SensorID, PM25, NO2) |>
  summarise(
    PM25_NA = mean(is.na(PM25)),
    NO2_NA = mean(is.na(NO2))
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

#check presence of NAs
master_tf_join |>
  ungroup() |>
  distinct(sensorID) |>
  count()

master_aq_join |>
  ungroup() |>
  distinct(SensorID) |>
  count()

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
  filter(sensor_id %in% master_aq_join$SensorID)

tf_sensor_df_ll_RDD <- tf_sensor_df_ll |>
  filter(sensorID %in% master_tf_join_hourly$sensorID)

basemap +
  geom_polygon(
    data    = df_ring,
    aes(x    = lon, y = lat, group = group),
    fill  = alpha("grey90", 0.6),
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
    aes(x   = lon, y = lat, shape = category),
    size = 3, alpha = 0.8
  ) +
  geom_text_repel(
    data    = aq_sensor_df_ll_RDD,
    aes(x = lon, y = lat, label = sensor_id),            # e.g. just above
    size    = 3,             # font size
    colour  = "black"
  ) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) 

basemap +
  # the 500 m exterior ring
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
    size    = 0.3
  ) +
  geom_point(
    data = tf_sensor_df_ll_RDD,
    aes(x   = lon, y = lat, color = category),
    size = 2, alpha = 0.8
  ) +
  scale_color_manual(values = c("Inside CAZ" = "#9c335a", 
                                "CAZ Adjacent" = "darkblue", 
                                "Other" = "grey50")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank())



# 
# #QC overview for both pollutants (excludes SCC)
# master_aq_join |>
#   filter(DateTime > '2022-08-27 00:00:00' & DateTime < '2023-08-27 00:00:00') |>
#   filter(PM25_QC != 0) |>
#   group_by(PM25_QC, SensorID) |>
#   count()
# 
# master_aq_join |>
#   filter(DateTime > '2022-08-27 00:00:00' & DateTime < '2023-08-27 00:00:00') |>
#   filter(NO2_QC != 0) |>
#   group_by(NO2_QC, SensorID) |>
#   count()
# 
# #check for NA values not included in the QC codes
# master_aq_join |>
#   filter(DateTime > '2022-08-27 00:00:00' & DateTime < '2023-08-27 00:00:00') |>
#   filter(type != 'NO2') |>
#   filter(is.na(PM25)) |>
#   group_by(SensorID) |>
#   count()
# 
# #plot template to investigate QC issues
# 
# master_aq_join |>
#   ggplot(aes(x = DateTime, y = PM25)) +
#   geom_line(color = 'darkgreen') +
#   geom_vline(
#     xintercept = as.POSIXct('2023-02-27 00:00:00', tz = "UTC"),
#     linetype   = "dashed",
#     color      = "black"
#   ) +
#   facet_wrap(~SensorID) +
#   labs(title = "PM2.5 Readings", x = "DateTime", y = "NO2 Value") +
#   theme_minimal()
# 
# master_aq_join |>
#   filter(SensorID == 'SCC_GH6') |>
#   mutate(value = ifelse(is.na(NO2), 1, 0)) |>
#   ggplot(aes(x = DateTime, y = value)) +
#   geom_line(color = 'black') +
#   geom_vline(
#     xintercept = as.POSIXct('2023-02-27 00:00:00', tz = "UTC"),
#     linetype   = "dashed",
#     color      = "black"
#   ) +
#   facet_wrap(~SensorID) +
#   labs(title = "NA periods - GH6_NO2", x = "DateTime", y = "NO2 Value") +
#   theme_minimal()
# 


# # # SCC data uptime plots ---------------------------------------------------
# # 
# # plot_sensor_uptime <- function(df,
# #                                pollutant,
# #                                downtime_hours = 24,
# #                                date_breaks   = "6 months",
# #                                date_labels   = "%d/%m/%Y") {
# # 
# #   # 1. Pivot to long, blank-out non-positive values
# #   status <- df |>
# #     pivot_longer(
# #       cols       = -DateTime,
# #       names_to   = c("location","pollutant", NA, NA),
# #       names_sep  = ", ",
# #       values_to  = "value"
# #     ) |>
# #     mutate(
# #       value = ifelse(value <= 0, NA_real_, value)
# #     ) |>
# #     group_by(location, DateTime) |>
# #     summarise(
# #       active = any(!is.na(value)),
# #       .groups = "drop"
# #     )
# # 
# #   # 2. Compute runs of up/down
# #   status_runs <- status |>
# #     arrange(location, DateTime) |>
# #     group_by(location) |>
# #     mutate(run = cumsum(active != lag(active, default = first(active)))) |>
# #     group_by(location, run, active) |>
# #     summarise(
# #       start     = first(DateTime),
# #       end       = last(DateTime),
# #       length_hr = as.numeric(difftime(end, start, units = "hours")) + 1,
# #       .groups   = "drop"
# #     ) |>
# #     mutate(
# #       status = case_when(
# #         active == TRUE                             ~ "up",
# #         active == FALSE & length_hr >= downtime_hours ~ "down",
# #         TRUE                                       ~ "up"
# #       )
# #     )
# # 
# #   # 3. Gantt-style plot
# #   ggplot(status_runs, aes(y = location)) +
# #     geom_linerange(aes(xmin = start, xmax = end, color = status),
# #                    size = 6, alpha = 0.8) +
# #     scale_color_manual(values = c(up = "forestgreen", down = "firebrick")) +
# #     scale_x_datetime(date_labels = date_labels, date_breaks = date_breaks) +
# #     labs(
# #       x     = "Time",
# #       y     = "Sensor location",
# #       color = "Status",
# #       title = paste0(pollutant, " — uptimes/downtimes")
# #     ) +
# #     theme_minimal()
# # }
# # 
# # 
# # #plot for all SCC sensors
# # plot_sensor_uptime(S, 'NO')
# # plot_sensor_uptime(SCC_NO2, 'NO2')
# # plot_sensor_uptime(SCC_PM25, 'PM25')
# # plot_sensor_uptime(SCC_PM10, 'PM10')
# # 
# # 
# # # # AQmesh uptime plots  ------------------------------------------------------------
# # # AM_sensor_names <- list.files("Data/AMfixed", pattern = "\\.csv$", full.names = TRUE)
# # # names(AM_sensor_names) <- basename(AM_sensor_names)
# # #
# # # #load in all AQM csvs into a single dataframe
# # # AM_aq_raw <- map_dfr(AM_sensor_names, read_csv, .id = "file_id")
# # #
# # # AM_aq_raw |>
# # #   filter(`amon.PM25:QC` != 0) |>
# # #   group_by(`amon.PM25:QC`) |>
# # #   count()
# # #
# # # AM_aq_raw |>
# # #   filter(`amon.NO2:QC` != 0) |>
# # #   group_by(`amon.NO2:QC`) |>
# # #   count()
# # #
# # #
# # # #tag based on pollutant measured
# # # AM_aq_raw <- AM_aq_raw |>
# # #   group_by(amon.sensor) |>
# # #   mutate(
# # #     measures_PM25 = any(!is.na(amon.PM25)),
# # #     measures_NO2  = any(!is.na(amon.NO2)),
# # #     pollutant_type = case_when(
# # #       measures_PM25 & measures_NO2  ~ "both",
# # #       measures_PM25 & !measures_NO2 ~ "PM25 only",
# # #       !measures_PM25 & measures_NO2 ~ "NO2 only",
# # #       TRUE                           ~ "none"
# # #     )
# # #   ) |>
# # #   ungroup()
# # #
# # # # build the “status_runs” table (as before) but drop the old status‐color mapping
# # # threshold <- 10
# # # AM_aq_status <- AM_aq_raw |>
# # #   rename(sensor = amon.sensor,
# # #          time   = amon.universalTime) |>
# # #   pivot_longer(
# # #     cols      = c(amon.PM25, amon.NO2),
# # #     names_to  = "pollutant",
# # #     values_to = "value"
# # #   ) |>
# # #   mutate(
# # #     value = ifelse(is.na(value) | value <= 0, NA_real_, value)
# # #   ) |>
# # #   group_by(sensor, time) |>
# # #   summarise(
# # #     active         = any(!is.na(value)),
# # #     pollutant_type = first(pollutant_type),
# # #     .groups        = "drop"
# # #   ) |>
# # #   arrange(sensor, time) |>
# # #   group_by(sensor, pollutant_type) |>
# # #   mutate(run = cumsum(active != lag(active, default = first(active)))) |>
# # #   group_by(sensor, pollutant_type, run, active) |>
# # #   summarise(
# # #     start  = first(time),
# # #     end    = last(time),
# # #     length = n(),
# # #     .groups = "drop"
# # #   )
# # #
# # # # 2. single Gantt‐style plot, coloring each sensor’s bars by pollutant_type
# # # ggplot(AM_aq_status, aes(y = sensor)) +
# # #   geom_linerange(
# # #     aes(xmin = start, xmax = end, color = pollutant_type),
# # #     size = 6, alpha = 0.8
# # #   ) +
# # #   geom_vline(
# # #     xintercept = as.POSIXct("2023-02-27 00:00:00", tz = "Europe/London"),
# # #     linetype   = "dotted"
# # #   ) +
# # #   scale_color_manual(
# # #     values = c(
# # #       "PM25 only" = "steelblue",
# # #       "NO2 only"  = "darkorange",
# # #       "both"      = "purple"
# # #     )
# # #   ) +
# # #   scale_x_datetime(
# # #     date_labels = "%d/%m/%Y",
# # #     date_breaks = "6 months"
# # #   ) +
# # #   labs(
# # #     title = "AQmesh Sensor uptime/downtime by pollutant type",
# # #     x     = "Time",
# # #     y     = "Sensor",
# # #     color = "Measures"
# # #   ) +
# # #   theme_minimal()
# # #
# # #
# # #
# # # # Envirowatch uptime plots ------------------------------------------------
# # #
# # # ENV_sensor_names <- list.files("Data/eWatch", pattern = "\\.csv$", full.names = TRUE)
# # # names(ENV_sensor_names) <- basename(ENV_sensor_names)
# # #
# # # #load in all AQM csvs into a single dataframe
# # # ENV_aq_raw <- map_dfr(ENV_sensor_names, read_csv, .id = "file_id") |>
# # #   select(c(data.sensor, data.universalTime,data.NO2))
# # #
# # # ENV_status <- ENV_aq_raw |>
# # #   rename(
# # #     sensor = data.sensor,
# # #     time   = data.universalTime,
# # #     no2    = data.NO2
# # #   ) |>
# # #   arrange(sensor, time) |>
# # #   group_by(sensor) |>
# # #   mutate(
# # #     active = !is.na(no2),
# # #     run    = cumsum(active != lag(active, default = first(active)))
# # #   ) |>
# # #   group_by(sensor, run, active) |>
# # #   summarise(
# # #     start  = first(time),
# # #     end    = last(time),
# # #     length = n(),
# # #     .groups = "drop"
# # #   ) |>
# # #   mutate(
# # #     status = ifelse(!active & length > threshold, "down", "up")
# # #   )
# # #
# # # # plot Gantt‐style uptime/downtime chart
# # # ggplot(ENV_status, aes(y = sensor)) +
# # #   geom_linerange(
# # #     aes(xmin = start, xmax = end, color = status),
# # #     size = 6, alpha = 0.8
# # #   ) +
# # #   geom_vline(
# # #     xintercept = as.POSIXct("2023-02-27 00:00:00", tz = "Europe/London"),
# # #     linetype   = "dotted"
# # #   ) +
# # #   scale_x_datetime(
# # #     date_labels = "%d/%m/%Y",
# # #     date_breaks = "3 months"
# # #   ) +
# # #   labs(
# # #     title = "Ewatch Sensor uptime (NO only)",
# # #     x     = "Time",
# # #     y     = "Sensor",
# # #     color = "Status"
# # #   ) +
# # #   theme_minimal() +
# # #   theme(legend.position = 'none')
# # #
# # #
# # #
# # #
# # #
# # #
