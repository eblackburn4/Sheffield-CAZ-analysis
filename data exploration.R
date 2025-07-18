
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

# read in AQ data from each sensor family --------------------------------------------------------

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
  select(SensorID, DateTime, NO2, PM25)

# eWatch data
EW_aq_raw <- map_dfr(EW_sensor_names, read_csv, .id = "file_id") |>
  select(data.sensor,
         data.time, 
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
         amon.time, 
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

#QC overview
AM_aq_raw |>
  filter(NO2_QC != 0) |>
  group_by(NO2_QC) |>
  count()


#various plots

SCC_aq_raw |>
  filter(SensorID == 'Pond Hill') |>
  filter(DateTime > '2023-02-10 01:00:00' & DateTime < '2023-03-27 00:00:00') |>
  ggplot(aes(x = DateTime, y = NO2)) +
  geom_line() +
  labs(title = "SCC NO2 Readings", x = "DateTime", y = "NO2 Value") +
  scale_x_datetime(
    date_labels = "%d/%m/%Y",
    date_breaks = "5 days"
  ) +
  theme_minimal()




#combine into one master df
master_aq_raw <- bind_rows(
  SCC_aq_raw |> mutate(family = "SCC"),
  EW_aq_raw   |> mutate(family = "eWatch"),
  defra_aq_raw |> mutate(family = "Defra"),
  AM_aq_raw   |> mutate(family = "AQMesh")
)

# read in traffic flow data -----------------------------------------------
traffic_sensor_names <- list.files("Data/Traffic", pattern = "\\.csv$", full.names = TRUE)
names(traffic_sensor_names) <- basename(traffic_sensor_names)

#load in all traffic csvs into a single dataframe
traffic_raw <- map_dfr(traffic_sensor_names, read_csv, .id = "file_id") |>
  select(-file_id)

traffic_raw |>
  filter(`flows.flow:QC` != 0) |>
  group_by(`flows.flow:QC`) |>
  count()


# SCC data uptime plots ---------------------------------------------------

plot_sensor_uptime <- function(df, 
                               pollutant,
                               downtime_hours = 24,
                               date_breaks   = "6 months",
                               date_labels   = "%d/%m/%Y") {
  
  # 1. Pivot to long, blank-out non-positive values
  status <- df |>
    pivot_longer(
      cols       = -DateTime,
      names_to   = c("location","pollutant", NA, NA),
      names_sep  = ", ",
      values_to  = "value"
    ) |>
    mutate(
      value = ifelse(value <= 0, NA_real_, value)
    ) |>
    group_by(location, DateTime) |>
    summarise(
      active = any(!is.na(value)),
      .groups = "drop"
    )
  
  # 2. Compute runs of up/down
  status_runs <- status |>
    arrange(location, DateTime) |>
    group_by(location) |>
    mutate(run = cumsum(active != lag(active, default = first(active)))) |>
    group_by(location, run, active) |>
    summarise(
      start     = first(DateTime),
      end       = last(DateTime),
      length_hr = as.numeric(difftime(end, start, units = "hours")) + 1,
      .groups   = "drop"
    ) |>
    mutate(
      status = case_when(
        active == TRUE                             ~ "up",
        active == FALSE & length_hr >= downtime_hours ~ "down",
        TRUE                                       ~ "up"
      )
    )
  
  # 3. Gantt-style plot
  ggplot(status_runs, aes(y = location)) +
    geom_linerange(aes(xmin = start, xmax = end, color = status),
                   size = 6, alpha = 0.8) +
    scale_color_manual(values = c(up = "forestgreen", down = "firebrick")) +
    scale_x_datetime(date_labels = date_labels, date_breaks = date_breaks) +
    labs(
      x     = "Time",
      y     = "Sensor location",
      color = "Status",
      title = paste0(pollutant, " — uptimes/downtimes")
    ) +
    theme_minimal()
}

SCC_aq_long <- SCC_aq_raw |>
  select(-file_id) |>            
  group_by(DateTime) |>         
  summarise(
    across(
      .cols   = everything(),
      .fns    = ~ first(na.omit(.x))
    ),
    .groups = "drop"
  )

#split into one table per pollutant
pollutants <- c("NO", "NO2", "PM25", "PM10")

make_tbl <- function(df, pollutant){
  df |>
    select(DateTime, matches(pollutant)) |>
    rename_with(~ str_remove(.x, "\\s*,.*"), -DateTime)
}

SCC_NO   <- make_tbl(SCC_aq_long, "NO, 000\\[M\\]")
SCC_NO2  <- make_tbl(SCC_aq_long, "NO2, DIF\\[M\\]")
SCC_PM25 <- make_tbl(SCC_aq_long, "PM25, 000\\[M\\]")
SCC_PM10 <- make_tbl(SCC_aq_long, "PM10, 000\\[M\\]")


#plot for all SCC sensors
plot_sensor_uptime(SCC_NO, 'NO')
plot_sensor_uptime(SCC_NO2, 'NO2')
plot_sensor_uptime(SCC_PM25, 'PM25')
plot_sensor_uptime(SCC_PM10, 'PM10')


# AQmesh uptime plots  ------------------------------------------------------------
AM_sensor_names <- list.files("Data/AMfixed", pattern = "\\.csv$", full.names = TRUE)
names(AM_sensor_names) <- basename(AM_sensor_names)

#load in all AQM csvs into a single dataframe
AM_aq_raw <- map_dfr(AM_sensor_names, read_csv, .id = "file_id")

AM_aq_raw |>
  filter(`amon.PM25:QC` != 0) |>
  group_by(`amon.PM25:QC`) |>
  count()

AM_aq_raw |>
  filter(`amon.NO2:QC` != 0) |>
  group_by(`amon.NO2:QC`) |>
  count()


#tag based on pollutant measured
AM_aq_raw <- AM_aq_raw |>
  group_by(amon.sensor) |>
  mutate(
    measures_PM25 = any(!is.na(amon.PM25)),
    measures_NO2  = any(!is.na(amon.NO2)),
    pollutant_type = case_when(
      measures_PM25 & measures_NO2  ~ "both",
      measures_PM25 & !measures_NO2 ~ "PM25 only",
      !measures_PM25 & measures_NO2 ~ "NO2 only",
      TRUE                           ~ "none"
    )
  ) |>
  ungroup()

# build the “status_runs” table (as before) but drop the old status‐color mapping
threshold <- 10
AM_aq_status <- AM_aq_raw |>
  rename(sensor = amon.sensor,
         time   = amon.universalTime) |>
  pivot_longer(
    cols      = c(amon.PM25, amon.NO2),
    names_to  = "pollutant",
    values_to = "value"
  ) |>
  mutate(
    value = ifelse(is.na(value) | value <= 0, NA_real_, value)
  ) |>
  group_by(sensor, time) |>
  summarise(
    active         = any(!is.na(value)),
    pollutant_type = first(pollutant_type),
    .groups        = "drop"
  ) |>
  arrange(sensor, time) |>
  group_by(sensor, pollutant_type) |>
  mutate(run = cumsum(active != lag(active, default = first(active)))) |>
  group_by(sensor, pollutant_type, run, active) |>
  summarise(
    start  = first(time),
    end    = last(time),
    length = n(),
    .groups = "drop"
  )

# 2. single Gantt‐style plot, coloring each sensor’s bars by pollutant_type
ggplot(AM_aq_status, aes(y = sensor)) +
  geom_linerange(
    aes(xmin = start, xmax = end, color = pollutant_type),
    size = 6, alpha = 0.8
  ) +
  geom_vline(
    xintercept = as.POSIXct("2023-02-27 00:00:00", tz = "Europe/London"),
    linetype   = "dotted"
  ) +
  scale_color_manual(
    values = c(
      "PM25 only" = "steelblue",
      "NO2 only"  = "darkorange",
      "both"      = "purple"
    )
  ) +
  scale_x_datetime(
    date_labels = "%d/%m/%Y",
    date_breaks = "6 months"
  ) +
  labs(
    title = "AQmesh Sensor uptime/downtime by pollutant type",
    x     = "Time",
    y     = "Sensor",
    color = "Measures"
  ) +
  theme_minimal()



# Envirowatch uptime plots ------------------------------------------------

ENV_sensor_names <- list.files("Data/eWatch", pattern = "\\.csv$", full.names = TRUE)
names(ENV_sensor_names) <- basename(ENV_sensor_names)

#load in all AQM csvs into a single dataframe
ENV_aq_raw <- map_dfr(ENV_sensor_names, read_csv, .id = "file_id") |>
  select(c(data.sensor, data.universalTime,data.NO2))

ENV_status <- ENV_aq_raw |>
  rename(
    sensor = data.sensor,
    time   = data.universalTime,
    no2    = data.NO2
  ) |>
  arrange(sensor, time) |>
  group_by(sensor) |>
  mutate(
    active = !is.na(no2),
    run    = cumsum(active != lag(active, default = first(active)))
  ) |>
  group_by(sensor, run, active) |>
  summarise(
    start  = first(time),
    end    = last(time),
    length = n(),
    .groups = "drop"
  ) |>
  mutate(
    status = ifelse(!active & length > threshold, "down", "up")
  )

# plot Gantt‐style uptime/downtime chart
ggplot(ENV_status, aes(y = sensor)) +
  geom_linerange(
    aes(xmin = start, xmax = end, color = status),
    size = 6, alpha = 0.8
  ) +
  geom_vline(
    xintercept = as.POSIXct("2023-02-27 00:00:00", tz = "Europe/London"),
    linetype   = "dotted"
  ) +
  scale_x_datetime(
    date_labels = "%d/%m/%Y",
    date_breaks = "3 months"
  ) +
  labs(
    title = "Ewatch Sensor uptime (NO only)",
    x     = "Time",
    y     = "Sensor",
    color = "Status"
  ) +
  theme_minimal() +
  theme(legend.position = 'none') 

