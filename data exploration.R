
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

# read in SCC AQ data --------------------------------------------------------

SCC_sensor_names <- list.files("Data/SCC", pattern = "\\.csv$", full.names = TRUE)
names(SCC_sensor_names) <- basename(SCC_sensor_names)

#load in all SCC csvs into a single dataframe
SCC_aq_raw <- map_dfr(SCC_sensor_names, read_csv, .id = "file_id")


#reshape so one row = one timestamp
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

#reshape for plotting

make_tbl_long <- function(df){
  df |>
    pivot_longer(
      cols      = -DateTime,
      names_to  = c("location", "pollutant", NA, NA),
      names_sep = ", ",
      values_to = "value"
    ) |>
    select(DateTime, location, value)
}

SCC_NO_long   <- make_tbl_long(SCC_NO)
SCC_NO2_long  <- make_tbl_long(SCC_NO2)
SCC_PM25_long <- make_tbl_long(SCC_PM25)
SCC_PM10_long <- make_tbl_long(SCC_PM10)

#plot them

ggplot(SCC_NO_long, aes(x = DateTime, y = value, color = location)) +
  geom_line() +
  facet_wrap(~location, scales = "free_y") +
  labs(
    title = "NO concentrations over time",
    x = "Date",
    y = "NO Concentration (ppb)"
  ) +
  scale_color_viridis_d() +
  theme_tufte()

ggplot(SCC_NO2_long, aes(x = DateTime, y = value, color = location)) +
  geom_line() +
  facet_wrap(~location, scales = "free_y") +
  labs(
    title = "NO2 concentrations over time",
    x = "Date",
    y = "NO Concentration (ppb)"
  ) +
  scale_color_viridis_d() +
  theme_tufte()

ggplot(SCC_PM25_long, aes(x = DateTime, y = value, color = location)) +
  geom_line() +
  facet_wrap(~location, scales = "free_y") +
  labs(
    title = "PM 2.5 concentrations over time",
    x = "Date",
    y = "NO Concentration (ppb)"
  ) +
  scale_color_viridis_d() +
  theme_tufte()

ggplot(SCC_PM10_long, aes(x = DateTime, y = value, color = location)) +
  geom_line() +
  facet_wrap(~location, scales = "free_y") +
  labs(
    title = "PM10 concentrations over time",
    x = "Date",
    y = "NO Concentration (ppb)"
  ) +
  scale_color_viridis_d() +
  theme_tufte()


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

#plot for all SCC sensors
plot_sensor_uptime(SCC_NO, 'NO')
plot_sensor_uptime(SCC_NO2, 'NO2')
plot_sensor_uptime(SCC_PM25, 'PM25')
plot_sensor_uptime(SCC_PM10, 'PM10')


# AQmesh uptime plots  ------------------------------------------------------------
AM_sensor_names <- list.files("Data/AMfixed", pattern = "\\.csv$", full.names = TRUE)
names(AM_sensor_names) <- basename(AM_sensor_names)

#load in all AQM csvs into a single dataframe
AM_aq_raw <- map_dfr(AM_sensor_names, read_csv, .id = "file_id") |>
  select(c(amon.sensor, amon.universalTime,amon.PM25,amon.NO2))

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

# 1. build the “status_runs” table (as before) but drop the old status‐color mapping
threshold <- 10
AM_aq_status <- AM_aq_raw %>%
  rename(sensor = amon.sensor,
         time   = amon.universalTime) %>%
  pivot_longer(
    cols      = c(amon.PM25, amon.NO2),
    names_to  = "pollutant",
    values_to = "value"
  ) %>%
  mutate(
    value = ifelse(is.na(value) | value <= 0, NA_real_, value)
  ) %>%
  group_by(sensor, time) %>%
  summarise(
    active         = any(!is.na(value)),
    pollutant_type = first(pollutant_type),
    .groups        = "drop"
  ) %>%
  arrange(sensor, time) %>%
  group_by(sensor, pollutant_type) %>%
  mutate(run = cumsum(active != lag(active, default = first(active)))) %>%
  group_by(sensor, pollutant_type, run, active) %>%
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

ENV_status <- ENV_aq_raw %>%
  rename(
    sensor = data.sensor,
    time   = data.universalTime,
    no2    = data.NO2
  ) %>%
  arrange(sensor, time) %>%
  group_by(sensor) %>%
  mutate(
    active = !is.na(no2),
    run    = cumsum(active != lag(active, default = first(active)))
  ) %>%
  group_by(sensor, run, active) %>%
  summarise(
    start  = first(time),
    end    = last(time),
    length = n(),
    .groups = "drop"
  ) %>%
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




