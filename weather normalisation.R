
## ---------------------------
## Purpose of script: weather normalisation
## Author: Ned Blackburn
## Date Created: 2025-07-12

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes) #nice graphs
library(rmweather) #for meteorological normalisation
library(ranger) #dependency for rmweather
library(weathermetrics) #for calculating humidity from dewpoint

## ---------------------------

# load in ERA5 weather data -----------------------------------------------

#define functions to calc wind speed and dir from u and v components
windDir <-function(u,v){
  (270-atan2(u,v)*180/pi)%%360 
}

windSpd <-function(u,v){
  sqrt(u^2+v^2)
}

#load in data and derive additional vars
ERA5_source <- list.files("Data/ERA5_weather", full.names = TRUE)
ERA5_data <- ERA5_source |>
  map(read_csv) |>
  reduce(full_join, by = c("valid_time",'latitude','longitude')) |>
  select(
    DateTime = valid_time,
    total_rain = tp, 
    u10,
    v10,
    surface_pressure = sp,
    temp_2m = t2m,
    dewpoint_2m = d2m,
    solar_rads = ssrd) |>
  mutate(temp_2m = temp_2m - 273.15, dewpoint_2m = dewpoint_2m - 273.15) |>
  mutate(
    wd_spd = windSpd(u10, v10),
    wd_dir = windDir(u10, v10),
    rel_hum = dewpoint.to.humidity(dewpoint_2m, temp_2m, temperature.metric = 'celsius'))|>
  mutate(rel_hum = ifelse(is.na(rel_hum), 100, rel_hum)) |>
  select(-u10, -v10, -dewpoint_2m) 
  
#Join ERA5 data to AQ sensor data
master_aq_ERA5 <- left_join(master_aq_join, ERA5_data, by = "DateTime") |>
  rename(date = DateTime)

#function for rmweather normalisation pipeline

weather_norm <- function(sensor, pollutant){
  master_aq_ERA5 |> 
    filter(SensorID == sensor) |>
    rmw_prepare_data(value = pollutant, na.rm = TRUE) |> 
    rmw_do_all(
      variables = c(
        "date_unix", 
        "day_julian", 
        "weekday",
        'hour', 
        "total_rain", 
        "surface_pressure", 
        "temp_2m", 
        "solar_rads",
        "wd_spd",
        'wd_dir',
        'rel_hum'
      ),
      n_trees = 300,
      n_samples = 300,
      verbose = TRUE
    )
}


#function for rmweather diagnostic plots

norm_plots <- function(sensor_norm) {
  sensor_norm$model |>
    rmw_model_importance() |> 
    rmw_plot_importance() +
    theme_ipsum_rc()
  
  rmw_predict_the_test_set(
    model = sensor_norm$model,
    df = sensor_norm$observations
  ) |>
    rmw_plot_test_prediction()
  
  rmw_plot_normalised(sensor_norm$normalised) +
    geom_vline(xintercept = as.POSIXct('2023-02-27 00:00:00', tz = "UTC"), linetype = 'dashed', color = 'red')
  
  ggplot(data = sensor_norm$observations, aes(x = date, y = value)) +
    geom_line()
  
  rmw_partial_dependencies(
    model = sensor_norm$model, 
    df = sensor_norm$observations,
    variable = NA
  ) |>
    filter(variable != "date_unix") |>
    rmw_plot_partial_dependencies() 
}

#generate normalised time series and diagnostic plots for each pollutant/sensor pair

SCCGH4_NO2 <- weather_norm("SCC_GH4", "NO2")
norm_plots(SCCGH4_NO2)


