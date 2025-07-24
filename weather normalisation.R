
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


# example from docs -------------------------------------------------------

# Have a look at rmweather's example data, from london
ldn <- data_london 

# Prepare data for modelling
# Only use data with valid wind speeds, no2 will become the dependent variable
data_london_prepared <- ldn %>% 
  filter(variable == "no2",
         !is.na(ws)) %>% 
  rmw_prepare_data(na.rm = TRUE)

# Grow/train a random forest model and then create a meteorological normalised trend 
list_normalised <- rmw_do_all(
  data_london_prepared,
  variables = c(
    "date_unix", "day_julian", "weekday", "air_temp", "rh", "wd", "ws",
    "atmospheric_pressure"
  ),
  n_trees = 300,
  n_samples = 300,
  verbose = TRUE
)

# What units are in the list? 
names(list_normalised)

# Check model object's performance
rmw_model_statistics(list_normalised$model)

# Plot variable importances
list_normalised$model %>% 
  rmw_model_importance() %>% 
  rmw_plot_importance()

# Check if model has suffered from overfitting
rmw_predict_the_test_set(
  model = list_normalised$model,
  df = list_normalised$observations
) %>% 
  rmw_plot_test_prediction()

# How long did the process take? 
list_normalised$elapsed_times

# Plot normalised trend
rmw_plot_normalised(list_normalised$normalised)

# Investigate partial dependencies, if variable is NA, predict all
data_pd <- rmw_partial_dependencies(
  model = list_normalised$model, 
  df = list_normalised$observations,
  variable = NA
)

# Plot partial dependencies
data_pd %>% 
  filter(variable != "date_unix") %>% 
  rmw_plot_partial_dependencies()



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
  




