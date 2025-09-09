
## ---------------------------
## Purpose of script: weather normalisation
## 
## Key outputs: 
## - master_aq_ERA5: air quality data with joined ERA5 weather data
## - AQ_norm_list: list of normalised air quality dataframes for each sensor/pollutant pair with rmweather model info
## - traffic_norm_list: list of normalised traffic dataframes for each road with rmweather model info
## 
## Author: Ned Blackburn
## Date Created: 2025-07-12

set.seed(9999) #for reproducibility of random forest model
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
      verbose = TRUE,
      se = FALSE
    )
}

#generate normalised time series for each pollutant/sensor pair

GH4_NO2<- weather_norm("SCC_GH4", "NO2") 
GH3_NO2 <- weather_norm("SCC_GH3", "NO2")
GH3_PM25 <- weather_norm("SCC_GH3", "PM25")
GH6_NO2 <- weather_norm("SCC_GH6", "NO2")
GH6_PM25 <- weather_norm("SCC_GH6", "PM25")
DFR1027_NO2 <- weather_norm("DFR_1027A", "NO2")
DFR1027_PM25 <- weather_norm("DFR_1027A", "PM25")
DFR1063_NO2 <- weather_norm("DFR_1063A", "NO2")
DFR1063_PM25 <- weather_norm("DFR_1063A", "PM25")
AMF245_NO2 <- weather_norm("AMF_2450229", "NO2")
AMF245_PM25 <- weather_norm("AMF_2450229", "PM25")

#put into a list for easier batch processing later on

AQ_norm_list <- list(  GH4_NO2  = GH4_NO2,
                       GH3_NO2  = GH3_NO2,
                       GH3_PM25 = GH3_PM25,
                       GH6_NO2  = GH6_NO2,
                       GH6_PM25 = GH6_PM25,
                       DFR1027_NO2  = DFR1027_NO2,
                       DFR1027_PM25 = DFR1027_PM25,
                       DFR1063_NO2  = DFR1063_NO2,
                       DFR1063_PM25 = DFR1063_PM25,
                       AMF245_PM25 = AMF245_PM25
                       )

#aggregate normalised data to daily means for RDD analysis and set running variable
CAZ_start <- as_date("2023-02-27")

daily_avg <- function(sensor) {
  sensor$normalised <- sensor$normalised |>
    mutate(day = as_date(date)) |>
    group_by(day) |>
    summarise(mean_value = mean(value_predict),
              median_value = median(value_predict), 
              .groups = "drop") |>
    complete(day = seq(min(day), max(day), by = "day")) |>
    mutate(t = as.numeric(day - CAZ_start))
  sensor
} 

#apply daily average to all sensors

AQ_norm_list <- map(AQ_norm_list, daily_avg)
list2env(AQ_norm_list, envir = .GlobalEnv)


# Traffic normalisation ---------------------------------------------------
#apply rmweather to averaged readings within sensor groups

traffic_norm <- function(road){
  master_tf_join_hourly |>
    filter(ref == road) |>
    group_by(date) |>
    summarise(cars_per_hour = mean(cars_per_hour, na.rm = TRUE), 
              .groups = 'drop') |>
    rmw_prepare_data(value = 'cars_per_hour', na.rm = TRUE) |> 
    rmw_do_all(
      variables = c(
        "date_unix", 
        "day_julian", 
        "weekday",
        'hour'),
      n_trees = 300,
      n_samples = 300,
      verbose = TRUE,
      se = FALSE
    )
}

#generate averaged, normalised time series for each road
refs <- master_tf_join_hourly |> distinct(ref) |> pull(ref)

traffic_norm_list <- refs |>
  set_names() |>                
  map(traffic_norm)   

#aggregate data to be daily sum of cars rather than hourly, add running variable and impute
daily_avg_traffic <- function(road) {
  road$normalised_daily <- road$normalised |>
    mutate(day = as_date(date)) |>
    group_by(day) |>
    summarise(cars_per_day = mean(value_predict, na.rm = TRUE),
              .groups = "drop") |>
    complete(day = seq(min(day), max(day), by = "day")) |>
    mutate(t = as.numeric(day - CAZ_start)
    )       
  road
}

#apply daily average to all roads
traffic_norm_list <- map(traffic_norm_list, daily_avg_traffic)

