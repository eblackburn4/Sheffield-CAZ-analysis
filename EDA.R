
## ---------------------------
## Purpose of script: EDA
## Author: Ned Blackburn
## Date Created: 2025-08-24

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)



# Example plots of norm vs non norm timeseries -----------------------------

GH4_NO2$normalised |>
  ggplot(aes(x = date, y = value_predict)) +
  geom_line() +
  theme_ipsum_rc(grid = 'XY',axis_title_size = 10, axis_text_size = 10) +
  labs(x = 'Date', y = 'Hourly normalised pollutant concentration (ug/m3)')

GH4_NO2$observations |>
  ggplot(aes(x = date, y = value)) +
  geom_line() +
  theme_ipsum_rc(grid = 'XY',axis_title_size = 10, axis_text_size = 10) +
  labs(x = 'Date', y = 'Hourly pollutant concentration (ug/m3)')



# Summary statistics ------------------------------------------------------

aq_EDA_stats <- AQ_norm_list |>
  imap_dfr(~{
    df <- .x$normalised |> select(day, value = mean_value)
    
    pre  <- df |> filter(day <  CAZ_start)
    post <- df |> filter(day >= CAZ_start)
    
    tibble(
      sensor        = str_remove(.y, "_[^_]+$"),
      pollutant     = str_extract(.y, "[^_]+$"),
      pre_mean      = mean(pre$value,  na.rm = TRUE),
      pre_sd        = sd(pre$value,    na.rm = TRUE),
      post_mean     = mean(post$value, na.rm = TRUE),
      post_sd       = sd(post$value,   na.rm = TRUE),
      pct_missing   = 100 * mean(is.na(df$value))
    ) |>
      mutate(
        pct_change_mean = if_else(is.finite(pre_mean) & pre_mean != 0,
                                  100 * (post_mean - pre_mean) / pre_mean,
                                  NA)
      )
  }) |>
  arrange(sensor, pollutant)

tf_EDA_stats <- traffic_norm_list |>
  imap_dfr(~{
    df <- .x$normalised_daily |> select(day, value = cars_per_day)
    
    pre  <- df |> filter(day <  CAZ_start)
    post <- df |> filter(day >= CAZ_start)
    
    tibble(
      sensor        = str_remove(.y, "_[^_]+$"),
      road     = coalesce(str_extract(.y, "[^_]+$"), "TRAFFIC"),
      pre_mean      = mean(pre$value,  na.rm = TRUE),
      pre_sd        = sd(pre$value,    na.rm = TRUE),
      post_mean     = mean(post$value, na.rm = TRUE),
      post_sd       = sd(post$value,   na.rm = TRUE),
      pct_missing   = 100 * mean(is.na(df$value))
    ) |>
      mutate(
        pct_change_mean = if_else(is.finite(pre_mean) & pre_mean != 0,
                                  100 * (post_mean - pre_mean) / pre_mean,
                                  NA)
      )
  }) |>
  arrange(sensor, road)



# pre/post boxplots -------------------------------------------------------

#nice labels for the charts
var_labels_NO2 <- c(
  DFR1027 = "DFR_1027A",
  DFR1063 = "DFR_1063A",
  GH3 = "SCC_GH3",
  GH4 = "SCC_GH4",
  GH6 = "SCC_GH6"
)

var_labels_PM25 <- c(
  AMF245 = "AMF_2450229",
  DFR1027 = "DFR_1027A",
  DFR1063 = "DFR_1063A",
  GH3 = "SCC_GH3",
  GH6 = "SCC_GH6"
)

#order sensors to match tables in thesis
sensor_order_NO2  <- c("GH4",'DFR1027','GH6','DFR1063','GH3')
sensor_order_PM25 <- c("AMF245",'DFR1027','GH6','DFR1063','GH3')
sensor_order_TF <- c("CAZ interior road","A61 Ring Road","A61","A6109","A57",
                     'A6135','A625','A621','B6069','B6388','B6547')

# plot grouped boxplot of pollutant concentrations before and after CAZ start

plot_caz_boxplot <- function(AQ_norm_list, start_date = CAZ_start, pollutant){
  pol  <- toupper(pollutant)
  lst  <- AQ_norm_list[str_detect(toupper(names(AQ_norm_list)), paste0("(^|_)", pol, "$"))]
  labels <- if (pol == 'NO2') var_labels_NO2 else var_labels_PM25
  sensor_order <- if (pol == 'NO2') sensor_order_NO2 else sensor_order_PM25
  
  df <- imap_dfr(lst, ~ tibble(
    sensor = str_remove(.y, "_[^_]+$"),
    day    = .x$normalised$day,
    value  = .x$normalised$mean_value
  )) |>
    mutate(period = if_else(day < start_date, "Pre-CAZ", "Post-CAZ")) |>
    mutate(period = factor(period, levels = c("Pre-CAZ", "Post-CAZ")))
  
  ggplot(df, aes(x = sensor, y = value, fill = period)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    labs(
      x = "Sensor",
      y = 'Mean normalised daily concentration (ug/m3)',
      fill = "Period"
    ) +
    theme_ipsum_rc(grid = 'Yy', axis_text_size = 10, axis_title_size = 10) +
    theme(panel.border = element_rect(color = "grey40",
                                      fill = NA,
                                      size = 0.8)) +
    scale_x_discrete(limits = sensor_order, labels = labels,
                     guide = guide_axis(angle = 90))
         
}
# plot the grouped boxplot
plot_caz_boxplot(AQ_norm_list, pollutant = 'NO2')
plot_caz_boxplot(AQ_norm_list, pollutant = 'PM25')


#grouped boxplot of traffic sensor readings before and after CAZ start

plot_caz_traffic_boxplot <- function(traffic_norm_list, start_date = CAZ_start) {
  traffic_norm_list |>
    imap_dfr(function(road, nm) {
      road$normalised_daily |>
        mutate(
          road   = nm,
          period = if_else(t > 0, "Post-CAZ", "Pre-CAZ")
        ) |>
        select(road, period, cars_per_day) |>
        mutate(period = factor(period,
                               levels = c("Pre-CAZ", "Post-CAZ")))
    }) |>
    ggplot(aes(x = road, y = cars_per_day, fill = period)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    labs(
      x     = "Road",
      y     = "Daily mean cars/hour",
      fill  = "Period"
    ) +
    theme_ipsum_rc(grid = 'Yy') +
    theme(panel.border = element_rect(color = "grey40",
                                      fill = NA,
                                      size = 0.8)) +
    scale_x_discrete(limits = sensor_order_TF, labels = var_labels_tf,
                     guide = guide_axis(angle = 90))
}


# plot the grouped boxplot
plot_caz_traffic_boxplot(traffic_norm_list)


#normalised time series


#facet plot of normalised air quality time series from AQ_norm_list
facet_plot_by_pollutant <- function(AQ_norm_list, pollutant, pol_label, order){
  
  pol <- toupper(pollutant)
  nm  <- names(AQ_norm_list)
  sensor_order <- order
  
  
  # keep only elements whose name ends with "_<pollutant>"
  keep_idx <- str_detect(toupper(nm), regex(paste0("(^|_)", pol, "$"), ignore_case = TRUE))
  lst <- AQ_norm_list[keep_idx]
  
  df <- imap_dfr(lst, ~ tibble(
    sensor = str_remove(.y, "_[^_]+$"),    
    t      = .x$normalised$t,
    value  = .x$normalised$mean_value
  )) |>
    mutate(sensor = factor(sensor, levels = sensor_order)) |>
    mutate(Location = if_else(sensor %in% c("GH4", "DFR1027", 'GH6','AMF245'), "CAZ", "Spillover"))
    
  
  ggplot(df, aes(t, value, group = sensor, color = Location)) +
    geom_line() +
    facet_wrap(~sensor, scales = "free_y", ncol = 5,
               labeller = as_labeller(pol_label)) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, .05))) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey30') +
    labs(x = 'Days before/after CAZ introduction',
         y = 'Normalised pollutant concentration (ug/m3)') +
    scale_colour_manual(values = c("CAZ" = "darkred", "Spillover" = "darkblue")) +
    theme_ipsum_rc(axis_title_size = 9, axis_text_size = 8) +
    theme(legend.position = "bottom", legend.title = element_blank())
}

facet_plot_by_pollutant(AQ_norm_list, "NO2", var_labels_NO2, sensor_order_NO2)
facet_plot_by_pollutant(AQ_norm_list, "PM25", var_labels_PM25, sensor_order_PM25)


#facet plot of normalised traffic time series

facet_plot_traffic <- function(traffic_norm_list){
  
  df <- imap_dfr(traffic_norm_list, ~ tibble(
    road = .y,
    t    = .x$normalised_daily$t,
    value = .x$normalised_daily$cars_per_day
  )) |>
    mutate(road = factor(road, levels = sensor_order_TF)) |>
    mutate(Location = if_else(road %in% c("CAZ interior road", 'A61 Ring Road'), 
                              "CAZ", "Spillover"))
  
  ggplot(df, aes(t, value, group = road, color = Location)) +
    geom_line() +
    facet_wrap(~ road, scales = "free_y") +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, .05))) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey30') +
    labs(x = 'Days before/after CAZ introduction',
         y = 'Normalised daily mean cars/hour') +
    scale_colour_manual(values = c("CAZ" = "darkred", "Spillover" = "darkblue")) +
    theme_ipsum_rc(axis_title_size = 9, axis_text_size = 8) +
    theme(legend.position = "bottom", legend.title = element_blank())
}

facet_plot_traffic(traffic_norm_list)


# # combined weather diagnostics: AQ  --------

#define nice labels for graph

var_labels <- c(
  date_unix       = "Trend (Unix time)",
  day_julian      = "Day of year",
  weekday         = "Day of week",
  hour            = "Hour of day",
  total_rain      = "Rainfall (m)",
  surface_pressure= "Surface pressure (hPa)",
  temp_2m         = "Temperature (C)",
  solar_rads      = "Solar radiation (J/m2)",
  wd_spd          = "Wind speed (m/s)",
  wd_dir          = "Wind direction (degrees)",
  rel_hum         = "Relative humidity (%)"
)

var_labels_tf <- c(
  date_unix       = "Trend (Unix time)",
  day_julian      = "Day of year",
  weekday         = "Day of week",
  hour            = "Hour of day")


rmw_aggregated_NO2 <- master_aq_ERA5 |>
  rmw_prepare_data(value = 'NO2', na.rm = TRUE) |> 
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

#variable importance plot (NO2)
rmw_aggregated_NO2$model |>
  rmw_model_importance() |> 
  rmw_plot_importance() +
  theme_ipsum_rc(axis_title_size = 12) +
  scale_y_discrete(labels = var_labels) 

#partial dependencies plot (NO2)
NO2_partial <- rmw_partial_dependencies(
  model = rmw_aggregated_NO2$model,
  df = rmw_aggregated_NO2$observations,
  variable = NA
) |> 
  rmw_plot_partial_dependencies() 

NO2_partial$facet$params$labeller <- ggplot2::as_labeller(var_labels)

NO2_partial + 
  labs(x = 'Variable', y = 'Partial dependency') + 
  theme_ipsum_rc(axis_title_size = 12, axis_text_size = 8)
  


#same for PM2.5

rmw_aggregated_PM25 <- master_aq_ERA5 |> 
  filter(SensorID != 'AMF_2450229') |>
  rmw_prepare_data(value = 'PM25', na.rm = TRUE) |> 
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

#variable importance plot (PM25)
rmw_aggregated_PM25$model |>
  rmw_model_importance() |> 
  rmw_plot_importance() +
  theme_ipsum_rc(axis_title_size = 12) +
  scale_y_discrete(labels = var_labels) 
  

#partial dependencies plot (PM25)
PM25_partial <- rmw_partial_dependencies(
  model = rmw_aggregated_PM25$model,
  df = rmw_aggregated_PM25$observations,
  variable = NA
) |> 
  rmw_plot_partial_dependencies() 

PM25_partial$facet$params$labeller <- ggplot2::as_labeller(var_labels)

PM25_partial + 
  labs(x = 'Variable', y = 'Partial dependency') + 
  theme_ipsum_rc(axis_title_size = 12, axis_text_size = 8)




#facet plot of predicted vs actuals
pred_df_aq <- imap_dfr(AQ_norm_list, function(sensor_norm, nm) {
  out <- rmw_predict_the_test_set(
    model = sensor_norm$model,
    df    = sensor_norm$observations
  )
  out$sensor <- nm
  out
}) |>
  mutate(
    pollutant = case_when( 
      grepl("PM25", sensor, ignore.case = TRUE) ~ "PM2.5",
      grepl("NO2",  sensor, ignore.case = TRUE) ~ "NO2",
      TRUE ~ "Other"
    ))
  
#NO2 plot
pred_df_aq |>
  filter(pollutant == 'NO2') |>
  ggplot(aes(x = value, y = value_predict)) +
  geom_hex(fill = 'darkblue') +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~ SensorID, scales = 'free', ncol = 5) +
  labs(
    x = "Actuals",
    y = "Predicted") +
  theme_ipsum_rc(axis_title_size = 10, axis_text_size = 8) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10))

pred_df_aq |>
  filter(pollutant == 'PM2.5') |>
  filter(value < 350) |>
  ggplot(aes(x = value, y = value_predict)) +
  geom_hex(fill = 'darkred') +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~ SensorID, scales = 'free', ncol = 5) +
  labs(
    x = "Actuals",
    y = "Predicted") +
  theme_ipsum_rc(axis_title_size = 10, axis_text_size = 8) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10))

# # combined weather diagnostics: traffic  --------
#use a sample to avoid crashing the model

rmw_aggregated_traffic <- master_tf_join_hourly |> 
  group_by(ref) |>
  slice_sample(prop = 0.2) |>
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

#variable importance plot (traffic)
rmw_aggregated_traffic$model |>
  rmw_model_importance() |> 
  rmw_plot_importance() +
  theme_ipsum_rc(axis_title_size = 12) +
  scale_y_discrete(labels = var_labels) 
  

#partial dependencies plot (traffic) 
tf_partial <- rmw_partial_dependencies(
  model = rmw_aggregated_traffic$model,
  df = rmw_aggregated_traffic$observations,
  variable = NA
) |> 
  rmw_plot_partial_dependencies()

tf_partial$facet$params$labeller <- ggplot2::as_labeller(var_labels_tf)

tf_partial + 
  labs(x = 'Variable', y = 'Partial dependency') + 
  theme_ipsum_rc(axis_title_size = 12, axis_text_size = 8)


#predicted vs actuals: traffic 

pred_df_tf <- imap_dfr(traffic_norm_list, function(sensor_norm, nm) {
  out <- rmw_predict_the_test_set(
    model = sensor_norm$model,
    df    = sensor_norm$observations
  )
  out$sensor <- nm
  out
}) 

pred_df_tf |>
  ggplot(aes(x = value, y = value_predict)) +
  geom_hex(fill = 'darkgreen') +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~sensor, scales = 'free', ncol = 5) +
  labs(
    x = "Actuals",
    y = "Predicted") +
  theme_ipsum_rc(axis_title_size = 10, axis_text_size = 8) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10))



  

