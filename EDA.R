
## ---------------------------
## Purpose of script: EDA
## Author: Ned Blackburn
## Date Created: 2025-08-24

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)



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

traffic_EDA <- traffic_norm_list |>
  imap_dfr(~{
    df <- .x$normalised_daily |> select(day, value = cars_per_day)
    
    pre  <- df |> filter(day <  CAZ_start)
    post <- df |> filter(day >= CAZ_start)
    
    tibble(
      sensor        = str_remove(.y, "_[^_]+$"),
      pollutant     = coalesce(str_extract(.y, "[^_]+$"), "TRAFFIC"),
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



# pre/post boxplots -------------------------------------------------------


# plot grouped boxplot of pollutant concentrations before and after CAZ start

plot_caz_boxplot <- function(AQ_norm_list, start_date = CAZ_start, pollutant){
  pol  <- toupper(pollutant)
  lst  <- AQ_norm_list[str_detect(toupper(names(AQ_norm_list)), paste0("(^|_)", pol, "$"))]
  
  df <- imap_dfr(lst, ~ tibble(
    sensor = str_remove(.y, "_[^_]+$"),
    day    = .x$normalised$day,
    value  = .x$normalised$mean_value
  )) |>
    mutate(period = if_else(day < start_date, "Pre-CAZ", "Post-CAZ")) |>
    mutate(period = factor(period, levels = c("Pre-CAZ", "Post-CAZ")))
  
  ggplot(df, aes(x = sensor, y = value, fill = period)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(
      x = "Sensor",
      y = 'Mean normalised daily concentration',
      fill = "Period"
    ) +
    theme_ipsum_rc(grid = 'Yy') +
    theme(panel.border = element_rect(color = "grey40",
                                      fill = NA,
                                      size = 0.8))
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
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(
      x     = "Road",
      y     = "Daily mean cars/hour",
      fill  = "Period"
    ) +
    theme_ipsum_rc(grid = 'Yy') +
    theme(panel.border = element_rect(color = "grey40",
                                      fill = NA,
                                      size = 0.8))
}


# plot the grouped boxplot
plot_caz_traffic_boxplot(traffic_norm_list)


#normalised time series
#facet plot of normalised air quality time series from AQ_norm_list
facet_plot_by_pollutant <- function(AQ_norm_list, pollutant){
  
  pol <- toupper(pollutant)
  nm  <- names(AQ_norm_list)
  
  # keep only elements whose name ends with "_<pollutant>"
  keep_idx <- str_detect(toupper(nm), regex(paste0("(^|_)", pol, "$"), ignore_case = TRUE))
  lst <- AQ_norm_list[keep_idx]
  if (length(lst) == 0) stop("No matches for pollutant: ", pollutant)
  
  df <- imap_dfr(lst, ~ tibble(
    sensor = str_remove(.y, "_[^_]+$"),     # drop pollutant suffix
    t      = .x$normalised$t,
    value  = .x$normalised$mean_value
  )) |>
    drop_na(t, value)
  
  ggplot(df, aes(t, value, group = sensor)) +
    geom_line() +
    facet_wrap(~ sensor, scales = "free_y") +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, .05))) +
    labs(x = 'days before/after CAZ',
         y = 'Normalised pollutant concentration (ug/m3)') +
    theme_ipsum_rc(axis_title_size = 11, axis_text_size = 8)
}

facet_plot_by_pollutant(AQ_norm_list, "NO2")
facet_plot_by_pollutant(AQ_norm_list, "PM25")


#facet plot of normalised traffic time series

facet_plot_traffic <- function(traffic_norm_list){
  
  df <- imap_dfr(traffic_norm_list, ~ tibble(
    road = .y,
    t    = .x$normalised_daily$t,
    value = .x$normalised_daily$cars_per_day
  )) 
  
  ggplot(df, aes(t, value, group = road)) +
    geom_line() +
    facet_wrap(~ road, scales = "free_y") +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, .05))) +
    labs(x = 'days before/after CAZ',
         y = 'Normalised daily mean cars/hour)') +
    theme_ipsum_rc(axis_title_size = 11, axis_text_size = 8)
}

facet_plot_traffic(traffic_norm_list)


# # combined weather diagnostics: AQ  --------

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
  labs(
    title = paste('Variable Importance for NO2')
  )

#partial dependencies plot (NO2)
rmw_partial_dependencies(
  model = rmw_aggregated_NO2$model,
  df = rmw_aggregated_NO2$observations,
  variable = NA
) |> 
  rmw_plot_partial_dependencies() 

#same for PM2.5

rmw_aggregated_PM25 <- master_aq_ERA5 |> 
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
  labs(
    title = paste('Variable Importance for PM25')
  )

#partial dependencies plot (PM25)
rmw_partial_dependencies(
  model = rmw_aggregated_PM25$model,
  df = rmw_aggregated_PM25$observations,
  variable = NA
) |> 
  rmw_plot_partial_dependencies() 

# # combined weather diagnostics: traffic  --------
rmw_aggregated_traffic <- master_tf_join_hourly |> 
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
  labs(
    title = paste('Variable Importance for Traffic')
  )

#partial dependencies plot (traffic)
rmw_partial_dependencies(
  model = rmw_aggregated_traffic$model,
  df = rmw_aggregated_traffic$observations,
  variable = NA
) |> 
  rmw_plot_partial_dependencies()

