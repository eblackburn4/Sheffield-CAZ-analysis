# data imputation
# mutate(
#   mean_value   = na.approx(mean_value,   x = day, na.rm = FALSE),
#   median_value = na.approx(median_value, x = day, na.rm = FALSE),
#   t = as.numeric(day - CAZ_start))



## ---------------------------
## Purpose of script: RDD analysis of AQ and traffic data
## Author: Ned Blackburn
## Date Created: 2025-07-26

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)
library(rdrobust)
library(ggpattern)
library(rdpower)
library(zoo)
library(scico)


# RDD estimators for AQ data ----------------------------------------------
#aggregate data to be daily avg rather than hourly, add running variable and impute 
# missing days using linear interpolation

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


#calculate optimal bandwidth, and remove any sensors that have more than 10% of NAs within bandwidth
screen_by_bw_na <- function(AQ_norm_list, max_na_pct = 10) {
  
  res <- imap_dfr(AQ_norm_list, function(sensor, nm) {
    t <- sensor$normalised$t
    y <- sensor$normalised$mean_value
    
    # bandwidth estimated on complete cases only
    cc <- which(!is.na(t) & !is.na(y))
    bw <- tryCatch(
      if (length(cc) >= 20) {
        rdrobust::rdbwselect(y = y[cc], x = t[cc])$bws[1, 1]
      } else NA_real_,
      error = function(e) NA_real_
    )
    
    # compute %NA within ±bw using original series
    if (is.na(bw)) {
      tibble(
        sensor = nm,
        bandwidth = NA_real_,
        total_in_bw = NA_integer_,
        na_in_bw = NA_integer_,
        pct_na_in_bw = NA_real_
      )
    } else {
      df <- tibble(t = t, y = y) |> filter(!is.na(t), abs(t) <= bw)
      total_in_bw <- nrow(df)
      na_in_bw    <- sum(is.na(df$y))
      pct_na_in_bw <- if (total_in_bw > 0) 100 * na_in_bw / total_in_bw else NA_real_
      
      tibble(
        sensor = nm,
        bandwidth = bw,
        total_in_bw = total_in_bw,
        na_in_bw = na_in_bw,
        pct_na_in_bw = pct_na_in_bw
      )
    }
  })
  
  keep_ids <- res |>
    filter(!is.na(bandwidth), !is.na(pct_na_in_bw), total_in_bw > 0, pct_na_in_bw <= max_na_pct) |>
    pull(sensor)
  
  clean_list <- AQ_norm_list[keep_ids]
  
  list(clean_list = clean_list, summary = res)
}

AQ_norm_list_clean <- screen_by_bw_na(AQ_norm_list, max_na_pct = 5)

aq_bw_all <- map(AQ_norm_list, function(sensor) {
  x <- rdbwselect(y = sensor$normalised$mean_value, x = sensor$normalised$t)
  x$bws[1,1]
}) |> as_tibble() |> pivot_longer(cols = everything(), names_to = "sensor", values_to = "bandwidth") |>
  rename(sensor_pollutant = sensor) |>
  mutate(
    sensor    = str_replace(sensor_pollutant, "_(NO2|PM2?5)$", ""),
    pollutant = str_extract(sensor_pollutant, "(NO2|PM2?5)")
  )


# #function to calculate MDE for all sensors
# RDD_mde <- map(AQ_norm_list, function(sensor) {
#   x <- rdmde(data = cbind(sensor$normalised$mean_value, sensor$normalised$t),
#               h = 26,
#               cutoff = 0,
#               alpha = 0.05,
#               beta = 0.8
#   )
#   x$mde
# })
# 
# test <- rdmde(data = cbind(GH3_PM25$normalised$median_value, GH3_PM25$normalised$t), 
#               p = 1, alpha = 0.05)
# 
# 

#function to run core RDD analysis on AQ sensors and plot graps

RDD_AQ_fn <- function(sensor, h = NULL, donut_hole = 0, p = 1) {
  sensor_name <- deparse(substitute(sensor))
  df <- sensor$normalised
  
  # Initialize output containers
  sensor$RDD <- list()
  sensor$RDD_donut <- list()
  
  ##Standard RD -----------------------------------------------
  
  bw_select <- rdbwselect(y = df$mean_value, x = df$t)
  sensor$RDD$bw_select <- bw_select
  
  # Estimation
  est_args <- list(y = df$mean_value, x = df$t, p = p, all = TRUE, kernel = "triangular")
  est_args$h <- ifelse(!is.null(h), h, bw_select$bws[1, ])
  est <- do.call(rdrobust, est_args)
  sensor$RDD$est <- est
  sensor$RDD$est_summary <- summary(est)
  
  # Fixed-bandwidth diagnostic
  h_fixed <- ifelse(!is.null(h), h, bw_select$bws[1, ])
  
  sensor$RDD$rdplot_bw <- rdplot(
    y       = df$mean_value,
    x       = df$t,
    subset  = df$t >= -h_fixed & df$t <= h_fixed,
    h       = h_fixed,
    p       = p,
    nbins   = 2 * h_fixed,
    kernel  = "uniform",           # confidence interval
    title   = paste(sensor_name, "Sharp RDD"),
    x.label = "Days from introduction",
    y.label = "mean daily concentration"
  ) 
  
  rd_bw_plot <- sensor$RDD$rdplot_bw$rdplot 
  
  #Donut RD
  if (donut_hole > 0) {
    df_donut <- subset(df, abs(t) > donut_hole)
    
    #bw_d <- rdbwselect(y = df_donut$mean_value, x = df_donut$t)
    #sensor$RDD_donut$bw_select <- bw_d
    
    h_sharp <- ifelse(!is.null(h), as.numeric(h), as.numeric(bw_select$bws[1, ]))
    h_d     <- h_sharp + donut_hole  # ensures same # points per side as sharp RDD
    
    est_args_d <- list(
      y = df_donut$mean_value,
      x = df_donut$t,
      p = p,
      all = TRUE,
      kernel = "uniform",
      h = h_d,
      bwcheck = h_sharp# CHANGED
    )
    
    est_d <- do.call(rdrobust, est_args_d)
    sensor$RDD_donut$est <- est_d
    sensor$RDD_donut$est_summary <- summary(est_d)
    
    h_fixed_d <- h_d             # CHANGED
    
    sensor$RDD_donut$rdplot_bw <- rdplot(
      y       = df_donut$mean_value,
      x       = df_donut$t,
      h       = h_fixed_d,       # CHANGED
      subset  = df_donut$t >= -h_fixed_d & df_donut$t <= h_fixed_d, # CHANGED
      p       = p,
      nbins   = 2 * h_fixed_d,
      kernel  = "uniform",
      title   = paste(sensor_name, "Donut data inside bandwidth"),
      x.label = "Days from CAZ introduction",
      y.label = "mean daily concentration (ppm)"
    ) 
  
  
  rd_donut_bw_plot <- sensor$RDD_donut$rdplot_bw$rdplot 
  
  sensor$RDD_donut$rdplot_bw$rdplot <- rd_donut_bw_plot +
    geom_rect(
      aes(xmin = -donut_hole, xmax = donut_hole,
          ymin = -Inf,        ymax = Inf),
      fill         = "orange",
      alpha        = 0.3,
      color        = 'grey60',
      linetype     = "dashed",
      inherit.aes  = FALSE
    ) +
    labs(
      title = paste(sensor_name, "RDD donut estimation"),
      x     = "Days from CAZ introduction",
      y     = "Mean daily concentration (ppm)"
    ) 
  }
  # return the object and print the bandwidth plots
  print(sensor$RDD$rdplot_bw$rdplot)
  if (donut_hole > 0) {
    print(sensor$RDD_donut$rdplot_bw$rdplot)
  }
  return(sensor)
}


# Run RDD analysis for each sensor with and without donut holes
GH4_NO2_RDD <- RDD_AQ_fn(GH4_NO2, donut_hole = 0)
DFR1027_NO2_RDD <- RDD_AQ_fn(DFR1027_NO2, donut_hole = 0)
DFR1027_PM25_RDD <- RDD_AQ_fn(DFR1027_PM25, donut_hole = 0)
DFR1063_NO2_RDD <- RDD_AQ_fn(DFR1063_NO2, donut_hole = 0)
DFR1063_PM25_RDD <- RDD_AQ_fn(DFR1063_PM25, donut_hole = 0)
GH6_NO2_RDD <- RDD_AQ_fn(GH6_NO2, donut_hole = 0)
GH6_PM25_RDD <- RDD_AQ_fn(GH6_PM25, donut_hole = 0)
GH3_NO2_RDD <- RDD_AQ_fn(GH3_NO2, donut_hole = 0)
GH3_PM25_RDD <- RDD_AQ_fn(GH3_PM25, donut_hole = 0)
AMF_PM25_RDD <- RDD_AQ_fn(AMF245_PM25, donut_hole = 0)


#pack into a list
AQ_RDD_list <- list(GH4_NO2_RDD = GH4_NO2_RDD,
                    DFR1027_NO2_RDD = DFR1027_NO2_RDD,
                    DFR1027_PM25_RDD = DFR1027_PM25_RDD,
                    DFR1063_NO2_RDD = DFR1063_NO2_RDD,
                    DFR1063_PM25_RDD = DFR1063_PM25_RDD,
                    GH6_NO2_RDD = GH6_NO2_RDD,
                    GH6_PM25_RDD = GH6_PM25_RDD,
                    GH3_NO2_RDD = GH3_NO2_RDD,
                    GH3_PM25_RDD = GH3_PM25_RDD,
                    AMF_PM25_RDD = AMF_PM25_RDD)

#function that pulls out the key coefs for each sensor-pollutant pair into one table
extract_coefs <- function(list) {
  map_dfr(names(list), function(name) {
    est_obj <- list[[name]]$RDD$est
      est_val     <- est_obj$coef[3]
      ci_lower    <- est_obj$ci[3, 1]
      ci_upper    <- est_obj$ci[3, 2]
      bw_val      <- est_obj$bws[1, 1]
      p_val       <- est_obj$pv[3]
      se          <- est_obj$se[3]
      
      tibble(
        id        = name,
        estimator = est_val,
        ci_lower  = ci_lower,
        ci_upper  = ci_upper,
        bandwidth = bw_val,
        p_value   = p_val,
        sig_95    = p_val < 0.05,
        se        = se
      )
    }
  )
}

#extract the coefs

AQ_RDD_summary <- extract_coefs(AQ_RDD_list) |>
  mutate(pollutant = str_extract(id, "NO2|PM25"),
         sensor = str_extract(id, "GH4|GH3|GH6|DFR1027|DFR1063|AMF"),
         type = case_when(
           sensor == 'GH4' | sensor == 'GH6' | sensor == 'DFR1063' | sensor == 'AMF245'  ~ 'CAZ',
           TRUE ~ 'Spillover'
         )) |>
  arrange(desc(type), desc(estimator)) |>
  mutate(sensor = factor(sensor, levels = unique(sensor))) |>
  rename(coef = estimator)

#forest plot of the estimators

ggplot(AQ_RDD_summary, aes(y = sensor, x = estimator, colour = type, shape = sig_95)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), 
                width = 0.2, size = 0.8) +  # width controls cap length
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ pollutant, scales = "free_y") +
  scale_colour_manual(labels = c('CAZ', 'Spillover'), 
                      values = c("CAZ" = "darkred", "Spillover" = "darkblue")) +
  labs(
    x = "RDD Estimate (Bias-corrected)",
    y = "Sensor",
    colour = "Sensor Type"
  ) +
  theme_ipsum_rc() +
  labs(title = 'Sharp RDD estimates for AQ sensors',
       subtitle = 'Bars show 95% confidence intervals') +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 10),
  ) 


#sensitivity analysis - produce output table of robust coefficients, p values and CIs for donut holes of 4, 5, 6 weeks respectively

RDD_donut_sensitivity <- function(sensors, donut_weeks = c(4, 5, 6), h = NULL) {
  
  if (is.null(names(sensors)) || any(names(sensors) == "")) {
    call_sensors <- match.call(expand.dots = FALSE)$sensors
    if (is.call(call_sensors) && call_sensors[[1]] == "list") {
      inferred_names <- sapply(call_sensors[-1], deparse)
      names(sensors) <- inferred_names
    }
  }
  donut_days = donut_weeks * 7
  metrics <- c("coef", "p_value", "significant", "ci_lower", "ci_upper", "se")
  
  sensor_results <- lapply(sensors, function(sensor) {
    df <- sensor$normalised
    bw_select <- rdbwselect(y = df$mean_value, x = df$t)
    
    res_vec <- unlist(lapply(donut_days, function(donut_hole) {
      df_donut <- subset(df, abs(t) > donut_hole)
      est_args <- list(y = df_donut$mean_value, x = df_donut$t, p = 1,
                       all = TRUE, kernel = "uniform")
      est_args$h <- ifelse(!is.null(h), h + donut_hole,
                           bw_select$bws[1, ] + donut_hole)
      est_d <- do.call(rdrobust, est_args)
      coef <- est_d$coef["Robust", "Coeff"]
      pval <- est_d$pv["Robust", "P>|z|"]
      ci_low <- est_d$ci["Robust", "CI Lower"]
      ci_high <- est_d$ci["Robust", "CI Upper"]
      se <- est_d$se["Robust", "Std. Err."]
      significant <- pval < 0.05
      c(coef, pval, significant, ci_low, ci_high, se)
    }))
    
    names(res_vec) <- unlist(lapply(donut_weeks, function(w)
      paste(metrics, paste0(w, "w"), sep = "_")))
    res_vec
  })
  
  result_df <- data.frame(sensor_results, check.names = FALSE)
  rownames(result_df) <- names(sensor_results[[1]])
  colnames(result_df) <- names(sensor_results)
  return(result_df)
}

#pull out coefficients and reshape for plotting

donut_sensitivity_aq <- RDD_donut_sensitivity(AQ_norm_list) |>
                        rownames_to_column("metric") |>
                        pivot_longer(-metric, names_to = "sensor_pollutant", values_to = "value") |>
                        pivot_wider(names_from = metric, values_from = value) |>
                        mutate(pollutant = str_extract(sensor_pollutant, "NO2|PM25"),
                               sensor = str_extract(sensor_pollutant, "GH4|GH3|GH6|DFR1027|DFR1063|AMF"),
                               type = case_when(
                                 sensor == 'GH4' | sensor == 'GH6' | sensor == 'DFR1063' | sensor == 'AMF245'  ~ 'CAZ',
                                 TRUE ~ 'Spillover'
                               )) |>
                        pivot_longer(
                          cols = starts_with("coef_4w"):starts_with("se_6w"),
                          names_to = c(".value", "donut"),
                          names_pattern = "(coef|p_value|significant|ci_lower|ci_upper|se)_(4w|5w|6w)"
                          ) |>
                        mutate(sensor = factor(sensor, levels = unique(sensor)),
                               donut = factor(donut, levels = c("4w", "5w", "6w")))
                          

#produce forest plot for donuts
shape_map <- c('ns' = 20, "p < 0.05" = 8)

forest_plot_donut_aq <- function(df, pol) {
  df |>
    filter(pollutant == pol) |>
    mutate(sig_lab = ifelse(significant == 1, "p < 0.05", "ns"),
           donut   = factor(donut, levels = c("4w","5w","6w"))) |>
    ggplot(aes(x = coef, y = sensor, colour = type, shape = sig_lab)) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, linewidth = 1) +
    geom_point(size = 2, stroke = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ donut, nrow = 1, scales = 'free_x') +
    scale_colour_viridis_d() +
    scale_shape_manual(values = shape_map) +
    labs(title = paste("RDD Donut Estimates –", pol),
         x = "RD Estimate (robust) with 95% CI",
         y = "Sensor", colour = "Area Type", shape = "Significance") +
    theme_ipsum_rc() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

forest_plot_donut_aq(donut_sensitivity_aq, "NO2")
forest_plot_donut_aq(donut_sensitivity_aq, "PM25")





# RDD estimators for traffic data -----------------------------------------

#aggregate data to be daily sum of cars rather than hourly, add running variable and impute
daily_avg_traffic <- function(road) {
  road$normalised_daily <- road$normalised |>
    mutate(day = as_date(date)) |>
    group_by(day) |>
    summarise(cars_per_day = mean(value_predict, na.rm = TRUE),
              .groups = "drop") |>
    complete(day = seq(min(day), max(day), by = "day")) |>
    mutate(
      cars_per_day   = na.approx(cars_per_day, x = day, na.rm = FALSE),
      t = as.numeric(day - CAZ_start)
    )       
  road
}

#apply daily average to all roads
traffic_norm_list <- map(traffic_norm_list, daily_avg_traffic)

#calculate optimal bandwidth across all sensors and take the median
tf_bw_all <- map(traffic_norm_list, function(road) {
  x <- rdbwselect(y = road$normalised_daily$cars_per_day, x = road$normalised_daily$t)
  x$bws[1,1]
})


#function to run RDD/donut analysis on traffic sensors

RDD_traffic_fn <- function(road, h = NULL, donut_hole = 0) {
  road_name <- deparse(substitute(road))
  df <- road$normalised_daily
  
  # Initialize output containers
  road$RDD <- list()
  road$RDD_donut <- list()
  
  ##Standard RD -----------------------------------------------
  
  bw_select <- rdbwselect(y = df$cars_per_day, x = df$t)
  road$RDD$bw_select <- bw_select
  
  # Estimation
  est_args <- list(y = df$cars_per_day, x = df$t, p = 1, all = TRUE, kernel = "uniform")
  est_args$h <- ifelse(!is.null(h), h, bw_select$bws[1, ])
  est <- do.call(rdrobust, est_args)
  road$RDD$est <- est
  road$RDD$est_summary <- summary(est)
  
  # Fixed-bandwidth diagnostic
  h_fixed <- ifelse(!is.null(h), h, bw_select$bws[1, ])
  
  road$RDD$rdplot_bw <- rdplot(
    y       = df$cars_per_day,
    x       = df$t,
    subset  = df$t >= -h_fixed & df$t <= h_fixed,
    h       = h_fixed,
    p       = 1,
    nbins   = 2 * h_fixed,
    kernel  = "uniform",
    title   = paste(road_name, "Sharp RDD"),
    x.label = "Days from introduction",
    y.label = "Cars per day"
  )
  rd_bw_plot <- road$RDD$rdplot_bw$rdplot
  
  #Donut RD
  if (donut_hole > 0) {
    df_donut <- subset(df, abs(t) > donut_hole)
  
    bw_d <- rdbwselect(y = df_donut$cars_per_day, x = df_donut$t)
    road$RDD_donut$bw_select <- bw_d
    
    est_args_d <- list(y = df_donut$cars_per_day, x = df_donut$t, p = 1, all = TRUE, kernel = "uniform")
    est_args_d$h <- ifelse(!is.null(h), h + donut_hole, bw_select$bws[1, ] + donut_hole)
    
    est_d <- do.call(rdrobust, est_args_d)
    road$RDD_donut$est <- est_d
    road$RDD_donut$est_summary <- summary(est_d)
    
    h_fixed_d <- ifelse(!is.null(h), est_d$bws[1, ] + donut_hole, bw_d$bws[1, ] +donut_hole)
    
    road$RDD_donut$rdplot_bw <- rdplot(
      y       = df_donut$cars_per_day,
      x       = df_donut$t,
      h       = est_args_d$h,
      subset  = df_donut$t >= -est_args_d$h & df_donut$t <= est_args_d$h,
      p       = 1,
      nbins   = 2 * est_args_d$h,
      kernel  = "uniform",
      ci      = 0.95,           
      shade   = TRUE,
      title   = paste(road_name, "Donut data inside bandwidth"),
      x.label = "Days from CAZ introduction",
      y.label = "Cars per day"
    ) 
  
  rd_donut_bw_plot <- road$RDD_donut$rdplot_bw$rdplot
  road$RDD_donut$rdplot_bw$rdplot <- rd_donut_bw_plot +
    geom_rect(
      aes(xmin = -donut_hole, xmax = donut_hole,
          ymin = -Inf,        ymax = Inf),
      fill         = "lightgreen",
      alpha        = 0.3,
      color        = 'grey60',
      linetype     = "dashed",
      inherit.aes  = FALSE
    ) +
    labs(
      title = paste(road_name, "RDD donut estimation"),
      x     = "Days from CAZ introduction",
      y     = "Cars per day"
    )
  }
  return(road)
}

# Run RDD analysis for each road without donut holes with map

traffic_RDD_list <- map(traffic_norm_list, function(road) {
  RDD_traffic_fn(road, donut_hole = 0)
})

#make a dataframe of robust coefficients and CIs for all traffic sensors, join to metadata for inside-CAZ/road type info
traffic_RDD_df <- map_dfr(traffic_RDD_list, function(road) {
  est <- road$RDD$est
  tibble(
    coef     = est$coef["Robust","Coeff"],
    p_value  = est$pv["Robust","P>|z|"],
    ci_lower = est$ci["Robust","CI Lower"],
    ci_upper = est$ci["Robust","CI Upper"],
    se = est$se["Robust", 'Std. Err.']
  )
}, .id = "sensor") |>
  left_join(sensor_to_road_RDD, by = c("sensor" = "ref")) |> 
  distinct(sensor, .keep_all = TRUE) |>
  select(-c(sensorID, road_id, match_method, match_dist_m, name))

#forest plot
ggplot(traffic_RDD_df, aes(x = coef, y = sensor)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), 
                width = 0.2, size = 0.8) +  # width controls cap length
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(
    x = "RDD Estimate (Bias-corrected)",
    y = "Sensor"
  ) +
  theme_ipsum_rc() +
  labs(title = 'Sharp RDD estimates for traffic sensors',
       subtitle = 'Bars show 95% confidence intervals') +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 10),
  )

# Function to run RDD analysis for each road with donut holes of 4, 5, 6 weeks respectively
RDD_donut_sensitivity_tf <- function(sensors, donut_weeks = c(4, 5, 6), h = NULL) {
  # Infer sensor names from the call if not provided
  if (is.null(names(sensors)) || any(names(sensors) == "")) {
    call_sensors <- match.call(expand.dots = FALSE)$sensors
    if (is.call(call_sensors) && call_sensors[[1]] == "list") {
      inferred_names <- sapply(call_sensors[-1], deparse)
      names(sensors) <- inferred_names
    }
  }
  
  # convert weeks to days for hole sizes
  donut_days <- donut_weeks * 7
  metrics <- c("coef", "p_value", "significant", "ci_lower", "ci_upper","se")
  
  sensor_results <- lapply(sensors, function(sensor) {
    df <- sensor$normalised_daily
    # bandwidth from full sample (as in original function)
    bw_select <- rdbwselect(y = df$cars_per_day, x = df$t)
    
    res_vec <- unlist(lapply(donut_days, function(donut_hole) {
      df_donut <- subset(df, abs(t) > donut_hole)
      est_args <- list(y = df_donut$cars_per_day, x = df_donut$t, p = 1,
                       all = TRUE, kernel = "uniform")
      est_args$h <- ifelse(!is.null(h), h + donut_hole,
                           bw_select$bws[1, ] + donut_hole)
      est_d <- do.call(rdrobust, est_args)
      coef <- est_d$coef["Robust", "Coeff"]
      pval <- est_d$pv["Robust", "P>|z|"]
      ci_low <- est_d$ci["Robust", "CI Lower"]
      ci_high <- est_d$ci["Robust", "CI Upper"]
      significant <- pval < 0.05
      se <- est_d$se["Robust", "Std. Err."]
      c(coef, pval, significant, ci_low, ci_high, se)
    }))
    
    names(res_vec) <- unlist(lapply(donut_weeks, function(w)
      paste(metrics, paste0(w, "w"), sep = "_")))
    res_vec
  })
  
  result_df <- data.frame(sensor_results, check.names = FALSE)
  rownames(result_df) <- names(sensor_results[[1]])
  colnames(result_df) <- names(sensor_results)
  return(result_df)
}


#apply to the roads list
tf_donut_sensitivity <- RDD_donut_sensitivity_tf(traffic_RDD_list) |>
  rownames_to_column("metric") |>
  pivot_longer(-metric, names_to = "sensor", values_to = "value") |>
  pivot_wider(names_from = metric, values_from = value) |>
  mutate(pollutant = "Cars",
         type = "Traffic") |>
  pivot_longer(
    cols = starts_with("coef_4w"):starts_with("se_6w"),
    names_to = c(".value", "donut"),
    names_pattern = "(coef|p_value|significant|ci_lower|ci_upper|se)_(4w|5w|6w)"
  ) |>
  mutate(sensor = factor(sensor, levels = unique(sensor)),
         donut = factor(donut, levels = c("4w", "5w", "6w"))) |>
  mutate(type = str_extract(sensor, "^A|^B")) |>
  mutate(road_type = case_when(type == 'A' ~ 'Primary road',
                               type == 'B' ~ 'Secondary road',
                               is.na(type) ~ 'CAZ Road',
                               type == 'A61 Ring Road' ~ 'CAZ Road'),
         category = case_when(road_type == 'Primary road' | road_type == 'Secondary road' ~ 'Spillover',
                              road_type == 'CAZ Road' ~ 'CAZ') 
         )|>
  arrange(desc(road_type), desc(coef)) |>
  mutate(sensor = factor(sensor, levels = unique(sensor))) 
#forest plot for traffic donut

forest_plot_donut_tf <- function(df) {
  df |>
    mutate(sig_lab = ifelse(significant == 1, "p < 0.05", "ns"),
           donut   = factor(donut, levels = c("4w","5w","6w"))) |>
    ggplot(aes(x = coef, y = sensor, colour = road_type, shape = sig_lab)) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, linewidth = 1) +
    geom_point(size = 2, stroke = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ donut, nrow = 1, scales = 'free_x') +
    scale_colour_viridis_d() +
    scale_shape_manual(values = shape_map) +
    labs(x = "RD Estimate (robust) with 95% CI",
         y = "Road", colour = "Road Type", shape = "Significance") +
    theme_ipsum_rc() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

forest_plot_donut_tf(traffic_donut_sensitivity_df)



# tables of aggregated coefficients: AQ ---------------------------------------

pool_simple <- function(df){               
  n   <- nrow(df)
  tau <- mean(df$coef)
  se  <- sqrt(sum(df$se^2)) / n    
  z   <- tau / se
  p   <- 2 * (1 - pnorm(abs(z)))
  tibble(est_mean = tau,
         se_mean  = se,
         lo_mean  = tau - 1.96*se,
         hi_mean  = tau + 1.96*se,
         z_mean = z,
         p_mean   = p)
}

pool_donut <- function(df){               
  n   <- nrow(df)
  tau <- mean(df$coef)
  se  <- sqrt(sum(df$se^2)) / n    
  z   <- tau / se
  p   <- 2 * (1 - pnorm(abs(z)))
  tibble(est_mean = tau,
         se_mean  = se,
         lo_mean  = tau - 1.96*se,
         hi_mean  = tau + 1.96*se,
         z_mean = z,
         p_mean   = p)
}


#for Sharp AQ RDD
# overall means by pollutant 
aq_pooled <-
  AQ_RDD_summary |>
  group_by(pollutant) |>
  group_modify(~ pool_simple(.x)) |>
  ungroup()

# within CAZ / Spillover by pollutant
aq_CAZvsNot_pooled <-
  AQ_RDD_summary |>
  group_by(pollutant, type) |>
  group_modify(~ pool_simple(.x)) |>
  ungroup()

aq_all_pooled <- bind_rows(aq_pooled, aq_CAZvsNot_pooled)
rm(aq_pooled, aq_CAZvsNot_pooled)

#for donut RDD
# overall means by pollutant and donut hole size
aq_donut_pooled <-
  donut_sensitivity_aq |>
  group_by(pollutant, donut) |>
  group_modify(~ pool_donut(.x)) |>
  ungroup()

aq_donut_CAZvNot_pooled <-
  donut_sensitivity_aq |>
  group_by(pollutant, donut, type) |>
  group_modify(~ pool_donut(.x)) |>
  ungroup()

aq_donut_all_pooled <- bind_rows(aq_donut_pooled, aq_donut_CAZvNot_pooled)
rm(aq_donut_pooled, aq_donut_CAZvNot_pooled)

# tables of aggregated coefficients: traffic ------------------------------

#for Sharp traffic RDD
tf_pooled <-
  traffic_RDD_df |>
  pool_simple() |>
  mutate(category = "All traffic sensors")

#by inside_caz/outside caz
tf_CAZvNot_pooled <-
  traffic_RDD_df |>
  group_by(category) |>
  group_modify(~ pool_simple(.x))

#by primary/secondary road
tf_primVsec_pooled <-
  traffic_RDD_df |>
  group_by(highway) |>
  group_modify(~ pool_simple(.x))

tf_all_pooled <- bind_rows(tf_pooled, tf_CAZvNot_pooled, tf_primVsec_pooled)
rm(tf_pooled, tf_CAZvNot_pooled, tf_primVsec_pooled)


#for donut RDD

tf_donut_pooled <-
  tf_donut_sensitivity |>
  pool_donut()

tf_donut_CAZvNot_pooled <-
  tf_donut_sensitivity |>
  group_by(category, donut) |>
  group_modify(~ pool_donut(.x)) |>
  ungroup()

tf_donut_primVsec_pooled <-
  tf_donut_sensitivity |>
  group_by(road_type, donut) |>
  group_modify(~ pool_donut(.x)) |>
  ungroup()

tf_donut_all_pooled <- bind_rows(tf_donut_pooled, tf_donut_CAZvNot_pooled, tf_donut_primVsec_pooled)
rm(tf_donut_pooled, tf_donut_CAZvNot_pooled, tf_donut_primVsec_pooled)

