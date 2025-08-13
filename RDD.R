
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


# RDD estimators for AQ data ----------------------------------------------
#aggregate data to be daily avg rather than hourly and add running variable

CAZ_start <- as_date("2023-02-27")

daily_avg <- function(sensor) {
  sensor$normalised <- sensor$normalised |>
    mutate(day = as_date(date)) |>
    group_by(day) |>
    summarise(mean_value = mean(value_predict),
              median_value = median(value_predict), 
              .groups = "drop") |>
    mutate(t = as.numeric(day - CAZ_start))            
  sensor
} 

#apply daily average to all sensors

AQ_norm_list <- list(  GH4_NO2  = GH4_NO2,
                       GH4_PM25 = GH4_PM25,
                       GH3_NO2  = GH3_NO2,
                       GH3_PM25 = GH3_PM25,
                       GH6_NO2  = GH6_NO2,
                       GH6_PM25 = GH6_PM25,
                       DFR_NO2  = DFR_NO2,
                       DFR_PM25 = DFR_PM25 )

AQ_norm_list <- map(AQ_norm_list, daily_avg)
list2env(AQ_norm_list, envir = .GlobalEnv)

#example normalised time series plot
ggplot(GH6_NO2$normalised, aes(x = day)) +
  geom_line(aes(y = mean_value), color = 'blue') +
  labs(title = "GH4 NO2 Daily Mean Values",
       x = "Date",
       y = "Mean NO2 Concentration") +
  theme_minimal()

#functiont to plot normalised time series 5 weeks either side of bandwidth
norm_plot <- function(sensor, h = 35) {
  sensor_name <- deparse(substitute(sensor))
  
  df <- sensor$normalised |>
        filter(t >= -h & t <= h) 
        # 
        # ggplot(aes(x = t, y = mean_value)) +
        #   geom_line(color = 'blue') +
        #   labs(title = paste(sensor_name, "Normalised Time Series"),
        #        x = "Date",
        #        y = "Mean Daily Concentration") +
        #   theme_minimal()
  df
}

#calculate optimal bandwidth across all sensors and take the median
aq_bw_all <- map(AQ_norm_list, function(sensor) {
  x <- rdbwselect(y = sensor$normalised$mean_value, x = sensor$normalised$t)
  x$bws[1,1]
})

#take median optimal bandwidth for consistency across sensors
aq_median_bw <- round(median(unlist(aq_bw_all)))

#function to calculate MDE for all sensors
RDD_mde <- map(AQ_norm_list, function(sensor) {
  x <- rdmde(data = cbind(sensor$normalised$mean_value, sensor$normalised$t),
              h = 26,
              cutoff = 0,
              alpha = 0.05,
              beta = 0.8
  )
  x$mde
})


#count NA values across all normalised sensor time series



test <- rdmde(data = cbind(GH6_NO2$normalised$mean_value, GH6_NO2$normalised$t), 
              p = 1, 
              h = aq_median_bw)



#function to run RDD analysis on AQ sensors
RDD_AQ_fn <- function(sensor, h = NULL, donut_hole = 0) {
  sensor_name <- deparse(substitute(sensor))
  df <- sensor$normalised
  
  # Initialize output containers
  sensor$RDD <- list()
  sensor$RDD_donut <- list()
  
  ##Standard RD -----------------------------------------------
  
  bw_select <- rdbwselect(y = df$mean_value, x = df$t)
  sensor$RDD$bw_select <- bw_select
  
  # Estimation
  est_args <- list(y = df$mean_value, x = df$t, p = 1, all = TRUE, kernel = "uniform")
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
    p       = 1,
    nbins   = 2 * h_fixed,
    kernel  = "uniform",
    title   = paste(sensor_name, "Sharp RDD"),
    x.label = "Days from introduction",
    y.label = "mean daily concentration"
  ) 
  
rd_bw_plot <- sensor$RDD$rdplot_bw$rdplot 
  
  #Donut RD
  if (donut_hole > 0) {
    df_donut <- subset(df, abs(t) > donut_hole)
  
    bw_d <- rdbwselect(y = df_donut$mean_value, x = df_donut$t)
    sensor$RDD_donut$bw_select <- bw_d
    
    est_args_d <- list(y = df_donut$mean_value, x = df_donut$t, p = 1, all = TRUE, kernel = "uniform")
    #est_args_d$h <- ifelse(!is.null(h), h + donut_hole, bw_d$bws[1, ] + donut_hole)
    est_args_d$h <- ifelse(!is.null(h), h + donut_hole, bw_select$bws[1, ] + donut_hole)
    
    est_d <- do.call(rdrobust, est_args_d)
    sensor$RDD_donut$est <- est_d
    sensor$RDD_donut$est_summary <- summary(est_d)
    
    h_fixed_d <- ifelse(!is.null(h), est_d$bws[1, ] + donut_hole, bw_d$bws[1, ] +donut_hole)
    
    sensor$RDD_donut$rdplot_bw <- rdplot(
      y       = df_donut$mean_value,
      x       = df_donut$t,
      h       = est_args_d$h,
      subset  = df_donut$t >= -est_args_d$h & df_donut$t <= est_args_d$h,
      p       = 1,
      nbins   = 2 * est_args_d$h,
      kernel  = "uniform",
      ci      = 0.95,           
      shade   = TRUE,
      title   = paste(sensor_name, "Donut data inside bandwidth"),
      x.label = "Days from CAZ introduction",
      y.label = "mean daily concentration (ppm)"
    ) 
  }
  
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

#return the object and print the bandwidth plots
  print(sensor$RDD$rdplot_bw$rdplot)
  if (donut_hole > 0) {
    print(sensor$RDD_donut$rdplot_bw$rdplot)
  }
  return(sensor)
}

# Run RDD analysis for each sensor with and without donut holes
GH4_NO2_RDD <- RDD_AQ_fn(GH4_NO2, donut_hole = 28)
GH4_PM25_RDD <- RDD_AQ_fn(GH4_PM25, donut_hole = 28)
DFR_NO2_RDD <- RDD_AQ_fn(DFR_NO2, donut_hole = 28)
DFR_PM25_RDD <- RDD_AQ_fn(DFR_PM25, donut_hole = 28)
GH6_NO2_RDD <- RDD_AQ_fn(GH6_NO2, donut_hole = 28)
GH6_PM25_RDD <- RDD_AQ_fn(GH6_PM25, donut_hole = 28)
GH3_NO2_RDD <- RDD_AQ_fn(GH3_NO2, donut_hole = 28)
GH3_PM25_RDD <- RDD_AQ_fn(GH3_PM25, donut_hole = 28)

#sensitivity analysis - produce output table of robust coefficients, p values and CIs for donut holes of 4, 5, 6 weeks respectively

RDD_donut_sensitivity <- function(sensors, donut_weeks = c(4, 5, 6), h = NULL) {
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
  metrics <- c("coef", "p_value", "significant", "ci_lower", "ci_upper")
  
  sensor_results <- lapply(sensors, function(sensor) {
    df <- sensor$normalised
    # bandwidth from full sample (as in original function)
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
      significant <- pval < 0.05
      c(coef, pval, significant, ci_low, ci_high)
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

donut_sensitivity_aq <- RDD_donut_sensitivity(list(GH4_NO2_RDD, 
                                                   GH4_PM25_RDD, 
                                                   DFR_NO2_RDD, 
                                                   DFR_PM25_RDD, 
                                                   GH6_NO2_RDD, 
                                                   GH6_PM25_RDD, 
                                                   GH3_NO2_RDD, 
                                                   GH3_PM25_RDD))



#function to calculate MDE for all sensors


