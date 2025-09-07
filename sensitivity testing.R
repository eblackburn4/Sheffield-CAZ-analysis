
## ---------------------------
## Purpose of script: sensitivity analysis for sharp RDD: traffic and AQ
## Author: Ned Blackburn
## Date Created: 2025-08-19

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)
library(rdrobust)

## ---------------------------

#Bandwidth test: air quality
#function to run RDD and extract estimates
get_rd_est <- function(y, x, h) {
  
  fit <- rdrobust(y = y, x = x, h = h, kernel = "uniform")
  
  # safer extraction by names
  est   <- unname(fit$coef["Robust", 1])
  se    <- unname(fit$se["Robust", 1])
  ci_lo <- unname(fit$ci["Robust", 1])
  ci_hi <- unname(fit$ci["Robust", 2])
  pv <- unname(fit$pv["Robust", 1])
  
  tibble(
    estimate = est,
    conf.low = ci_lo,
    conf.high = ci_hi,
    se = se,
    pvalue = pv
  )
}

#run bw_tests for various values of h
run_bw_test_aq <- function(AQ_norm_list, h_grid = c(7, 14, 21, 28, 35)) {
  map_dfr(names(AQ_norm_list), function(nm) {
    s <- AQ_norm_list[[nm]]
    y <- s$normalised$mean_value
    x <- s$normalised$t
    
    map_dfr(h_grid, function(h) {
      get_rd_est(y, x, h) |>
        mutate(sensor_pollutant = nm, h = h, .before = 1) |>
        mutate(significant = ifelse(pvalue < 0.05, 1, 0),
               sensor_pollutant = str_replace(sensor_pollutant, "_", " "))
    })
  })
}

AQ_bw_test_summary <- run_bw_test_aq(AQ_norm_list) |>
  mutate(
    sensor    = str_extract(sensor_pollutant, "(GH4|GH3|GH6|DFR1027|DFR1063|AMF)"),
    pollutant = str_extract(sensor_pollutant, "(NO2|PM25)")
  ) |>
  mutate(sensor = if_else(sensor == 'AMF', 'AMF245', sensor)) |>
  left_join(
    aq_EDA_stats |> select(sensor, pollutant, pre_mean),
    by = c("sensor", "pollutant")
  ) |>
  mutate(pct_change = 100 * estimate / pre_mean)
  
#bandwidth test: traffic
run_bw_test_tf <- function(traffic_norm_list, h_grid = c(7, 14, 21, 28, 35)) {
  map_dfr(names(traffic_norm_list), function(nm) {
    s <- traffic_norm_list[[nm]]
    y <- s$normalised_daily$cars_per_day
    x <- s$normalised_daily$t
    
    map_dfr(h_grid, function(h) {
      get_rd_est(y, x, h) |>
        mutate(road = nm, h = h, .before = 1) |>
        mutate(significant = ifelse(pvalue < 0.05, 1, 0))
    })
  })
}

#generate summary tablea and join to EDA table for pre-CAZ percentages
TF_bw_test_summary <- run_bw_test_tf(traffic_norm_list) |>
  rename(sensor = road) |>
  left_join(tf_EDA_stats |> select(sensor, pre_mean), by = "sensor") |>
  mutate(pct_change = 100 * estimate / pre_mean) 

#plot heatmaps summarising the results (per pollutant)

plot_heatmap <- function(df_pollutant, title_lab = "") {
  # order sensors
  sensor_order <- df_pollutant %>%
    group_by(sensor) %>%
    summarise(mn = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mn)) %>% pull(sensor)
  
  # factors
  dfp <- df_pollutant %>%
    mutate(
      sensor = factor(sensor, levels = sensor_order),
      h      = factor(h, levels = sort(unique(h)))
    )
  
  ggplot(dfp, aes(x = h, y = sensor, fill = pct_change)) +
    geom_tile() +
    geom_text(
      data = dfp %>% filter(significant == 1),
      aes(label = "*"),
      size = 12,
      color = 'grey80',
      vjust = 0.8
    ) +
    scale_fill_scico(palette = 'vik', midpoint = 0, direction = -1) +
    labs(title = title_lab, x = "Bandwidth (days)", y = "Sensor", fill = "") +
    theme_ipsum_rc(axis_title_size = 10, axis_text_size = 9) +
    theme(panel.grid = element_blank(), legend.position = "right")
}

# split and plot with y axis labels to match sensor names

plot_heatmap(AQ_bw_test_summary |> filter(pollutant == 'NO2')) +
   scale_y_discrete(limits = c('GH3', 'DFR1063','GH6','DFR1027','GH4'),
                    labels = c(
                    "GH3"    = "SCC_GH3",
                    "GH6"    = "SCC_GH6",
                    "DFR1027"= "DFR_1027A",
                    "DFR1063"= "DFR_1063A",
                    "GH4"    = "SCC_GH4"))
       
  
plot_heatmap(AQ_bw_test_summary |> filter(pollutant == 'PM25')) +
  scale_y_discrete(limits = c('GH3', 'DFR1063','GH6','DFR1027','AMF245'),
                   labels = c(
                     "GH3"    = "SCC_GH3",
                     "GH6"    = "SCC_GH6",
                     "DFR1027"= "DFR_1027A",
                     "DFR1063"= "DFR_1063A",
                     "AMF245"    = "AMF_2450229"))

#heatmap for traffic
plot_heatmap(TF_bw_test_summary) +
  scale_y_discrete(limits = rev(sensor_order_TF),
                   labels = rev(sensor_order_TF)) +
  labs(y = 'Road')



#calculate MDE for each RDD model

calculate_mde <- function(sensor_list, power = 0.8, alpha = 0.05) {
  map_dfr(names(sensor_list), function(name) {
    sensor <- sensor_list[[name]]
    df     <- sensor$normalised
    h      <- sensor$RDD$est$bws[1, 1]
    data   <- cbind(df$mean_value, df$t)
    samph  <- c(h, h)
    
    power_diff <- function(tau) {
      rdpower(data = data, tau = tau, samph = samph, alpha = alpha)$power.rbc - power
    }
    
    upper <- sd(df$mean_value, na.rm = TRUE)
    while (power_diff(upper) < 0) {
      upper <- upper * 2
    }
    
    tau_mde  <- uniroot(power_diff, lower = 0, upper = upper)$root
    pre_mean <- mean(df$mean_value[df$t < 0], na.rm = TRUE)
    
    tibble(
      id        = name,
      pollutant = str_extract(name, "NO2|PM25"),
      sensor    = str_extract(name, "GH4|GH3|GH6|DFR1027|DFR1063|AMF"),
      mde       = tau_mde,
      mde_pct   = 100 * tau_mde / pre_mean,
      bw = samph[1]
    )
  })
}

AQ_MDE_table <- calculate_mde(AQ_RDD_list) |>
  mutate(sensor = factor(sensor, levels = c("GH4", "DFR1027", "GH6", "DFR1063", 'GH3', "AMF"))) |>
  arrange(pollutant, sensor)


