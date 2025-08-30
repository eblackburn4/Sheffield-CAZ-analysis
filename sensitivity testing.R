
## ---------------------------
## Purpose of script: sensitivity analysis for sharp RDD: traffic and AQ
## Author: Ned Blackburn
## Date Created: 2025-08-19

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)
library(rdrobust)

## ---------------------------


get_rd_est <- function(y, x, h) {
  
  fit <- rdrobust(y = y, x = x, h = h, kernel = "triangular")
  
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

run_bw_test <- function(AQ_norm_list, h_grid = c(7, 14, 21, 28, 35)) {
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


# Run the grid analysis
AQ_bw_test_summary <- run_bw_test(AQ_norm_list) |>
  mutate(
    sensor    = str_replace(sensor_pollutant, "_(NO2|PM2?5)$", ""),
    pollutant = str_extract(sensor_pollutant, "(NO2|PM2?5)")
  )
  

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
  
  ggplot(dfp, aes(x = h, y = sensor, fill = estimate)) +
    geom_tile() +
    geom_text(
      data = dfp %>% filter(significant == 1),
      aes(label = "*"),
      size = 8
    ) +
    scale_fill_gradient2(midpoint = 0, na.value = "grey85") +
    labs(title = title_lab, x = "Bandwidth h (days)", y = "Sensor", fill = "RD coef") +
    theme_ipsum_rc() +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 9))
}

# --- main ------------------------------------------------------------------

# split and plot

plot_heatmap((AQ_bw_test_summary |> filter(pollutant == 'NO2')), "RD estimates by h — NO2")
plot_heatmap((AQ_bw_test_summary |> filter(pollutant == 'PM25')), "RD estimates by h — PM25")

p_pm25 <- plot_heatmap(df_PM25, lines_PM25, "RD estimates by h — PM₂.₅")

p_no2
p_pm25


