
## ---------------------------
## Purpose of script: code for generation of maps and graphs in diss proposal
## Author: Ned Blackburn
## Date Created: 2025-03-19

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)
library(ggthemes)
library(MASS)  # for robust regression
library(ggmap)
library(jsonlite)
library(sf)
library(chattr)

## ---------------------------

# simulated RDD chart -----------------------------------------------------


#simulate data for RDD graph
# Set seed for reproducibility of simulated noise
set.seed(123)

# Generate fake time series(600 days)
time <- seq.Date(from = as.Date("2023-01-01"), 
                 by = "day", length.out = 1000)

# Define change point
change_point <- as.Date("2024-06-20")  # Later date due to longer series

# Function to generate autocorrelated noise (AR(1) process)
generate_ar1_noise <- function(n, phi = 0.8, sigma = 0.5) {
  eps <- rnorm(n, sd = sigma) # White noise
  ar1_noise <- numeric(n)
  ar1_noise[1] <- eps[1]
  
  for (i in 2:n) {
    ar1_noise[i] <- phi * ar1_noise[i - 1] + eps[i]
  }
  return(ar1_noise)
}

# Phase 1: Nearly flat line at high y-value (longer duration)
n1 <- sum(time < change_point)
phase1 <- data.frame(
  date = time[1:n1],
  value = 80 + 0.0003 * (1:n1) + generate_ar1_noise(n1) # Small slope, high value
)

# Phase 2: Gentle downward slope, starting below Phase 1 end (longer duration)
n2 <- sum(time >= change_point)
gap_size <- 1  # Increased gap between the two phases

# Adjusted start for phase 2 to model discontinuity
phase2_time <- time[(n1 + gap_size + 1):(n1 + gap_size + n2)]
phase2 <- data.frame(
  date = phase2_time,
  value = 70 - 0.0002*(1:length(phase2_time)) + generate_ar1_noise(length(phase2_time)) # Gentle downward slope
)

#initialize robust linear model for each phase

model1 <- rlm(value ~ as.numeric(date), data = phase1)  # Phase 1 fit line
model2 <- rlm(value ~ as.numeric(date), data = phase2)  # Phase 2 fit line

# Predict fitted values
phase1$fitted <- predict(model1, newdata = phase1)
phase2$fitted <- predict(model2, newdata = phase2)

# Plot separately for each phase to prevent line joining
ggplot() +
  geom_line(data = phase1, aes(x = date, y = value), color = "grey85") +
  geom_line(data = phase2, aes(x = date, y = value), color = "grey85") +
  geom_line(data = phase1, aes(x = date, y = fitted), color = "purple", linewidth= 0.8) +
  geom_line(data = phase2, aes(x = date, y = fitted), color = "purple", linewidth = 0.8) +
  geom_vline(xintercept = change_point, linetype = "dashed", color = "grey30", linewidth = 0.9) +
  labs(x = "Time",
       y = "Pollutant concentration") +
  theme_ipsum(grid = '', axis = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12, hjust = 0.95),
        axis.title.y = element_text(size = 12, vjust = -3))



# map of sensor locations -------------------------------------------------

#make base map

sheffield <- c(left = - 1.495, bottom = 53.367, right = - 1.445, top = 53.395)
basemap <- get_stadiamap(sheffield, zoom = 15, maptype="stamen_toner_lite") |> ggmap()

#read in CAZ json

polygon_data <- fromJSON("mapfinal.geojson")

# Extract coordinates

coords <- polygon_data$features$geometry$coordinates[[1]]
polygon_coords <- coords[1, , ]

# Convert the resulting 24 x 2 matrix into a data frame.
df_polygon <- as.data.frame(polygon_coords)
colnames(df_polygon) <- c("lon", "lat")

#plot on base map

basemap + 
  geom_polygon(data = df_polygon, aes(x = lon, y = lat),
             fill = 'pink',        
             color = "red",
             alpha = 0.5,   
             size = 1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
  
