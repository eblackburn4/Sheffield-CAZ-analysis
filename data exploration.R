
## ---------------------------
## Purpose of script: data exploration and preparation code for dissertation
## Author: Ned Blackburn
## Date Created: 2025-06-27

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)
library(rmweather)

## ---------------------------

#load in sensor uptime data

sensor_uptime <- read_delim("sensoruptime.txt", delim = '|') |>
  rename(metrics = measurements) |>
  mutate(metrics = str_remove_all(metrics,'\\\\')) |>
  mutate(NOx = ifelse(str_detect(metrics, "NO"), 1, 0),
         PMx = ifelse(str_detect(metrics, "PM"), 1, 0)) |>
  mutate(
    tBegin = as.Date(tBegin),
    tEnd   = if_else(tEnd == "Active",
                     as.Date(NA),
                     as.Date(tEnd))
  )

#plot sensor uptimes

NOx_up <- sensor_uptime |>
  filter(NOx == 1) |>
  bind_rows(
  NOx_up %>% select(date = tBegin) %>% mutate(delta =  1),
  NOx_up %>% select(date = tEnd) %>% mutate(delta = -1)
) %>%
  group_by(date) %>%
  summarise(delta = sum(delta), .groups="drop") %>%
  arrange(date) %>%
  mutate(active_sensors = cumsum(delta))

PMx_up <- sensor_uptime |>
  filter(PMx == 1) |>
  bind_rows(
    PMx_up %>% select(date = tBegin) %>% mutate(delta =  1),
    PMx_up %>% select(date = tEnd) %>% mutate(delta = -1)
  ) %>%
  group_by(date) %>%
  summarise(delta = sum(delta), .groups="drop") %>%
  arrange(date) %>%
  mutate(active_sensors = cumsum(delta))

ggplot() +
  geom_step(data = NOx_up, aes(x = date, y = active_sensors), color = 'red') +
  geom_step(data = PMx_up, aes(x = date, y = active_sensors), color = 'blue') +
  coord_cartesian(
    xlim = c(min(events$date), yesterday)
  ) +
  labs(
    x = "Date",
    y = "Number of sensors up",
    title = "Active sensors over time"
  ) +
  theme_minimal() 
