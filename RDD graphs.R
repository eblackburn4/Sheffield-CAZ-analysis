
## ---------------------------
## Purpose of script: erample RDD plot for thesis
## Author: Ned Blackburn
## Date Created: 2025-08-23

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)

## ---------------------------

library(tidyverse)
set.seed(42)

n <- 200
t <- seq(-n/2 + 1, n/2, by = 1)  # centered at 0: -199..200
w <- n/6                         # fit within ± n/4

pre  <- t[t < 0]
post <- t[t >= 0]

# denser, noisier, more distinct trends and level shift
y_pre  <- 60 + 0.018*pre  + rnorm(length(pre),  sd = 6)
y_post <- 10 - 0.035*post + rnorm(length(post), sd = 6)

df <- tibble(
  t = c(pre, post),
  period = c(rep("pre", length(pre)), rep("post", length(post))),
  y = c(y_pre, y_post)
)

fit_pre  <- df |> filter(t >= -w, t < 0)
fit_post <- df |> filter(t >= 0,  t <= w)

ggplot(df, aes(x = t, y = y)) +
  geom_line(color = "grey80", linewidth = 0.5) +
  geom_smooth(data = fit_pre,  method = "lm", se = FALSE, color = "#7B1FA2", linewidth = 1.2) +
  geom_smooth(data = fit_post, method = "lm", se = FALSE, color = "#7B1FA2", linewidth = 1.2) +
  geom_vline(xintercept = 0, linewidth = 0.8, linetype = 'dashed', colour = "grey10") +
  labs(x = "Day", y = "Pollutant concentration") +
  theme_ipsum_rc(grid = '', axis_title_size = 12.5) +
  theme(axis.line = element_line(color = "grey30", linewidth = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())


library(tidyverse)
set.seed(30)

n <- 400
t <- seq(-n/2 + 1, n/2, by = 1)  # -199..200 centered at 0
w <- n/6                         # fit length (same as base)
g <- n/6                         # half-gap width around 0

pre_t  <- t[t < 0]
post_t <- t[t >= 0]

# noisy series with different trends and a level shift
y_pre  <- 60 + 0.018*pre_t  + rnorm(length(pre_t),  sd = 6)
y_post <- 10 - 0.035*post_t + rnorm(length(post_t), sd = 6)

df <- tibble(
  t = c(pre_t, post_t),
  period = c(rep("pre", length(pre_t)), rep("post", length(post_t))),
  y = c(y_pre, y_post)
)

# truncate base series to leave a hole (-g, g); keep sides separate to avoid bridging
df_pre_trunc  <- df |> filter(t <= -g)
df_post_trunc <- df |> filter(t >=  g)

# fit windows of length w, ending at -g and starting at +g
fit_pre  <- df |> filter(t >= -g - w + 1, t <= -g)
fit_post <- df |> filter(t >=  g, t <=  g + w - 1)

ggplot() +
  geom_line(data = df_pre_trunc,  aes(t, y), color = "grey80", linewidth = 0.5) +
  geom_line(data = df_post_trunc, aes(t, y), color = "grey80", linewidth = 0.5) +
  geom_smooth(data = fit_pre,  aes(t, y), method = "lm", se = FALSE, color = "#7B1FA2", linewidth = 1.2) +
  geom_smooth(data = fit_post, aes(t, y), method = "lm", se = FALSE, color = "#7B1FA2", linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, colour = "grey10") +
  labs(x = "Day", y = "Pollutant concentration") +
  theme_ipsum_rc(grid = FALSE, axis_title_size = 12.5) +
  theme(axis.line = element_line(color = "grey30", linewidth = 0.5),
                axis.text.x = element_blank(),
                axis.text.y = element_blank())
