
## ---------------------------
## Purpose of script: master script to run other scripts + load packages 
## Author: Ned Blackburn
## Date Created: 2025-08-02

options(scipen = 6, digits = 5) 
library(tidyverse)
library(hrbrthemes)
library(rdrobust)
library(rmweather) #for meteorological normalisation
library(ranger) #dependency for rmweather
library(weathermetrics) #for calculating humidity from dewpoint
library(patchwork)
library(ggthemes)
library(ggmap)
library(jsonlite)
library(sf)
library(ggrepel)

## ---------------------------

#set random seed to ensure reproducibility of weather normalisation/RDD models
set.seed(9999)

source("Mapping code.R")
source("data exploration.R")
source('weather normalisation.R')
source('RDD.R')