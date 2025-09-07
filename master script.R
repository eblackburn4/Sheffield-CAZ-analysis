
## ---------------------------
## Purpose of script: master script to run other scripts + load packages 
## Author: Ned Blackburn
## Date Created: 2025-08-02

options(scipen = 6, digits = 5) 
library(ggmap)
library(ggpattern)
library(ggrepel)
library(ggsflabel)
library(ggthemes)
library(hrbrthemes)
library(jsonlite)
library(osmdata)
library(patchwork)
library(ranger)
library(rdpower)
library(rdrobust)
library(rmweather)
library(scico)
library(sf)
library(tidyverse)
library(weathermetrics)
library(zoo)

## ---------------------------

source("Mapping code.R")
source("data exploration.R")
source('weather normalisation.R')
source('RDD.R')