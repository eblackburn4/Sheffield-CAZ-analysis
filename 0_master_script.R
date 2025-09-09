
## ---------------------------
## Purpose of script: master script to run other scripts + load required packages. 
## Running this script using 'source' will automatically run all the parts of the analysis in order
## Author: Ned Blackburn
## Date Created: 2025-08-02

options(scipen = 6, digits = 5) 
library(ggmap) #for ggplotting maps
library(ggrepel) #for map labels
library(ggsflabel) #for map labels
library(ggthemes) # nice ggplot themes
library(hrbrthemes) # more nice ggplot themes
library(jsonlite) #for reading in json files
library(osmdata) #for getting OpenStreetMap data
library(patchwork) #for combining ggplots
library(ranger) #for random forest models
library(rdpower) #for RDiT power calculations 
library(rdrobust) #for RDiT analysis
library(rmweather) #for random forest normalisation
library(scico) #for nice colour palettes
library(sf) #for processing spatial data
library(tidyverse) #for data manipulation and plotting
library(weathermetrics) #for converting weather units
library(zoo) #additional stats/timeseries functions

## ---------------------------

source("1_spatial_preprocessing.R") #load and preprocess sensor metadata and create CAZ/spillover maps
source("2_timeseries_preprocessing.R") #load and preprocess sensor pollutant/traffic data
source('3_timeseries_normalisation.R') #rmweather normalisation of pollutant/traffic data
source('4_EDA.R') #exploratory data analysis of normalised data
source('5_RDiT.R') #sharp and donut RDiT analysis of normalised data
source('6_sensitivity_testing.R') #sensitivity testing of RDiT results for bandwidth and MDEs

