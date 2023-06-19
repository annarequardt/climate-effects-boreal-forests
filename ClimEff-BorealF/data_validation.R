# VALIDATION OF DATA #

# load packages
require(rgdal)
library(dplyr)
library(broom)
library(gdalUtilities)
library(readxl)
library(ggplot2)
library(tidyr)







###### 1. comparison with tower-based albedo of clear-cut #####

cc <- readxl::read_excel(path = "input_data/reference-sites.xlsx", sheet = "clearcuts")

cc$mean_albedo_both <- (cc$mean_albedo_tbs + cc$mean_albedo_tbn)/2

mean_cc <- mean(cc$mean_albedo_both) #0.307
median_cc <- median(cc$mean_albedo_both) #0.145





##### 2. comparison with tower-based albedo of Svartberget ####

s <- readxl::read_excel(path = "input_data/reference-sites.xlsx", sheet = "svartberget")

mean_s <- mean(s$mean_albedo) #0.077
median_s <- median(s$mean_albedo) #0.078