# VALIDATION OF DATA #

# load packages
require(rgdal)
library(dplyr)
library(broom)
library(gdalUtilities)
library(readxl)
library(ggplot2)
library(tidyr)



#### 1. aggregation of S2 albedo to MODIS resolution ####

gdal_translate("data/output_data/Indiv_albe_V7/BluSA_mean_fullyear_allyears/BluSa_mean_fullyear_allyears.tif", "data/output_data/Indiv_albe_V7/BluSA_mean_fullyear_allyears/BluSa_mean_fullyear_allyears_aggregated.tif", r= "average", tr = c(500, 500), co = "COMPRESS=LZW")

# comparison of albedo values between S2 aggregated and MODIS product 
shape <- as.data.frame(readOGR(dsn = "data/output_data/raster_albedo_mean", layer = "albedo-differences-plots-2"))

df1$differences_S2_MOD <- shape$X_mean
df1_na <- na.omit(df1)
aggregate(df1_na$differences_S2_MOD, list(df1_na$SDS), FUN=mean)




###### 2. comparison with tower-based albedo of clear-cut #####

cc <- readxl::read_excel(path = "data/reference_sites_df.xlsx", sheet = "clearcuts")

cc$mean_albedo_both <- (cc$mean_albedo_tbs + cc$mean_albedo_tbn)/2

mean_cc <- mean(cc$mean_albedo_both) #0.307
median_cc <- median(cc$mean_albedo_both) #0.145





##### 3. comparison with tower-based albedo of Svartberget ####

s <- readxl::read_excel(path = "data/output_data/dataframes/reference_sites_df.xlsx", sheet = "svartberget")

mean_s <- mean(s$mean_albedo) #0.077
median_s <- median(s$mean_albedo) #0.078