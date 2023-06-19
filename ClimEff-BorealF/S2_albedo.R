# ALBEDO ESTIMATION ON S2 RESOLUTION #


# load required packages #
library(raster)
library(GeoLight)
library(stringr)
library(dplyr)
library(terra)
library(sp)


###### 1. Albedo estimation following Li et al. 2018 #####

#' Title
#'
#' @param img_s2 Sentinel image name
#' @param img_brdf_modis  image name
#' @param img_albedo_modis image name
#' @param s2_band characte, band name
#' @param modis_brdf_band named vector of 3, in the order iso, vol, geo
#' @param modis_albedo_band named vector of 2, bsa, wsa
#' @param hvol 
#' @param hgeo 
#'
#' @return named list of 2, bsa and wsa
#' @export
#'



######Narrow-to-broadband converter#####
#narrow-to-broadband conversion with coefficients from Li et al. 2018
#(eq. 10)
ntb_converter <- function(band2,
                          band3,
                          band4,
                          c_0 = -0.0048,
                          c_B2 = 0.5673,
                          c_B3 = 0.1407,
                          c_B4 = 0.2359) {
  a = c_0 + (c_B2*band2 + c_B3*band3 + c_B4*band4)
  return(a)
}   


######Degree-to-Radian converter#####
DegRadconverter <- function(deg){
  return (deg * pi/180)
}

###########calculator for BRDF kernels############

#calculation of Kvol and Kgeo as BRDF kernels as described in Lucht et al. 2000
Kvol=function(teta_s, teta_v, phi){
  xi=acos(cos(teta_s)*cos(teta_v)+sin(teta_s)*sin(teta_v)*cos(phi))         
  Krt=((pi/2-xi)*cos(xi)+sin(xi))/(cos(teta_s)+cos(teta_v))-pi/4           
  return(Krt)
} 

sec=function(x){1/cos(x)} #secant

Kgeo=function(teta_s, teta_v, phi){
  b_r=1;h_b=1 
  teta_sPrim=atan(b_r*tan(teta_s))  
  teta_vPrim=atan(b_r*tan(teta_v))  
  xiPrim=acos(cos(teta_sPrim)*cos(teta_vPrim)+sin(teta_sPrim)*sin(teta_vPrim)*cos(phi))   
  D=sqrt(tan(teta_sPrim)**2+tan(teta_vPrim)**2-2*tan(teta_sPrim)*tan(teta_vPrim)*cos(phi))     
  t=acos(h_b*((sqrt(D**2+(tan(teta_sPrim)*tan(teta_vPrim)*sin(phi))**2))/(sec(teta_sPrim)+sec(teta_vPrim))))     
  O=(1/pi)*(t-sin(t)*cos(t))*(sec(teta_sPrim)+sec(teta_vPrim))     
  kLSR=O-sec(teta_sPrim)-sec(teta_vPrim)+.5*(1+cos(xiPrim))*sec(teta_sPrim)*sec(teta_vPrim)     
  return(kLSR)
}




##### 1.1: calculation of S2 Blue Sky Albedo#####

S2_WSA_BSA_calculator <- function(img_s2, 
                                  img_brdf_modis, 
                                  img_albedo_modis, 
                                  s2_band, 
                                  modis_brdf_band, 
                                  modis_albedo_band,
                                  hvol = hvol,
                                  hgeo = hgeo,
                                  beta = 0.2,
                                  kvol = kvol,
                                  kgeo = kgeo){
  
  #approximations for isotropic, volumetric and geometric kernels for white-sky integral accoring to look-up table in Lucht et al. 2000
  Hiso <- 1 
  Hvol <- 0.189184
  Hgeo <- -1.377622 
  
  ## 1.1 calculation of MODIS BSA and WSA##
  
  #(eq. 6)
  bsa <- img_brdf_modis[[as.character(modis_brdf_band["iso"])]] + img_brdf_modis[[as.character(modis_brdf_band["vol"])]] *hvol + img_brdf_modis[[as.character(modis_brdf_band["geo"])]]*hgeo
  
  #(eq. 7)
  wsa <- img_brdf_modis[[as.character(modis_brdf_band["iso"])]] + img_brdf_modis[[as.character(modis_brdf_band["vol"])]] *Hvol + img_brdf_modis[[as.character(modis_brdf_band["geo"])]]*Hgeo
  
  
  
  
  
  ## 1.2 calculation of MODIS bidirectional reflectance##
  
  #(eq. 4)
  R <- img_brdf_modis[[as.character(modis_brdf_band["iso"])]] + img_brdf_modis[[as.character(modis_brdf_band["vol"])]] * kvol + img_brdf_modis[[as.character(modis_brdf_band["geo"])]]*kgeo
  
  
  
  
  
  
  ## 1.3 calculation of MODIS BSA and WSA AN-ratios##
  
  #(eq. 3a)
  alpha_bsa <- bsa/R
  #(eq. 3b)
  alpha_wsa <- wsa/R
  
  
  
  ## 1.4 resampling of AN ratios to S2 resolution##
  
  alpha_bsa_res <- terra::resample(alpha_bsa, img_s2[[s2_band]], method="near")
  alpha_wsa_res <- terra::resample(alpha_wsa, img_s2[[s2_band]], method="near")
  
  
  
  ## 1.5 calculation of S2 BSA and WSA##
  
  #(eq. 8a)
  BSA_S2 <- alpha_bsa_res * img_s2[[s2_band]]
  #(eq. 8b)
  WSA_S2 <- alpha_wsa_res * img_s2[[s2_band]]
  
  
  
  ## 1.6 calculation of S2 Blue Sky Albedo##
  
  #(eq. 9)
  BluSA_S2 <- (1-beta) * BSA_S2 + beta * WSA_S2
  return(list(bsa = BSA_S2, wsa = WSA_S2, bluSa = BluSA_S2))
}





######## 1.2: application of calculations to all images #######

Albedo_all_calculator <- function(path_s2,
                                  path_modis_brdf,
                                  path_modis_albedo,
                                  datetime,
                                  S2_time,
                                  output_dir,
                                  teta_s,
                                  teta_v2,
                                  teta_v3,
                                  teta_v4,
                                  teta_v8,
                                  phi2,
                                  phi3,
                                  phi4
){
  
  
  ## 2.1 read images as raster ##
  s2 <- terra::rast(path_s2)
  modis_brdf  <- terra::rast(path_modis_brdf)
  modis_albedo  <- terra::rast(path_modis_albedo)
  
  
  
  ## 2.2 calculate kernel approximations for BSA and WSA calculation
  
  #approximations for volumetric and geometric kernels for the BSA according to look-up table in Lucht et al. 2000 
  #(eq. 5a)
  hvol <- -0.007574 - 0.070987*(DegRadconverter(teta_s))^2 + 0.307588*DegRadconverter(teta_s)^3
  #(eq. 5b)
  hgeo <- -1.284909 - 0.166314*(DegRadconverter(teta_s))^2+0.041840*DegRadconverter(teta_s)^3
  
  
  
  
  
  
  ## 2.3 applying BSA and WSA calculations to each band ##
  
  #Band 2 Sentinel
  bsa_wsa_list2 <- S2_WSA_BSA_calculator(img_s2 = s2, 
                                         img_brdf_modis = modis_brdf, 
                                         img_albedo_modis = modis_albedo, 
                                         s2_band = "B2", 
                                         modis_brdf_band = c(iso = "BRDF_Albedo_Parameters_Band3_iso", vol = "BRDF_Albedo_Parameters_Band3_vol", geo = "BRDF_Albedo_Parameters_Band3_geo"), 
                                         modis_albedo_band = c(bsa = "Albedo_BSA_Band3", wsa= "Albedo_WSA_Band3"),
                                         hvol = hvol,
                                         hgeo = hgeo,
                                         kvol = Kvol(teta_s=DegRadconverter(teta_s), teta_v=DegRadconverter(teta_v2), phi=DegRadconverter(phi2)),
                                         kgeo = Kgeo(teta_s=DegRadconverter(teta_s), teta_v=DegRadconverter(teta_v2), phi=DegRadconverter(phi2)))
  
  #Band 3 Sentinel
  bsa_wsa_list3 <- S2_WSA_BSA_calculator(img_s2 = s2, 
                                         img_brdf_modis = modis_brdf, 
                                         img_albedo_modis = modis_albedo, 
                                         s2_band = "B3", 
                                         modis_brdf_band = c(iso = "BRDF_Albedo_Parameters_Band4_iso", vol = "BRDF_Albedo_Parameters_Band4_vol", geo = "BRDF_Albedo_Parameters_Band4_geo"), 
                                         modis_albedo_band = c(bsa = "Albedo_BSA_Band4", wsa= "Albedo_WSA_Band4"),
                                         hvol = hvol,
                                         hgeo = hgeo,
                                         kvol = Kvol(teta_s=DegRadconverter(teta_s), teta_v=DegRadconverter(teta_v3), phi=DegRadconverter(phi3)),
                                         kgeo = Kgeo(teta_s=DegRadconverter(teta_s), teta_v=DegRadconverter(teta_v3), phi=DegRadconverter(phi3)))
  
  
  #Band 4 Sentinel
  bsa_wsa_list4 <- S2_WSA_BSA_calculator(img_s2 = s2, 
                                         img_brdf_modis = modis_brdf, 
                                         img_albedo_modis = modis_albedo, 
                                         s2_band = "B4", 
                                         modis_brdf_band = c(iso = "BRDF_Albedo_Parameters_Band1_iso", vol = "BRDF_Albedo_Parameters_Band1_vol", geo = "BRDF_Albedo_Parameters_Band1_geo"), 
                                         modis_albedo_band = c(bsa = "Albedo_BSA_Band1", wsa= "Albedo_WSA_Band1"),
                                         hvol = hvol,
                                         hgeo = hgeo,
                                         kvol = Kvol(teta_s=DegRadconverter(teta_s), teta_v=DegRadconverter(teta_v4), phi=DegRadconverter(phi4)),
                                         kgeo = Kgeo(teta_s=DegRadconverter(teta_s), teta_v=DegRadconverter(teta_v4), phi=DegRadconverter(phi4)))
  
  
  
  ## 2.4 Stacking and saving BSA, WSA and Blue Sky Albedo ##
  
  bsa_all <- c(bsa_wsa_list2$bsa, bsa_wsa_list3$bsa, bsa_wsa_list4$bsa)
  output_dir_spec <- file.path(output_dir, "BSA")
  if(!dir.exists(output_dir_spec)) dir.create(output_dir_spec)
  terra::writeRaster(bsa_all, file.path(output_dir_spec, paste0("BSA_S2_", datetime, ".tif")))
  
  
  wsa_all <- c(bsa_wsa_list2$wsa, bsa_wsa_list3$wsa, bsa_wsa_list4$wsa ) 
  names(wsa_all) <- c("B2","B3","B4")
  output_dir_spec <- file.path(output_dir, "WSA")
  if(!dir.exists(output_dir_spec)) dir.create(output_dir_spec)
  terra::writeRaster(wsa_all, file.path(output_dir_spec, paste0("WSA_S2_", datetime, ".tif")))
  
  
  blusa_all <- c(bsa_wsa_list2$bluSa, bsa_wsa_list3$bluSa, bsa_wsa_list4$bluSa) 
  names(blusa_all) <- c("B2","B3","B4")
  output_dir_spec <- file.path(output_dir, "BluSA")
  if(!dir.exists(output_dir_spec)) dir.create(output_dir_spec)
  terra::writeRaster(blusa_all, file.path(output_dir_spec, paste0("BluSA_S2_", datetime, ".tif")))
  
  
  
  
  
  
  ## 2.5 Narrow-to-broadband conversion Blue Sky Albedo ##
  
  blusa_ntb <- ntb_converter(band2 = blusa_all$B2,
                             band3 = blusa_all$B3,
                             band4 = blusa_all$B4)
  
  names(blusa_ntb) <- "blue_sky_albedo"
  blusa_ntb[blusa_ntb$blue_sky_albedo>1] <- NA
  blusa_ntb[blusa_ntb$blue_sky_albedo<0] <- NA
  
  
  
  
  ## 2.6 Stacking and saving Broadband Blue Sky Albedo ##
  
  output_dir_spec_ntb <- file.path(output_dir, "BluSA_ntb")
  if(!dir.exists(output_dir_spec_ntb)) dir.create(output_dir_spec_ntb)
  terra::writeRaster(blusa_ntb, file.path(output_dir_spec_ntb, paste0("BluSA_ntb_S2_", datetime, ".tif")))
  
  #return(list(bsa_all, wsa_all,blusa_all))
  
}




##### 1.3: reading images and creating dataframe #####

S2_files_dir= "input_data/GEE_download/ANNA_S2_ALL"
S2_files <- list.files(S2_files_dir, pattern = "*.tif$")
S2_files_full <- list.files(S2_files_dir, pattern = "*.tif$", full.names = TRUE)

modis_albe_dir <- "input_data/GEE_download/ANNA_MODIS_ALBE_ALL"
modis_albe_files <- list.files(modis_albe_dir, pattern = "*.tif$")
modis_albe_files_full <- list.files(modis_albe_dir, pattern = "*.tif$", full.names = TRUE)

modis_brdf_dir <- "input_data/GEE_download/ANNA_MODIS_BRDF_ALL"
modis_brdf_files <- list.files(modis_brdf_dir, pattern = "*.tif$")
modis_brdf_files_full <- list.files(modis_brdf_dir, pattern = "*.tif$", full.names = TRUE)


df_files <- data.frame(S2 = S2_files, S2_full = S2_files_full, Mod_Albe = modis_albe_files, Mod_Albe_full = modis_albe_files_full, Mod_brdf = modis_brdf_files, Mod_brdf_full = modis_brdf_files_full)


df_files$datetime <- sapply(strsplit(df_files$S2, "_"), "[[", 2)


df_files <- df_files %>%
  mutate(year = str_sub(datetime,1,4),
         month = str_sub(datetime,5,6),
         day = str_sub(datetime,7,8),
         hour = str_sub(datetime,9,10),
         min = str_sub(datetime,11,12),
         sec = str_sub(datetime,13,14)) %>%
  mutate(S2_time = paste0(year,"-",month,"-", day, " ", hour, ":", min,":", sec))

df_files$sun_zenith <- as.numeric(sapply(strsplit(df_files$S2, "_"), "[[", 3))/10000
df_files$view_zenith_B2 <- as.numeric(sapply(strsplit(df_files$S2, "_"), "[[", 4))/10000
df_files$view_zenith_B3 <- as.numeric(sapply(strsplit(df_files$S2, "_"), "[[", 5))/10000
df_files$view_zenith_B4 <- as.numeric(sapply(strsplit(df_files$S2, "_"), "[[", 6))/10000

df_files$sun_azimuth <- as.numeric(sapply(strsplit(df_files$S2, "_"), "[[", 8))/10000
df_files$view_azimuth_B2 <- as.numeric(sapply(strsplit(df_files$S2, "_"), "[[", 9))/10000
df_files$view_azimuth_B3 <- as.numeric(sapply(strsplit(df_files$S2, "_"), "[[", 10))/10000
df_files$view_azimuth_B4 <- as.numeric(sapply(strsplit(df_files$S2, "_"), "[[", 11))/10000

phi_B2 <- (df_files$view_azimuth_B2-df_files$sun_azimuth)%%180
df_files$phi_B2 <- abs(phi_B2)

phi_B3 <- (df_files$view_azimuth_B3-df_files$sun_azimuth)%%180
df_files$phi_B3 <- abs(phi_B3)

phi_B4 <- (df_files$view_azimuth_B4-df_files$sun_azimuth)%%180
df_files$phi_B4 <- abs(phi_B4)




##### 1.4: application of calculation to dataset #####

lapply(as.integer(rownames(df_files)), FUN=function(x){
  Albedo_all_calculator(path_s2 = df_files[x,"S2_full"],
                        path_modis_brdf =  df_files[x,"Mod_brdf_full"],
                        path_modis_albedo =  df_files[x,"Mod_Albe_full"],
                        datetime = df_files[x,"datetime"],
                        S2_time =df_files[x,"S2_time"],
                        output_dir = "output_data/Indiv_albedo",
                        teta_s = df_files[x,"sun_zenith"],
                        teta_v2 = df_files[x,"view_zenith_B2"],
                        teta_v3 = df_files[x,"view_zenith_B3"],
                        teta_v4 = df_files[x,"view_zenith_B4"],
                        phi2 = df_files[x,"phi_B2"],
                        phi3 = df_files[x,"phi_B3"],
                        phi4 = df_files[x,"phi_B4"]
  )
  
})





##### 2. Calculation of yearly mean albedo  #####

wd1 <- 'output_data/Indiv_albedo/BluSA_ntb'
blusa_path <- file.path(wd1)
blusa <- list.files(blusa_path, pattern = "*.tif$")
blusa_full <- list.files(blusa_path, pattern = "*.tif$", full.names = TRUE)
blusa_df <- data.frame(BluSA = blusa, BluSa_full = blusa_full)


# extracting date from images #
blusa_df$BluSa_dates <- str_sub(sapply(strsplit(blusa, "_"),  "[[", 4), 1, 8)
blusa_df <- blusa_df %>%
  mutate(year = str_sub(BluSa_dates,1,4),
         month = str_sub(BluSa_dates,5,6),
         day = str_sub(BluSa_dates,7,8)) %>%
  mutate(BluSa_time = paste0(year,"-",month,"-", day, " "))


# group by date #
blusa_df_2018 <- subset(blusa_df, year == 2018)
blusa_df_2019 <- subset(blusa_df, year == 2019)


# stacking different years #
stacked2018 <- terra::rast(blusa_df_2018$BluSa_full)
stacked2019 <- terra::rast(blusa_df_2019$BluSa_full)
stacked_allyears <- terra::rast(blusa_df$BluSa_full)


# yearly mean ##
mean_2018 <- terra::app(stacked2018, mean, na.rm = TRUE)
mean_2019 <- terra::app(stacked2019, mean, na.rm = TRUE)
mean_allyears <- terra::app(stacked_allyears, mean, na.rm = TRUE)


# # export as .tif files #
# output_dir2 <- "output_data/Indiv_albedo"
# output_dir_spec_mean2018 <- file.path(output_dir2, "BluSA_mean_fullyear_2018")
# if(!dir.exists(output_dir_spec_mean2018)) dir.create(output_dir_spec_mean2018)
# terra::writeRaster(mean_2018, file.path(output_dir_spec_mean2018, paste0("BluSA_mean_fullyear_2018",".tif")))
# output_dir_spec_mean2019<- file.path(output_dir2, "BluSA_mean_fulyear_2019")
# if(!dir.exists(output_dir_spec_mean2019)) dir.create(output_dir_spec_mean2019)
# terra::writeRaster(mean_2019, file.path(output_dir_spec_mean2019, paste0("BluSA_mean_fullyear_2019",".tif")))
# output_dir_spec_mean_allyears<- file.path(output_dir2, "BluSA_mean_fullyear_allyears")
# if(!dir.exists(output_dir_spec_mean_allyears)) dir.create(output_dir_spec_mean_allyears)
# terra::writeRaster(mean_allyears, file.path(output_dir_spec_mean_allyears, paste0("BluSA_mean_fullyear_allyears",".tif")))

# plot-level mean values were calculated in QGIS