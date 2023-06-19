###### radiative forcing models ######

# load packages #
library(ggplot2)
library(readxl)
library(dplyr)
library(rgdal)
library(writexl)


# read data # 
df <- read_excel("data/plot_data/df.xlsx", sheet = "data") #plot-level stand data 
radiation <- read_excel("data/radiation_data/radiation_krycklan.xlsx", sheet = "data") #mean monthly and annual sw downward solar radiation over krycklan from NASA portal




# definitions # 
# A: radiative efficiency of CO2 (0.0198*10^-13)*1000 Wm^-2g^-1)
# f: CO2 fractions based on Frolking and Roulet (2007), table 1
# t1: time horizon of emission impact (t in Frolking and Roulet (2007) eq.1)
# t2: years of emission (t' in Frolking and Roulet (2007) eq. 1)
# tau: lifetime constants of fractions based on Frolking and Roulet (2007), table 1




###### 1. CO2-induced RF based on Frolking and Roulet (2007) eq. 1 #######
RF_calculator <- function(f,
                          tau,   
                          flux,  
                          t1,
                          t2,
                          A=(0.0198*10^-13)*1000 #radiative efficiency of CO2 in Wm^-2g^-1)
){
  return(A*f*flux*(44/12)*exp((t2-t1)/tau))
}


rf_pl_ym <- function(NEP
                     ){
  return(RF_calculator(0.176,10^8, NEP, t1=1, t2=2)+
           RF_calculator(0.138,421,NEP, t1=1, t2=2)+
           RF_calculator(0.186,70.6,NEP, t1=1, t2=2)+
           RF_calculator(0.242,21.4,NEP, t1=1, t2=2)+
           RF_calculator(0.259,3.42, NEP, t1=1, t2=2)
  )
}


df$NEP_RF <- mapply(rf_pl_ym, (df$NEP * (-1)))








###### 2. albedo-induced RF, based on Lohila et al. (2010) #####


### 2.1 convert radiation data into correct unit 

daylight_hrs <- 12.475 #mean annual daylight hours for 2019
radiation$ANN <- as.numeric(radiation$ANN) 
radiation$ANN <- radiation$ANN*1000 #convert from kW to W
radiation$sw_dwn_day <- radiation$ANN/daylight_hrs # mean daily sw downward radiation: 280.5611 W/m^2/day
radiation$sw_dwn_yr <- radiation$sw_dwn_day*365 #mean annual sw downward radiation: 102404.8 W/m^2/yr --> with this value calculate absorbed radiation and RF per plot


### 2.2 calculation of absorbed radiation per plot 

cc_albedo <- 0.239 #mean albedo of clear-cut reference site 
cc_abs <- radiation$sw_dwn_yr - (radiation$sw_dwn_yr*cc_albedo) #mean yearly absorbed solar radiation of clear-cut reference plot (W/m^2/yr)
df$sw_abs <- radiation$sw_dwn_yr - (radiation$sw_dwn_yr*df$albedo)



#### 2.3 calculate local radiative forcing 
df$albedo_RF_loc <- df$sw_abs - cc_abs #local RF caused by albedo change between cc and stand



#### 2.4 calculate global radiative forcing 
df$albedo_RF_glob <- df$albedo_RF_loc * (10/(5.1*10^14)) #global RF, scaled by earth's surface (see Lohila et al. 2010)


##### 2.5 calculate combined radiative forcing 
df$RF_combined <- df$albedo_RF_glob + df$NEP_RF



