# Climate effects of managed boreal forests
 This repository contains codes to estimate the radiative forcing of forests based on their albedo and carbon sequestration. The main code contains the calculation of high-resolution albedo values based on remotely sensed MODIS and S2 data, following the approach by Li et al. (2018).
 The repository is part of the M.Sc. Thesis project by Anna Valeria Requardt (SLU, Uppsala, Sweden, 2023).


## What this repository contains
This repository contains a Google Earth Engine code and an R project with which albedo and radiative forcing can be calculated. The Google Earth Engine code is for the data download (S2 reflection and MODIS BRDF data), whilst the R code is for data processing. 
The R project contains three codes for data processing: 
1) the albedo calculation based on S2 reflectance and MODIS BRDF data (following the approach by Li et al. 2018)
2) the calculation of Radiative Forcing caused by relative albedo differences and (pre-existing) Carbon-flux data (NEP) of the same forest area
3) the comparison of the computed albedo values with tower-based albedo measurements of reference sites
Moreover, the R project contains the used input data and produced output data. 
   

## How to use this repository
Users interested in high-resolution albedo calculation are recommended to firstly download data over the ROI using the GEE code and then applying the R code "S2_albedo.R" to the downloaded dataset. Users that are further interested in calculating the Radiative Forcing based on albedo changes / relative albedo differences are recommended to further process the albedo data using the code "RF_models.R". 

For questions and feedback contact the owner. 
