# STATISTICS AND PLOTS #

# load packages # 
library(rgdal)
library(broom)
library(maptools)
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggpubr)
library(corrr)
library(dplyr)
library(patchwork)
library(Metrics)



# read data #

pm <- as.data.frame(readOGR(
  dsn = "data/output_data/raster_albedo_mean/50_plots_VIS_albedo_mean_allyears",
layer = "50_plots_VIS_albedo_mean_allyears"
))

df1 <- read_xlsx("data/plot_data/df1.xlsx", sheet = "Sheet1")
df1$SDS <- as.factor(df1$SDS)


 #### 1. test data for normality and homogenity of variances #####
shapiro.test(df1$albedo_mean) #albedo is not normally distributed so always use spearman method


 #### 2. correlation tests ####

plot(df1$age, df1$albedo_mean)
cor.test(df1$age, df1$albedo_mean, method = "spearman") #tied data so use kendall method instead
cor.test(df1$age, df1$albedo_mean, method = "kendal") # p < 0.001, tau = -0.512
cor.test(df1$SDS, df1$albedo_mean, method = "kendal") # p<0.001, tau = -0.505
cor.test(df1$LAImax, df1$albedo_mean, method = "spearman") #p<0.001, rho = -0.742 -> relationship as expected (higher LAI for lower albedo)
cor.test(df1$NEP, df1$albedo_mean, method = "spearman") # p < 0.05, rho = -0.355
cor.test(df1$age, df1$NEP, method = "kendal")

cor.test(df1$tree_density, df1$albedo_mean, method = "kendal") # p < 0.04, tau = -0.208
cor.test(df1$SPM, df1$albedo_mean, method = "kendal") #not significant




#### 4. plots for significant correlations#####

## age and albedo

#as scatterplot
scatter <- ggplot(df1, aes(x=age, y=albedo))+
  geom_point(aes(color = SDS))+
  geom_smooth(method = "loess")+
  scale_color_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="stand age (yr)", y="albedo")+
  labs(color = "stand age class") + 
  stat_cor(method = "kendall", label.x = 100, label.y = 0.12, cor.coef.name = "tau")+
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.7, 0.84))
  
#as boxplot
box <- ggplot(df1, aes(x=age, y=albedo))+
  geom_boxplot(aes(fill = SDS))+
  scale_fill_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="stand age (yr)", y="albedo")+
  labs(fill = "stand age class") + 
  stat_cor(method = "kendall", label.x = 100, label.y = 0.12, cor.coef.name = "tau")+
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.7, 0.84))

scatter + box


## LAI and albedo
ggplot(df1, aes(x=LAImax, y=albedo))+
  geom_point(aes(color = SDS))+
  geom_smooth(method = "loess")+
  scale_color_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="LAI (m^2/m^2)", y="albedo")+
  labs(color = "stand age class") + 
  stat_cor(method = "spearman", label.x = 2, cor.coef.name = "rho")+
  theme_classic(base_size = 16)+
  theme(legend.background = element_rect(fill = "white", color = "black"), legend.position = c(0.75, 0.7), legend.title = element_text(size=12),legend.text = element_text(size=8))+
  ylim(0, 0.25)
  



## tree density and albedo
df_na <- na.omit(df1)
ggplot(df_na, aes(x=tree_density, y=albedo))+
  geom_point(aes(color = SDS))+
  geom_smooth(method = "loess", na.rm = TRUE)+
  scale_color_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="tree density", y="albedo",
       title = "albedo and tree density")+
  labs(subtitle= "kendall correlation between albedo and tree density")+
  labs(color = "stand age class") + 
  stat_cor(method = "kendall", label.x = 1500, cor.coef.name = "tau", na.rm = TRUE)




#### 3. pca ####

df_na$SDS <- as.numeric(df_na$SDS)
df_na$SPM <- as.numeric(df_na$SPM)
pca1 <- prcomp(df_na[, 5:11], center = TRUE, scale = TRUE) #pca excluding the dependent variable (plot id) and min, max stdv, variance of albedo

print(pca1)
pca1$rotation <- -1*pca1$rotation
pca1$rotation
pca1$x <- -1*pca1$x

biplot(pca1, scale = 0)
summary(pca1)
var_explained = pca1$sdev^2 / sum(pca1$sdev^2)

qplot(c(1:7), var_explained) + 
  geom_line() + 
   xlab("Principal Component") + 
   ylab("Variance Explained") +
   ggtitle("Scree Plot") +
   ylim(0, 1)
head(df_na[, 5:11])



####### 4. plotting NEP #####

#as scatterplot
ggplot(df1, aes(x=age, y=NEP))+
  geom_point(aes(color = SDS))+
  geom_smooth(method = "loess")+
  scale_color_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="stand age (yr)", y="NEP [gC/m^2/yr]")+
  labs(color = "stand age class") + 
  stat_cor(method = "kendall", label.x = 110, label.y = 300, cor.coef.name = "tau")+
  geom_hline(yintercept = 0)+
  theme_classic(base_size = 14)+
  theme(legend.position = c(0.74, 0.175),legend.background = element_rect(fill = "white", color = "black"), legend.title = element_text(size=12),legend.text = element_text(size=8))
  
  

#as boxplot
ggplot(df1, aes(x=age, y=NEP))+
  geom_boxplot(aes(fill = SDS))+
  scale_fill_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="stand age class", y="NEP [gC/m^2/yr]")+
  labs(fill = "stand age class")+
  geom_hline(yintercept = 0)+
  theme_classic(base_size = 16)+
  theme(legend.title = element_text(size=12),legend.text = element_text(size=8))


#mean for stand age classes
NEPmeans <- df1 %>%
  group_by(SDS) %>%
  summarise_at(vars(NEP), list(name = mean))
NEPmeans
aggregate(df1$NEP, list(df$SDS), FUN=mean) 

NEPmedian <- df1 %>%
  group_by(SDS) %>%
  summarise_at(vars(NEP), list(name = median))
NEPmedian


RFmedian <- df1 %>%
  group_by(SDS) %>%
  summarise_at(vars(albedo), list(name = median))
RFmedian


####### 5. plotting RF ######

options(scipen = 999)

## c-flux induced RF 
ggplot(df1, aes(x=age, y=(NEP_RF*1000000000)))+
  geom_point(aes(color = SDS))+
  geom_smooth(method = "loess", na.rm = TRUE)+
  scale_color_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="stand age [yr]", y="RF [nW/m^2/yr]")+
  labs(color = "stand age class")+
  geom_hline(yintercept = 0)+
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.74, 0.84),legend.title = element_text(size=12),legend.text = element_text(size=8))+
  ylim(-3,1.5)
  

ggplot(df1, aes(x=age, y=(NEP_RF*1000000000)))+
  geom_boxplot(aes(fill = SDS))+
  scale_fill_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="stand age [yr]", y="RF [nW/m^2/yr]",
       title = "C-flux induced radiative forcing for stand age classes")+
  labs(subtitle= "based on NEP of the respective plot")+
  labs(fill = "stand age class")+
  geom_hline(yintercept = 0)+
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.74, 0.84),legend.title = element_text(size=12),legend.text = element_text(size=8))+
  ylim(-3,1.5)



## albedo-induced RF local
ggplot(df1, aes(x=age, y= (albedo_RF_loc/1000)))+
  geom_point(aes(color = SDS))+
  geom_smooth(method = "loess", na.rm = TRUE)+
  scale_color_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="stang age [yr]", y="RF [kW/m^2/yr]",
       title = "albedo induced local radiative forcing")+
  labs(subtitle= "based on relative albedo between initiation stage and the respective plot")+
  labs(color = "stand age class")+
  geom_hline(yintercept = 0)



## albedo-induced RF global
ggplot(df1, aes(x=age, y= (albedo_RF_glob*1000000000)))+
  geom_point(aes(color = SDS))+
  geom_smooth(method = "loess", na.rm = TRUE)+
  scale_color_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="stang age [yr]", y="RF [nW/m^2/yr]")+
  labs(color = "stand age class")+
  geom_hline(yintercept = 0)+
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.8, 0.35),legend.title = element_text(size=12),legend.text = element_text(size=8))


ggplot(df1, aes(x=age, y=(albedo_RF_glob*1000000000)))+
  geom_boxplot(aes(fill = SDS))+
  scale_fill_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="stand age [yr]", y="RF [nW/m^2/yr]",
       title = "albedo induced global radiative forcing for stand age classes")+
  labs(fill = "stand age class")+
  geom_hline(yintercept = 0)




## comparing both RF parts

df1_long <- gather(df1, key = "component", value = "RF", 12,14)

plot1 <- ggplot(df1_long, aes(x=age, y=(RF*1000000000), color=component)) +
  geom_point()+
  geom_smooth(method = "loess", na.rm = TRUE)+
  scale_color_hue(labels = c("albedo", "C-flux"))+
  labs(x="stand age [yr]", y="RF [nW/m^2/yr]")+
  geom_hline(yintercept = 0)+
  theme_classic(base_size = 16)+
  theme(legend.background = element_rect(fill = "white", color = "black"), legend.position = c(0.75, 0.9),legend.title = element_text(size=12),legend.text = element_text(size=10))+
  ylim(-3,1.5)


## combined RF
plot2 <- ggplot(df1, aes(x=age, y= (RF_combined*1000000000)))+
  geom_point(aes(color = SDS))+
  geom_smooth(method = "loess", na.rm = TRUE)+
  scale_color_hue(labels = c("initiation", "young", "middle-aged", "mature", "old-growth"))+
  labs(x="stang age [yr]", y="RF [nW/m^2/yr]")+
  labs(color = "stand age class")+
  geom_hline(yintercept = 0)+
  theme_classic(base_size = 16)+
  theme(legend.background = element_rect(fill = "white", color = "black"), legend.position = c(0.7, 0.84),legend.title = element_text(size=12),legend.text = element_text(size=10))+
  ylim(-3,1.5)

plot1 + plot2





###### 6. fitting models to albedo data ####


# 6.1 age-albedo relationship

plot(df1$age, df1$albedo) # plot age vs. albedo
df1$LNage <-log(df1$age) # log-transform the age and store in a new variable
plot(df1$LNage, df1$albedo) # plot the transformed relationship - looks more linear

LM <- lm(albedo~LNage,data=df1) # create linear model
summary(LM) # show statistics. The "Estimate" column contains the estimated model parameters.

plot(df1$albedo, fitted(LM)) # plot the model's predicted values vs. the reference values
abline(0,1) # add 1:1 line


plot(df1$age, df1$albedo,xlab="Age",ylab="Albedo")
df1$fitted <- fitted(LM)
df1 <- df1[order(df1$age),]
lines(df1$age,df1$fitted)

summary(df1$fitted)
rmse <- rmse(df1$albedo, df1$fitted)
rmse/(0.207745 - (-0.007988))

plot(df1$age,fitted(LM),xlab="Age",ylab="Albedo")


ggplot() + 
  geom_point(df1, mappin = aes(x = age, y = albedo)) + 
  geom_line(df1, mapping = aes(x=age, y = fitted))+
  labs(x="age (yr)", y="albedo")+
  theme_classic(base_size = 16)+
  ylim(0, 0.25)





# 6.2 LAI-albedo relationship

plot(df1$LAImax, df1$albedo) # plot age vs. albedo
df1$LN_LAI <-log(df1$LAImax) # log-transform the age and store in a new variable
plot(df1$LN_LAI, df1$albedo) # plot the transformed relationship - looks more linear

LM2 <- lm(albedo~LN_LAI,data=df1) # create linear model
summary(LM2) # show statistics. The "Estimate" column contains the estimated model parameters.

plot(df1$albedo, fitted(LM2)) # plot the model's predicted values vs. the reference values
abline(0,1) # add 1:1 line


plot(df1$LAImax, df1$albedo,xlab="LAI (m^2/m^2)",ylab="Albedo")
df1$fitted2 <- fitted(LM2)
df1 <- df1[order(df1$LAImax),]
lines(df1$LAImax,df1$fitted2)

summary(df1$fitted2)
rmse_2 <- rmse(df1$albedo, df1$fitted2)
rmse/(0.214163 - (-0.003914))
rmse_2

ggplot() + 
  geom_point(df1, mappin = aes(x =LAImax, y = albedo)) + 
  geom_line(df1, mapping = aes(x=LAImax, y = fitted2))+
  labs(x="LAI (m^2/m^2)", y="albedo")+
  theme_classic(base_size = 16)+
  ylim(0, 0.25)






#### for NEP-age relationship######
#pre-fitted dataset from Peichl et al. (2023)

fitted <- read_xlsx("data/NEP_data/NEP_table/plot_level_own/NEP-fitted.xlsx")

ggplot() + 
  geom_point(df1, mapping = aes(x = age, y = NEP)) + 
  geom_line(fitted, mapping = aes(x=Age, y = NEP, color="red"))+
  labs(x="stand age (yrs)", y="NEP (gC/m^2/yr)")+
  theme_classic(base_size = 16)+
  geom_hline(yintercept = 0)+
  ylim(-200, 400)



###### for age- albedoRF - relationship
plot(df1$age, df1$albedo_RF_glob) # plot age vs. albedo
df1$LNage <-log(df1$age) # log-transform the age and store in a new variable
plot(df1$LNage, df1$albedo_RF_glob) # plot the transformed relationship - looks more linear

LM4 <- lm(albedo_RF_glob~LNage,data=df1) # create linear model
summary(LM4) # show statistics. The "Estimate" column contains the estimated model parameters.

plot(df1$albedo_RF_glob, fitted(LM4)) # plot the model's predicted values vs. the reference values
abline(0,1) # add 1:1 line


plot(df1$age, df1$albedo_RF_glob,xlab="Age",ylab="Albedo")
df1$fitted4 <- fitted(LM4)
df1 <- df1[order(df1$age),]
lines(df1$age,df1$fitted4)

summary(df1$fitted4)
rmse <- rmse(df1$albedo_RF_glob, df1$fitted4)
rmse


ggplot() + 
  geom_point(df1, mappin = aes(x = age, y = (albedo_RF_glob)*1000000000)) + 
  geom_line(df1, mapping = aes(x=age, y = (fitted4)*1000000000))+
  labs(x="age (yr)", y="RF [nW/m^2/yr]")+
  theme_classic(base_size = 16)