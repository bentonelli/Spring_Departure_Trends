library(dplyr)
library(MCMCvis)
library(ggplot2)
library(birdcolors)
library(ggbeeswarm)
library(tidyr)
library(terra)

#Model data
mdl_data <- readRDS("data/input/mdl_data/full_data_ss.rds")
indexing_tab <- readRDS("data/input/mdl_data/indexing_tab_data_ss.rds")

# Look at expected temperature changes
temp_abs_change1 <- readRDS("data/output/depature_projections/temp_abs_change_CESM2_SSP2_4p5.rds")
temp_abs_change2 <- readRDS("data/output/depature_projections/temp_abs_change_CNRM-ESM2-1_SSP2_4p5.rds")
temp_abs_change3 <- readRDS("data/output/depature_projections/temp_abs_change_MIROC6_SSP2_4p5.rds")

sp_it <- mdl_data %>% 
  select(c("sp_site","cell","elev",
           "lat","lon",
           "species","mean_dep",
           "mean_temp60","sd_temp_60",
           "mean_total_precip60","sd_total_precip60")) %>%
  unique()

sp_it$ssp4p5_abs_change_temp_eoc <- rowMeans(((temp_abs_change1+temp_abs_change2+temp_abs_change3)/3)[77:86])
sp_it$ssp4p5_anom_change_temp_eoc <- sp_it$ssp4p5_abs_change_temp_eoc/sp_it$sd_temp_60 

bc_pal <- bird_colors("Variable Sunbird")

area_avg <- sp_it %>%
  group_by(cell,elev,lon,lat) %>%
  summarise(abs_temp_change = mean(ssp4p5_abs_change_temp_eoc),
            anom_temp_change = mean(ssp4p5_anom_change_temp_eoc))

pp1 <- ggplot() +
  geom_point(data=sp_it,aes(x=mean_temp60-273.15,
                            y=ssp4p5_abs_change_temp_eoc,
                            col=abs(lat)),size=0.1,alpha=.5) +
  geom_hline(yintercept = 0,col="black",lty=2) +
  scale_color_bird(bc_pal) +
  ylab("Temperature Change - Absolute (C)") +
  xlab("Mean Historical Temperature (C)") +
  theme_classic(base_size = 16)

pp2 <- ggplot() +
  geom_point(data=sp_it,aes(x=mean_temp60-273.15,
                            y=ssp4p5_anom_change_temp_eoc,
                            col=abs(lat)),size=0.1,alpha=.5) +
  geom_hline(yintercept = 0,col="black",lty=2) +
  scale_color_bird(bc_pal) +
  ylab("Temperature Change - Anomaly (SD)") +
  xlab("Mean Historical Temperature") +
  theme_classic(base_size = 16)


# Look at expected precipitation changes
precip_abs_change1 <- readRDS("data/output/depature_projections/precip_abs_change_CESM2_SSP2_4p5.rds")
precip_abs_change2 <- readRDS("data/output/depature_projections/precip_abs_change_CNRM-ESM2-1_SSP2_4p5.rds")
precip_abs_change3 <- readRDS("data/output/depature_projections/precip_abs_change_MIROC6_SSP2_4p5.rds")

sp_it$ssp4p5_abs_change_precip_eoc <- rowMeans(((precip_abs_change1+precip_abs_change2+precip_abs_change3)/3)[77:86])
sp_it$ssp4p5_anom_change_precip_eoc <- sp_it$ssp4p5_abs_change_precip_eoc/sp_it$sd_total_precip60 

pp3 <- ggplot() +
  geom_point(data=sp_it,aes(x=mean_total_precip60*1000,
                            y=1000*ssp4p5_abs_change_precip_eoc,
                            col=abs(lat)),size=0.25,alpha=.75) +
  geom_hline(yintercept = 0,col="black",lty=2) +
  scale_color_bird(bc_pal) +
  scale_x_log10()+
  ylab("Precipitation Change - Absloute (mm/day)") +
  xlab("Mean Historical Precipitation (mm/day)") +
  theme_classic(base_size = 16)

pp4 <- ggplot() +
  geom_point(data=sp_it,aes(x=mean_total_precip60*1000,
                            y=ssp4p5_anom_change_precip_eoc,
                            col=abs(lat)),size=0.25,alpha=.75) +
  geom_hline(yintercept = 0,col="black",lty=2) +
  scale_color_bird(bc_pal) +
  scale_x_log10()+
  ylab("Precipitation Change - Anomaly (SD)") +
  xlab("Mean Historical Precipitation (mm/day)") +
  theme_classic(base_size = 16)

pp1
pp2
pp3
pp4

