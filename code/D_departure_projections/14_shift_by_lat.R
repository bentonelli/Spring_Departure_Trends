library(dplyr)

# Look at expected temperature changes

mdl_data <- readRDS("data/input/mdl_data/full_data_ss.rds")

dep_change1 <- readRDS("data/output/depature_projections/est_dep_CESM2_SSP2_4p5.rds")
dep_change2 <- readRDS("data/output/depature_projections/est_dep_CNRM-ESM2-1_SSP2_4p5.rds")
dep_change3 <- readRDS("data/output/depature_projections/est_dep_MIROC6_SSP2_4p5.rds")

sp_it <- mdl_data %>% 
  select(c("sp_site","cell","elev",
           "lat","lon",
           "species","mean_dep",
           "mean_temp60","sd_temp_60",
           "mean_total_precip60","sd_total_precip60")) %>%
  unique()

sp_it$dep_change_mean <- rowMeans(((dep_change1+dep_change2+dep_change3)/3)[77:86])
plot(sp_it$lat,sp_it$dep_change_mean,cex=.1,col="dodgerblue4",xlab="Latitude",ylab="Departure Change 2091-2100")
abline(h=0)
abline(v=25)

sp_it %>% 
  filter(abs(lat) > 23.4) %>%
  summarise(count_over_0 = sum(dep_change_mean>0),
            count_below_0 = sum(dep_change_mean<0),
            median_change = median(dep_change_mean),
            total= n())

sp_it %>% 
  filter(abs(lat) <= 23.4) %>%
  summarise(count_over_0 = sum(dep_change_mean>0),
            count_below_0 = sum(dep_change_mean<0),
            median_change = median(dep_change_mean),
            total= n())

sp_level_mean  <- sp_it %>% group_by(species) %>% summarise(mean_change = mean(dep_change_mean))
sum(sp_level_mean$mean_change > 0)
sum(sp_level_mean$mean_change < 0)
