library(dplyr)
library(ggplot2)
library(rstan)
library(dggridR)
library(sf)
library(terra)
library(MCMCvis)
library(rsoi)
library(shinystan)
library(mapdata)

# Read in phenology, enviro, and spatial data
sp_dep <- readRDS("data/input/phenology/spr_dep_core.rds")
era_climate_data <- readRDS("data/input/environmental/enviro_data_ts_60.rds")
grid_elev_bands <- st_read("data/input/spatial/grid_elev_bands.shp")

# Get unique sites in phenology dataset
unq_cl_elev <- sp_dep %>% dplyr::select(cell,elev) %>% unique()

# Merge with elev, cell
grid_elev_bands <- merge(grid_elev_bands,unq_cl_elev,by=c("elev","cell"))

# Get centroids for lat-lon data
elev_cell_centroids <- st_centroid(grid_elev_bands)
coordinates <- st_coordinates(elev_cell_centroids)
cell_elev_points <- data.frame(cell = elev_cell_centroids$cell,
                               elev=elev_cell_centroids$elev,
                               lon = coordinates[,1],
                               lat = coordinates[,2])

# Species, site combos with at least 3 data points, site level SD < 20 days
sp_dep <- sp_dep %>% 
  group_by(species,cell,elev) %>% 
  mutate(sp_s_n = n(),
         sd_s = sd(trailing_hm_mean)) %>% 
  filter(sp_s_n >= 3 & sd_s < 20) %>%
  ungroup()

# Remove species with less than 20 data points
sp_dep <- sp_dep %>% 
  group_by(species) %>% 
  mutate(count = n()) %>%
  filter(count >= 20) %>%
  ungroup()

nrow(sp_dep)


unq_cl_elev <- sp_dep %>% 
  dplyr::select(cell,elev) %>% 
  unique()

grid_elev_bands <- merge(grid_elev_bands,unq_cl_elev,by=c("elev","cell"))

sp_dep_temp <- sp_dep

#Get average departure dates at the species-site level
mean_dep_dates_sp_site <- sp_dep_temp %>% 
  group_by(species,cell,elev) %>%
  summarise(mean_dep = round(mean(trailing_hm_mean)))

sp_dep_temp <- merge(sp_dep_temp,mean_dep_dates_sp_site,by=c("species","cell","elev")) 

#Now add on climate stuff 
era_climate_long <- bind_rows(era_climate_data)

era_climate_long_cells <- era_climate_long %>% 
  group_by(cell,elev,day) %>%
  summarise(
    #Temp
    mean_temp60 = mean(rm_temp60,na.rm=TRUE),
    sd_temp_60 = sd(rm_temp60,na.rm=TRUE),
    #Precip
    mean_total_precip60 = mean(rm_precip60,na.rm=TRUE),
    sd_total_precip60 = sd(rm_precip60,na.rm=TRUE))

colnames(era_climate_long_cells)[3] <- "mean_dep"

sp_temp_climate <- merge(sp_dep_temp,era_climate_long_cells,by=c("cell","elev","mean_dep"))

sp_temp_climate <- merge(sp_temp_climate,era_climate_long,
                         by.x=c("year","cell","elev","mean_dep"),
                         by.y=c("year","cell","elev","day"))

#Clear memory
sp_temp_climate$ID <- sp_temp_climate$temp <- sp_temp_climate$precip <- sp_temp_climate$likely_migratory <- era_climate_long <- NULL

sp_temp_climate$temp60_stand <- (sp_temp_climate$rm_temp60 - sp_temp_climate$mean_temp60)
sp_temp_climate$precip60_stand  <- (sp_temp_climate$rm_precip60 - sp_temp_climate$mean_total_precip60)

sp_temp_climate$temp60_anom <- sp_temp_climate$temp60_stand/sp_temp_climate$sd_temp_60
sp_temp_climate$precip60_anom <- sp_temp_climate$precip60_stand/sp_temp_climate$sd_total_precip60

test_dataset <- sp_temp_climate %>% 
  filter(!is.na(temp60_anom)) %>%
  arrange(species,cell,elev,year)

#Make sure all of the species still have >= 20 data points
test_dataset <- test_dataset %>% 
  group_by(species) %>% 
  mutate(sp_n = n())

#And >= 2 sites
unq_site_check <- test_dataset %>% 
  select(species,cell,elev) %>%
  unique() %>%
  group_by(species) %>%
  summarise(unq_sp_site_count = n()) %>%
  filter(unq_sp_site_count >= 2)

test_dataset <- test_dataset %>% 
  filter(sp_n >= 20 & species %in% unq_site_check$species) %>%
  ungroup()

sp_dep_count <- test_dataset %>% 
  group_by(species) %>% 
  summarise(count = n())

sp_ind <- test_dataset %>% select(species) %>% unique()
sp_ind$sp <- 1:nrow(sp_ind)

sp_site_ind <- test_dataset %>% 
  select(species,cell,elev) %>% 
  unique()
sp_site_ind$sp_site <- 1:nrow(sp_site_ind)

site_ind <- test_dataset %>% 
  select(cell,elev) %>% 
  unique()
site_ind$site_ind <- 1:nrow(site_ind)

test_dataset <- merge(test_dataset,sp_ind,by=c("species"))
test_dataset <- merge(test_dataset,sp_site_ind,by=c("species","cell","elev"))
test_dataset <- merge(test_dataset,site_ind,by=c("cell","elev"))

test_dataset <- test_dataset %>% 
  arrange(species,cell,elev,year)

indexing_tab <- test_dataset %>% 
  select(sp,sp_site,site_ind,cell,elev,mean_dep,mean_temp60,mean_total_precip60,lat,lon) %>% 
  unique()

indexing_tab <- indexing_tab %>% 
  group_by(sp) %>% 
  mutate(mean_temp_niche = mean(mean_temp60),
         temp_niche = (mean_temp60 - mean(mean_temp60))/sd(mean_temp60),
         mean_precip_niche = mean(mean_total_precip60),
         precip_niche = (mean_total_precip60 - mean(mean_total_precip60))/sd(mean_total_precip60),
         mean_lat = mean(lat),
         rel_lat = lat - mean(lat),
         mean_lon = mean(lon),
         rel_lon = lon - mean(lon)) %>%
  ungroup()

sp_avg_climate_metrics <- indexing_tab %>%
  select(sp,mean_temp_niche,mean_precip_niche,mean_lat) %>%
  unique()

print(unique(test_dataset$species))
print(nrow(test_dataset))

#Clear memory
era_climate_long_cells <- era_climate_data <- NULL


### Model Data ###
data_mdl <- list(
  N = nrow(test_dataset),
  yy = test_dataset$trailing_hm_mean - mean(test_dataset$trailing_hm_mean),
  yy_unc = test_dataset$trailing_hm_sd,
  temp = test_dataset$temp60_anom,
  precip = test_dataset$precip60_anom,
  
  Nsp = length(unique(test_dataset$sp)),
  sp = test_dataset$sp,
  
  Nsp_c = length(unique(test_dataset$sp_site)),
  sp_c = test_dataset$sp_site,
  Nsite = length(unique(test_dataset$site_ind)),
  
  sp_match = indexing_tab$sp,
  sp_cl_match = indexing_tab$sp_site,
  site_match = indexing_tab$site_ind,
  
  lat = indexing_tab$rel_lat,
  
  TN = indexing_tab$temp_niche,
  PN = indexing_tab$precip_niche,
  
  MT = (sp_avg_climate_metrics$mean_temp_niche - mean(sp_avg_climate_metrics$mean_temp_niche)),
  MP = (sp_avg_climate_metrics$mean_precip_niche - mean(sp_avg_climate_metrics$mean_precip_niche))*1000,
  
  ML = (sp_avg_climate_metrics$mean_lat - mean(sp_avg_climate_metrics$mean_lat))
)

saveRDS(data_mdl,"data/input/mdl_data/list_format_data_ss_test.rds")
saveRDS(indexing_tab,"data/input/mdl_data/indexing_tab_data_ss_test.rds")
saveRDS(test_dataset,"data/input/mdl_data/full_data_ss_test.rds")

