# ERA5 processing for cells

library(dplyr)
library(terra)
library(sf)
library(lubridate)
library(ncdf4)
library(zoo)


grid_elev_bands <- st_read("data/input/spatial/grid_elev_bands.shp")
names(grid_elev_bands)[1] <- "elev"

sp_dep <- readRDS("data/input/phenology/spring_dep.rds")

sp_counts <- sp_dep %>% group_by(species) %>% summarise(count = n())
sp_include <- sp_counts$species[which(sp_counts$count >= 20)]
sp_dep_filt <- sp_dep %>% filter(species %in% sp_include)

unq_cl_elev <- sp_dep_filt %>% dplyr::select(cell,elev) %>% unique()

grid_elev_bands <- merge(grid_elev_bands,unq_cl_elev,by=c("elev","cell"))
elev_cell_centroids <- st_centroid(grid_elev_bands)
coordinates <- st_coordinates(elev_cell_centroids)

cell_elev_points <- data.frame(cell = elev_cell_centroids$cell,
                               elev=elev_cell_centroids$elev,
                               lon = coordinates[,1],
                               lat = coordinates[,2])

poly_shifted <- vect(grid_elev_bands)
poly_shifted <- terra::rotate(poly_shifted,longitude=0,left=FALSE,split=TRUE)
poly_shifted$elev <- grid_elev_bands$elev
poly_shifted$cell <- grid_elev_bands$cell

pdf("pg_era5_check.pdf")
plot(poly_shifted)
dev.off()

full_temp_list <- list()
count <- 0
for (year_in in 2001:2024){
  print(year_in)
  count <- count + 1
  year_rec <- c()
  comb_rast <- rast(paste("data/input/environmental/ERA5/temp_",year_in,"_",1,".nc",sep=""))
  for(month_in in 2:12){
    month_in_rast <- rast(paste("data/input/environmental/ERA5/temp_",year_in,"_",month_in,".nc",sep=""))
    comb_rast <- c(comb_rast,month_in_rast)
    
  }
  print(length(names(comb_rast)))
  names(comb_rast) <- 1:length(names(comb_rast))
  
  tictoc::tic()
  polygon_avg_temp <- terra::extract(comb_rast,poly_shifted,fun=mean,na.rm=TRUE)
  tictoc::toc()
  full_temp_list[[count]] <- polygon_avg_temp
}

saveRDS(full_temp_list,"data/input/environmental/intermediate_temp/polygon_global_temp.rds")


### Precip ####
full_precip_list <- list()
count <- 0
for (year_in in 2001:2024){
  print(year_in)
  count <- count + 1
  year_rec <- c()
  comb_rast <- rast(paste("data/input/environmental/ERA5/precip_",year_in,"_",1,".nc",sep=""))
  for(month_in in 2:12){
    month_in_rast <- rast(paste("data/input/environmental/ERA5/precip_",year_in,"_",month_in,".nc",sep=""))
    comb_rast <- c(comb_rast,month_in_rast)
    
  }
  print(length(names(comb_rast)))
  names(comb_rast) <- 1:length(names(comb_rast))
  
  tictoc::tic()
  polygon_avg_precip <- terra::extract(comb_rast,poly_shifted,fun=mean,na.rm=TRUE)
  tictoc::toc()
  full_precip_list[[count]] <- polygon_avg_precip
}

saveRDS(full_precip_list,"data/input/environmental/intermediate_precip/polygon_global_precip.rds")


### Calculate rolling metrics ####
full_temp_list <- readRDS("data/input/environmental/intermediate_temp/polygon_global_temp.rds")
full_precip_list <- readRDS("data/input/environmental/intermediate_precip/polygon_global_precip.rds")
number_of_sites <- nrow(full_temp_list[[1]])

#Temp
ID_full_ts <- list()
for(ID_in in 1:nrow(poly_shifted)){
  print(ID_in)
  #Create site-specific dataframe
  df_site_rec <- c()
  for (ii in 1:24){
    site_yr_temp <- full_temp_list[[ii]][ID_in,2:ncol(full_temp_list[[ii]])]
    site_yr_precip <- full_precip_list[[ii]][ID_in,2:ncol(full_precip_list[[ii]])]
    
    df_site <- data.frame(site = ID_in,
                          cell = poly_shifted$cell[ID_in],
                          elev = poly_shifted$elev[ID_in],
                          year=2000+ii,
                          day = 1:ncol(site_yr_precip),
                          temp = as.numeric(site_yr_temp),
                          precip = as.numeric(site_yr_precip))
    
    df_site_rec <- rbind(df_site_rec,df_site)
  }
  #Create rolling averages

  df_site_rec$rm_temp60 <- rollmean(df_site_rec$temp,k=60,align = "right",fill=NA)

  df_site_rec$rm_precip60 <- rollmean(df_site_rec$precip,k=60,align = "right",fill=NA)

  ID_full_ts[[ID_in]] <- df_site_rec
}

saveRDS(ID_full_ts,"data/input/environmental/enviro_data_ts_60.rds")

