# Look at climate correlations between:

# 1. May temperatures on the breeding grounds
# 2. Departure period temperature

#AND

# 1. May temperatures on the breeding grounds
# 2. Departure period *precipitation*

### PACKAGES ####
library(dplyr)
library(dggridR)
library(terra)
library(exactextractr)
library(tidyr)
library(sf)
library(geosphere)
library(readr)

#Set enviro_var of interest (on NB grounds) - Can change to "temp" or "precip"
enviro_var <- "precip"

### READ DATA ####
#Species names
bird_names <- read_csv("~/Documents/Coding/R/Main/Avian_Pheno_Database/data/sci_name_birdlife_merged.csv")

may_enviro <- rast("data/input/environmental/temp_may_global.nc")
may_enviro <- project(may_enviro,"EPSG:4326")

enviro_ts <- readRDS("data/input/environmental/enviro_data_ts_60.rds")
era_climate_long <- bind_rows(enviro_ts)

mdl_data <- readRDS("data/input/mdl_data/full_data_ss.rds")
unq_area_days <- mdl_data %>% select(cell,elev,species,mean_dep) %>% unique()

### Set up grid system ####
hexgrid6 <- dggridR::dgconstruct(res = 6)
grid_poly <- dgrectgrid(hexgrid6,minlat = -91,minlon = -181,maxlat = 91,maxlon = 181)
grid_poly <- sf::st_wrap_dateline(grid_poly, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)

# Get breeding range of species
sp_targ_list <- unique(mdl_data$species)

sp_corr_recs <- c()
all_corrs_list <- list()
ll <- 0
for(spec_target in sp_targ_list){
  print(spec_target)
  #spec_target <- "Gray Catbird"
  
  sci_name_target <- bird_names$sci_name_birdlife[which(bird_names$species == spec_target)]
  
  #Breeding polygons, avialable via BirdLife International
  bird_br <- vect(paste0("data/input/spatial/breeding_polygons/",gsub(" ","_",sci_name_target),".shp"))
  
  bird_br <- aggregate(bird_br)
  
  bird_br <- project(bird_br, crs(grid_poly))#project(bird_br,crs(grid_poly))
  
  fixed_grid <- makeValid(vect(grid_poly))
  
  cells_clipped <- intersect(fixed_grid, bird_br)
  cells_clipped$area_km2 <- expanse(cells_clipped, unit = "km")
  
  cells_clipped <- cells_clipped[which(cells_clipped$area_km2 >= 625)]
  
  #plot(cells_clipped)
  
  yr_vals <- terra::extract(may_enviro,cells_clipped,fun="mean",na.rm=TRUE)
  
  colnames(yr_vals)[2:26] <- 2001:2025
  yr_vals$ID <- cells_clipped$seqnum
  
  breed_temps_flipped <- t(yr_vals[,2:25])
  
  spec_unq_areas <- unq_area_days %>% filter(species == spec_target)
  spec_unq_areas$max_r_sq <- NA
  spec_unq_areas$up95_r_sq <- NA
  spec_unq_areas$mean_r_sq <- NA
  
  spec_unq_areas$max_r <- NA
  spec_unq_areas$mean_r <- NA
  
  #To speed things up
  era_climate_sp <- era_climate_long %>% filter(cell %in% spec_unq_areas$cell)
  
  for(mm in 1:nrow(spec_unq_areas)){
    ll <- ll + 1
    nb_ts <- era_climate_sp %>% filter(cell == spec_unq_areas$cell[mm] & 
                                           elev == spec_unq_areas$elev[mm] & 
                                           day == spec_unq_areas$mean_dep[mm])
    if(enviro_var == "temp"){
      all_corrs <- cor(nb_ts$rm_temp60,breed_temps_flipped,use = "pairwise.complete.obs")
    } else if(enviro_var == "precip"){
      all_corrs <- cor(nb_ts$rm_precip60,breed_temps_flipped,use = "pairwise.complete.obs")
      if(sum(nb_ts$rm_precip60,na.rm=TRUE)==0){
        stop()
      }
    } else {
      stop()
    }
    
    spec_unq_areas$max_r_sq[mm] <- max(all_corrs^2)
    spec_unq_areas$up95_r_sq[mm] <- quantile(all_corrs^2,0.95)
    spec_unq_areas$mean_r_sq[mm] <- mean(all_corrs^2)
    
    spec_unq_areas$max_r[mm] <- all_corrs[min(which(all_corrs^2 == spec_unq_areas$max_r_sq[mm]))]
    spec_unq_areas$mean_r[mm] <- mean(all_corrs)
    
    all_corrs_list[[ll]] <- all_corrs
  }
  sp_corr_recs <- rbind(sp_corr_recs,spec_unq_areas)
}

sp_summary <- sp_corr_recs %>% 
  group_by(species) %>% 
  summarise(mean_max_r_sq = mean(max_r_sq),
            mean_95_r_sq = mean(up95_r_sq),
            mean_mean_r_sq = mean(mean_r_sq),
            count_areas = n())
sp_summary

# Save intermediate step 
saveRDS(sp_summary,paste0("data/output/variance_explained/sp_summary_",enviro_var,".rds"))
saveRDS(all_corrs_list,paste0("data/output/variance_explained/all_corrs_list_",enviro_var,".rds"))

all_corrs_list_t <- readRDS("data/output/variance_explained/all_corrs_list_temp.rds")
all_corrs_list_p <- readRDS("data/output/variance_explained/all_corrs_list_precip.rds")

temp_max_all_pops <- unlist(sapply(all_corrs_list_t,function(x) max(abs(x), na.rm = TRUE)))
precip_max_all_pops <- unlist(sapply(all_corrs_list_p,function(x) max(abs(x), na.rm = TRUE)))

site_data <- readRDS("data/input/mdl_data/indexing_tab_data_ss.rds")

site_data$best_corr <- temp_max_all_pops
site_data$best_corr[which(temp_max_all_pops < precip_max_all_pops)] <- precip_max_all_pops[which(temp_max_all_pops < precip_max_all_pops)] 

sp_max_corr <- site_data %>% group_by(sp) %>%
  summarise(mean_best_R2 = mean(best_corr^2))

saveRDS(sp_max_corr,"data/output/variance_explained/site_data_best_corr.rds")

