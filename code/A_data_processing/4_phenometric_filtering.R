# Code to compile all GAM fits, then filter out bad estimates
library(dplyr)
library(tidyr)
library(ggplot2)
library(dggridR)
library(sf)
library(zoo)


### Compile report cards, GAM fits ####
#Get all "report cards" combined into one
all_rcs <- lapply(Sys.glob("data/input/phenology/report_cards/*.rds"), readRDS)
compiled_rcs <- do.call(bind_rows,all_rcs)

print(nrow(compiled_rcs))

#Bring in all species, year files
compiled_yrs <- c()
rc_all_years <- c()
for (each_year in 2002:2024){
  print(each_year)
  all_yr_fits <- lapply(Sys.glob(paste("data/input/phenology/phenometrics_cell_yr_files/",each_year,"/*.rds",sep="")), readRDS)
  compiled_yrs <- rbind(compiled_yrs,do.call(bind_rows,all_yr_fits))
  all_yr_fits <- c()
}

#Filter out bad fits
compiled_yrs <- filter(compiled_yrs,max_Rhat < 1.02 & min_neff >= 400 & 
                         mlmax == TRUE & plmax >=.99 & num_diverge == 0)

### Spring departure ####

s_dep <- compiled_yrs %>% filter(season=="spring" & arr_or_dep == "departure" & trailing_hm_sd <= 15)
nrow(s_dep)

#Filter out late fits, potentially problematic fall departures
s_dep <- s_dep %>% filter(arr_max_mean < 150)
nrow(s_dep)

hist(s_dep$trailing_hm_mean)

s_dep_m <- s_dep %>% select(c(species,cell,elev,year,trailing_hm_mean,trailing_hm_sd,arr_max_mean,arr_max_sd))
s_dep_m$likely_migratory = NA

### Run through additional filtering steps ####

grid_elev_bands <- st_read("data/input/spatial/grid_elev_bands.shp")
names(grid_elev_bands)[1] <- "elev"

elev_cell_centroids <- st_centroid(grid_elev_bands)
coordinates <- st_coordinates(elev_cell_centroids)

cell_elev_points <- data.frame(cell = elev_cell_centroids$cell,
                               elev=elev_cell_centroids$elev,
                               lon = coordinates[,1],
                               lat = coordinates[,2])

sp_dep <- merge(s_dep_m,cell_elev_points,by=c("cell","elev"))

# Filter to species-sites with at least 3 data points,
# And sites where the average error is less than 10 ()
sp_dep <- sp_dep %>% 
  group_by(species,cell,elev) %>% 
  mutate(sp_site_n = n(),
         mean_obs_err = mean(trailing_hm_sd)) %>%
  filter(sp_site_n >= 3 & mean_obs_err < 10) %>%
  ungroup()
nrow(sp_dep)
#Filter to species with at least 3 sites
unq_site_by_sp <- sp_dep %>% 
  select(species,cell,elev) %>%
  unique() %>%
  group_by(species) %>%
  summarise(unq_sp_site_count = n())

sp_dep <- merge(sp_dep,unq_site_by_sp,by="species")

#Get total data points for each species
#Get average sd on phenometric
sp_dep <- sp_dep %>% 
  group_by(species) %>% 
  mutate(sp_n = n()) %>%
  ungroup()

sp_dep_filt <- sp_dep %>% 
  filter(unq_sp_site_count >= 2 & sp_n >= 20) 

nrow(sp_dep_filt)
unique(sp_dep_filt$species)

#Filter out breeders that move south
south_breeders <- c("Brown-chested Martin","Double-collared Seedeater","Glittering-bellied Emerald",
                    "Lined Seedeater","Sick's Swift","White-crested Elaenia","Dollarbird","Greater Striped Swallow",
                    "Rainbow Bee-eater","Australian Reed Warbler","Australian Rufous Fantail","Black-faced Monarch",
                    "Blue-and-white Swallow","Bran-colored Flycatcher","Channel-billed Cuckoo","Chilean Swallow",
                    "Chivi Vireo","Dideric Cuckoo","Dusky Woodswallow","Fairy Martin","Fork-tailed Flycatcher",
                    "Gray-breasted Martin","Pacific Koel","Sacred Kingfisher","Sahul Brush Cuckoo",
                    "Sahul Cicadabird","Sooty Shearwater","Southern Martin","Streaked Flycatcher",
                    "Torresian Imperial-Pigeon","Tropical Kingbird","Variegated Flycatcher","Vermilion Flycatcher",
                    "White-winged Becard")

sp_dep_all_filts <- sp_dep_filt %>% filter(!(species %in% south_breeders))

nrow(sp_dep_all_filts)

saveRDS(sp_dep_all_filts,"data/input/phenology/spr_dep_core.rds")
