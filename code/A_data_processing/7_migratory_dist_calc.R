# Calculate migratory distances for the species in the analysis
library(dplyr)

sp_mig_data <- readRDS("data/input/mdl_data/full_data_ss.rds") %>% select(species) %>% unique()

mig_dist_bl <- readRDS("data/input/mig_dist/mig-dis-2025-11-26.rds")
sp_names <- read.csv("data/input/mig_dist/sp_names.csv")
colnames(sp_names) <- c("sci_name","species")

sp_mig_data <- merge(sp_mig_data,sp_names,by="species")

#Change species name conflicts
sp_mig_data$sci_name[which(sp_mig_data$species == "Evening Grosbeak")] <- "Hesperiphona vespertina"
sp_mig_data$sci_name[which(sp_mig_data$species == "Tundra Bean-Goose")] <- "Anser fabalis"
sp_mig_data$sci_name[which(sp_mig_data$species == "Sandhill Crane")] <- "Grus canadensis"
sp_mig_data$sci_name[which(sp_mig_data$species == "Pallas's Gull")] <- "Larus ichthyaetus"
sp_mig_data$sci_name[which(sp_mig_data$species == "Short-billed Gull")] <- "Larus canus"
sp_mig_data$sci_name[which(sp_mig_data$species == "Franklin's Gull")] <- "Larus pipixcan"
sp_mig_data$sci_name[which(sp_mig_data$species == "Double-crested Cormorant")] <- "Nannopterum auritus"
sp_mig_data$sci_name[which(sp_mig_data$species == "American Tree Sparrow")] <- "Passerella arborea"

sp_mig_data_all <- merge(sp_mig_data,mig_dist_bl,by="sci_name",all.x=TRUE)

saveRDS(sp_mig_data_all,"data/output/mig_dist/sp_mig_data_all.rds")
