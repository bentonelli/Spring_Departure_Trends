# Code to look at projections for departure dates

library(dplyr)
library(MCMCvis)
library(ggplot2)
library(tidyr)

#Model fit and data
mdl_fit <- readRDS("data/output/departure_model/departure_climate_fit_ss2.rds")

alphas <- rstan::extract(mdl_fit,pars=c("alpha"),permuted=TRUE)
alphas_df <- as.matrix(as.data.frame(alphas))
mean_dep_2019 <- MCMCsummary(alphas_df,pg0 = TRUE,probs = c(0.055,.5,0.945))

b1s <- rstan::extract(mdl_fit,pars=c("beta1"),permuted=TRUE)
b1s_df <- as.matrix(as.data.frame(b1s))
beta1s <- MCMCsummary(b1s_df,pg0 = TRUE,probs = c(0.055,.5,0.945))

b2s <- rstan::extract(mdl_fit,pars=c("beta2"),permuted=TRUE)
b2s_df <- as.matrix(as.data.frame(b2s))
beta2s <- MCMCsummary(b2s_df,pg0 = TRUE,probs = c(0.055,.5,0.945))

mdl_fit <- NULL # Clear memory
gc()

mdl_data <- readRDS("data/input/mdl_data/full_data_ss.rds")
indexing_tab <- readRDS("data/input/mdl_data/indexing_tab_data_ss.rds")

#Get data from ERA5 used for fitting model
era_climate_data <- readRDS("data/input/environmental/enviro_data_ts.rds")
era_climate_long <- bind_rows(era_climate_data)
era_climate_long <- era_climate_long %>% select(site,cell,elev,year,day,rm_temp60,rm_precip60)

#Get all the cell, elev, species, mean_dep
unq_site_sp <- mdl_data %>% 
  select(c("sp_site","cell","elev","species","mean_dep")) %>% 
  unique()

for (clim_proj_model_name in c("CESM2","CNRM-ESM2-1","MIROC6")){
  for(clim_proj_ssp_name in c("SSP1_2p6","SSP2_4p5","SSP3_7p0")){
    
    print(clim_proj_model_name)
    print(clim_proj_ssp_name)
    
    #Pick climate model run to use - for testing
    #clim_proj_model_name <- "CNRM-ESM2-1"
    #clim_proj_ssp_name <- "SSP2_4p5"
    
    clim_proj_temp <- readRDS(paste("data/input/environmental/climate_projections/",clim_proj_model_name,"/temp/",
                                    clim_proj_ssp_name,"/temp_",clim_proj_model_name,"_",clim_proj_ssp_name,".rds",sep=""))
    
    clim_proj_precip <- readRDS(paste("data/input/environmental/climate_projections/",clim_proj_model_name,"/precip/",
                                      clim_proj_ssp_name,"/precip_",clim_proj_model_name,"_",clim_proj_ssp_name,".rds",sep=""))
    #Put into long format
    clim_proj_temp_long <- bind_rows(clim_proj_temp)
    clim_proj_precip_long <- bind_rows(clim_proj_precip)
    
    clim_proj_temp <- clim_proj_precip <- NULL
    
    site_skeleton <- mdl_data %>% select(cell,elev,species,mean_dep) %>% unique()
    
    #Look at how climate model temperature estimates vary between CMIP6 and ERA5 Land,
    clim_proj_temp_long_hist <- clim_proj_temp_long %>% 
      filter(year %in% 2016:2024)
    
    clim_proj_precip_long_hist <- clim_proj_precip_long %>% 
      filter(year %in% 2016:2024)
    
    clim_historical <- era_climate_long %>% 
      filter(year %in% 2016:2024)
    
    site_clim <- merge(site_skeleton,clim_historical,
                       by.x=c("cell","elev","mean_dep"),
                       by.y =c("cell","elev","day"))
    
    site_clim <- merge(site_clim,clim_proj_temp_long_hist,
                       by.x=c("cell","elev","year","mean_dep"),
                       by.y =c("cell","elev","year","jday"))
    
    site_clim <- merge(site_clim,clim_proj_precip_long_hist,
                       by.x=c("cell","elev","year","mean_dep"),
                       by.y =c("cell","elev","year","jday"))
    
    plot(site_clim$rm_temp60,site_clim$mean_temp60,cex=.01)
    cor(site_clim$rm_temp60,site_clim$mean_temp60)
    plot(site_clim$rm_precip60,site_clim$total_precip60,cex=.01)
    cor(site_clim$rm_precip60,site_clim$total_precip60)
    
    # Data for differences at population level
    mdl_matching_diff <- site_clim %>% 
      group_by(cell,elev,species) %>% 
      summarise(mean_diff_temp = mean(mean_temp60-rm_temp60),
                mean_diff_precip = mean(total_precip60-rm_precip60))
    
    sp_it <- mdl_data %>% 
      select(c("sp_site","cell","elev",
               "lat","lon",
               "species","mean_dep",
               "mean_temp60","sd_temp_60",
               "mean_total_precip60","sd_total_precip60")) %>%
      unique()
    
    
    sp_it <- merge(sp_it,mdl_matching_diff,by=c("cell","elev","species")) %>% 
      arrange(sp_site)
    
    sp_it$a_sp <- mean_dep_2019$mean
    sp_it$b1_sp <- beta1s$mean
    sp_it$b2_sp <- beta2s$mean
    
    #Estimate future expectations for departure
    proj_temp_all <- c()
    proj_precip_all <- c()
    count <- 0
    for(each_year in 2015:2100){
      print(each_year)
      count <- count + 1
      yr_clim_proj_temp <- clim_proj_temp_long %>% filter(year == each_year) %>% select(-c("year"))
      yr_clim_proj_precip <- clim_proj_precip_long %>% filter(year == each_year) %>% select(-c("year"))
      colnames(yr_clim_proj_temp)[4] <- paste("proj_temp60_",each_year,sep="")
      colnames(yr_clim_proj_precip)[4] <- paste("proj_precip60_",each_year,sep="")
      
      proj_temp <- merge(sp_it,yr_clim_proj_temp,
                         by.x=c("cell","elev","mean_dep"),
                         by.y=c("cell","elev","jday"),
                         all.x=TRUE) %>%
        arrange(sp_site)
      
      proj_precip <- merge(sp_it,yr_clim_proj_precip,
                           by.x=c("cell","elev","mean_dep"),
                           by.y=c("cell","elev","jday"),
                           all.x=TRUE) %>%
        arrange(sp_site)
      
      proj_temp_all <- cbind(proj_temp_all,proj_temp[,ncol(proj_temp)])
      proj_precip_all <- cbind(proj_precip_all,proj_precip[,ncol(proj_precip)])
    }
    
    
    temp_change_abs <- (proj_temp_all-sp_it$mean_diff_temp)- sp_it$mean_temp60
    temp_anom_cor <- ((proj_temp_all-sp_it$mean_diff_temp) - sp_it$mean_temp60)/sp_it$sd_temp_60

    precip_change_abs <- (proj_precip_all-sp_it$mean_diff_precip) - sp_it$mean_total_precip60
    precip_anom_cor <- ((proj_precip_all-sp_it$mean_diff_precip) - sp_it$mean_total_precip60)/sp_it$sd_total_precip60

    year_mat <- c()
    count <- -4
    for (nn in 1:ncol(temp_anom_cor)){
      year_mat <- cbind(year_mat,rep(count,nrow(sp_it)))
      count <- count + 1
    }
    
    change_dep <- sp_it$b1_sp * temp_anom_cor + sp_it$b2_sp * precip_anom_cor
    change_dep <- as.data.frame(change_dep)
    
    saveRDS(as.data.frame(temp_change_abs),paste("data/output/depature_projections/temp_abs_change_",clim_proj_model_name,"_",clim_proj_ssp_name,".rds",sep=""))
    saveRDS(as.data.frame(precip_change_abs),paste("data/output/depature_projections/precip_abs_change_",clim_proj_model_name,"_",clim_proj_ssp_name,".rds",sep=""))
    saveRDS(change_dep,paste("data/output/depature_projections/est_dep_",clim_proj_model_name,"_",clim_proj_ssp_name,".rds",sep=""))
  }
}
