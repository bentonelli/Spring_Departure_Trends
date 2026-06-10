library(dplyr)
library(readr)
library(rstan)
library(rstanarm)
library(ggplot2)
library(dggridR)
library(sf)
library(terra)
library(MCMCvis)
library(mapdata)
library(birdcolors)

### Subplot 1 - Departure GAMs ####

# Code to read in zero-filled data, run GAMs by location, year.

# If testing, uncomment the below
# #Set cell to test on
year_in <- 2023
cell_in <- 702
spec_in <- "White-throated Sparrow"

thinning_level <- 5000
save_dfits <- FALSE
prior_smooth_val <- 1
targ_prop_pos <- 0.10

source("code/X_data_processing/0_phenometric_functions.R")

report_cards <- readRDS(paste("data/input/phenology/report_cards/rc_",year_in,".rds",sep="")) %>%
  filter(cell %in% cell_in & accepted == 1) %>% arrange(cell,spec)

report_cards <- filter(report_cards,spec == spec_in & season == "spring" & arr_or_dep == "departure")

spec_matrix_dt <- readRDS(paste("data/input/phenology/zf_files/",year_in,"/zf_dt_",year_in,"_",cell_in,".rds",sep=""))

spec_list <- unique(report_cards$spec)
Sys.time()
if(nrow(report_cards) > 0){
  
  print(paste("Processing: ",nrow(report_cards)," Metrics",sep=""))
  print(paste("Representing: ",length(unique(report_cards$spec))," species",sep=""))
  
  #Set up stan model conditions
  DELTA <- .95
  TREE_DEPTH <- 10
  ITER <- 1500
  CHAINS <- 4
  
  phenometric_df <- data.frame(species = report_cards$spec, 
                               year = report_cards$year, 
                               cell = report_cards$cell,
                               
                               num_og_records = NA,
                               num_og_pos = NA,
                               
                               num_used_records = NA,
                               num_used_pos = NA,
                               
                               season = NA,
                               arr_or_dep = NA,
                               max_Rhat = NA,
                               min_neff = NA,
                               mlmax = NA,
                               plmax = NA,
                               num_diverge = NA,
                               num_tree = NA,
                               num_BFMI = NA,
                               delta = NA,
                               tree_depth = NA,
                               t_iter = NA,
                               elev = NA,
                               
                               avg_over_max = NA,
                               lmin_before_max = NA,
                               
                               lmin_after_max = NA,
                               min_depth_over_max = NA,
                               
                               leading_hm_mean = NA,
                               leading_hm_sd = NA,
                               trailing_hm_mean = NA,
                               trailing_hm_sd = NA,
                               
                               arr_diff_hm_mean = NA,
                               arr_diff_hm_sd = NA,
                               
                               arr_max_mean = NA,
                               arr_max_sd = NA)
  counter <- 0
  for(n_spec in spec_list){
    gc()
    print(n_spec)
    
    cl_sp_data <- spec_matrix_dt %>% 
      dplyr::select(c(1,2,3,4,6,which(colnames(spec_matrix_dt)==n_spec)))
    
    
    colnames(cl_sp_data)[6] <- "detect"
    og_number_rec <- nrow(cl_sp_data)
    pos_number_rec <- sum(cl_sp_data[6])
    
    #Thin here - this drastically reduces runtimes for high volume cell-year combos
    
    metric_targets <- report_cards %>% filter(spec == n_spec)
    cl_sp_data <- filter(cl_sp_data, elev_factor %in% unique(metric_targets$elev_band))
    
    cl_sp_data_og <- cl_sp_data
    
    cl_sp_data_band_filt <- c()
    for (ee in unique(cl_sp_data$elev_factor)){
      
      #Get data from elev band
      cl_sp_data_eb <- filter(cl_sp_data,elev_factor == ee)
      band_rec_num_og <- nrow(cl_sp_data_eb)
      pos_rec_og <- sum(cl_sp_data_eb[,6])
      
      targ_records <- thinning_level
      og_prop_pos <- sum(cl_sp_data_eb[,6])/band_rec_num_og
      
      
      # If the proportion of positive records is less than the target
      # AND there are more than the total # of target records
      if(og_prop_pos < targ_prop_pos & band_rec_num_og > targ_records){
        
        # Ideal number of positive records
        pos_rec_target <- floor(targ_prop_pos*targ_records)
        
        #Sometimes there aren't enough records to make that work though
        # If there are not:
        if (pos_rec_target >= sum(cl_sp_data_eb[,6])){
          #Use all positive records available
          cl_sp_data_eb_pos <- cl_sp_data_eb[which(cl_sp_data_eb[,6]==1),]
          
          #Then use the number available to caluclate the number of negative samples to subsample
          neg_rec_target <- nrow(cl_sp_data_eb_pos)/targ_prop_pos - nrow(cl_sp_data_eb_pos)
          cl_sp_data_eb_neg <- cl_sp_data_eb[sample(which(cl_sp_data_eb[,6]==0),neg_rec_target),]
          
          #If there are enough records, subsample to the pos_rec_target
        } else {
          cl_sp_data_eb_pos <- cl_sp_data_eb[sample(which(cl_sp_data_eb[,6]==1),pos_rec_target),]
          
          #Then use that number to get the number of negative checklists
          neg_rec_target <- targ_records - pos_rec_target
          cl_sp_data_eb_neg <- cl_sp_data_eb[sample(which(cl_sp_data_eb[,6]==0),neg_rec_target),]
        }
        
        #Combine
        cl_sp_data_eb <- rbind(cl_sp_data_eb_pos,cl_sp_data_eb_neg)
        
        #When there are less than X,XXX records, use all postive records available
      } else if (og_prop_pos < targ_prop_pos){
        pos_rec_target <- sum(cl_sp_data_eb[,6])
        neg_rec_target <- floor(pos_rec_target/targ_prop_pos - pos_rec_target) 
        
        #Include all positive records
        cl_sp_data_eb_pos <- cl_sp_data_eb[which(cl_sp_data_eb[,6]==1),]
        
        cl_sp_data_eb_neg <- cl_sp_data_eb[sample(which(cl_sp_data_eb[,6]==0),neg_rec_target),]
        
        cl_sp_data_eb <- rbind(cl_sp_data_eb_pos,cl_sp_data_eb_neg)
        
        #If more postive records than target, 
      } else {
        # And less # records than target
        if(band_rec_num_og < targ_records){
          #Use all data
          cl_sp_data_eb <- cl_sp_data_eb
          
          # Otherwise subsample without preference to postive/negative
        } else {
          cl_sp_data_eb <- cl_sp_data_eb[sample(1:nrow(cl_sp_data_eb),targ_records),]
        }
        
      }
      
      #Print some useful stuff out
      if (nrow(cl_sp_data_eb) < band_rec_num_og){
        reduction <- (1-(nrow(cl_sp_data_eb)/band_rec_num_og))*100
        print(paste("Thinned to",nrow(cl_sp_data_eb),"records,","a",round(reduction,2),"% reduction"))
        print(paste("Total 1s:",sum(cl_sp_data_eb[,6])))
      } else {
        print(paste("Elev. Level:",as.numeric(ee),"not thinned"))
      }
      
      cl_sp_data_eb$shr <- scale(cl_sp_data_eb$EH, scale = FALSE)[,1]
      cl_sp_data_band_filt <- rbind(cl_sp_data_band_filt,cl_sp_data_eb)
    }
    
    cl_sp_data <- cl_sp_data_band_filt %>% arrange(elev_factor,julian_day)
    
    elev_level_to_run <- unique(cl_sp_data$elev_factor)
    
    for (each_elev in elev_level_to_run){
      
      DELTA <- 0.95
      TREE_DEPTH <- 15
      
      print(paste("Running elev. Level:",as.numeric(each_elev)))
      #Elevation data split
      elev_split <- cl_sp_data %>% filter(elev_factor == each_elev)
      
      fit2 <- rstanarm::stan_gamm4(detect ~ s(julian_day, k = 30) + shr,
                                   data = elev_split,
                                   family = binomial(link = "logit"),
                                   algorithm = 'sampling',
                                   prior_smooth = exponential(prior_smooth_val),
                                   iter = ITER,
                                   chains = CHAINS,
                                   cores = CHAINS,
                                   adapt_delta = DELTA,
                                   control = list(max_treedepth = TREE_DEPTH))
      
      
      #Save model output
      num_diverge <- rstan::get_num_divergent(fit2$stanfit)
      num_tree <- rstan::get_num_max_treedepth(fit2$stanfit)
      num_BFMI <- length(rstan::get_low_bfmi_chains(fit2$stanfit))
      max_Rhat <- round(max(summary(fit2)[, 'Rhat']), 2)
      min_neff <- min(summary(fit2)[, 'n_eff'])
      
      #If there are divergences, try, try again
      while (num_diverge > 0 & DELTA <= 0.98){
        DELTA <- DELTA + 0.01
        
        fit2 <- rstanarm::stan_gamm4(detect ~ s(julian_day, k = 30) + shr,
                                     data = elev_split,
                                     family = binomial(link = "logit"),
                                     algorithm = 'sampling',
                                     prior_smooth = exponential(prior_smooth_val),
                                     iter = ITER,
                                     chains = CHAINS,
                                     cores = CHAINS,
                                     adapt_delta = DELTA,
                                     control = list(max_treedepth = TREE_DEPTH))
        
        #Save model output
        num_diverge <- rstan::get_num_divergent(fit2$stanfit)
        num_tree <- rstan::get_num_max_treedepth(fit2$stanfit)
        num_BFMI <- length(rstan::get_low_bfmi_chains(fit2$stanfit))
        max_Rhat <- round(max(summary(fit2)[, 'Rhat']), 2)
        min_neff <- min(summary(fit2)[, 'n_eff'])
      }
      
      
      #Testing: Add predictions for each elevation level in the data
      #Only do this when species detections are >= 30 in an elevational band
      if (save_dfits){
        elev_bands <- elev_split %>% filter(detect == 1) %>% 
          group_by(elev_factor) %>% 
          count() %>% 
          filter(n >= 30) %>% 
          select(elev_factor) %>% 
          unlist() %>%
          as.numeric()
        
        #elev_bands <- unique(elev_split$elev_factor)
        predictDays <- range(elev_split$julian_day)[1]:range(elev_split$julian_day)[2]
        dfit_list <- vector(mode="list",length = 7)
        
        rand_iter_sample <- sample(1:3000,500,replace=FALSE)
        for (each_elev_band in elev_bands){
          if (length(table(elev_split$elev_factor)) > 1){
            #Want to get arrival at most common elevation
            newdata <- data.frame(julian_day = predictDays, shr = 0,elev_factor=each_elev_band)
          } else {
            #If elevation is only in one bucket:
            newdata <- data.frame(julian_day = predictDays, shr = 0)
          }
          dfit_list[[as.numeric(each_elev_band)+1]] <- rstanarm::posterior_linpred(fit2, newdata = newdata, transform = T)[rand_iter_sample,]
        }
        spec_name_gsub <- gsub(" ","_",n_spec)
        # If the user wants the whole GAM fit
        #saveRDS(dfit_list,paste("data/dfits/",year_in,"/",spec_name_gsub,"_",cell_in,"_",year_in,".rds",sep=""))
      }
      
      # For each of the metrics needed, run the next chunk of code to extract
      # those estimates
      
      metric_targets <- report_cards %>% filter(spec == n_spec & elev_band == each_elev) 
      
      for(each_metric_target in 1:nrow(metric_targets)){
        counter <- counter + 1
        #Get season and arr/dep designation
        phenometric_df$elev[counter] <- as.numeric(each_elev)
        
        season_in <- metric_targets$season[each_metric_target]
        arr_or_dep_in <- metric_targets$arr_or_dep[each_metric_target]
        
        #Get time windows
        window_start <- season_time_windows$time_window_start[season_time_windows$season == season_in & 
                                                                season_time_windows$arr_or_dep == arr_or_dep_in]
        
        window_end <- season_time_windows$time_window_end[season_time_windows$season == season_in & 
                                                            season_time_windows$arr_or_dep == arr_or_dep_in]
        
        
        phenometric_df$num_og_records[counter] <- sum(cl_sp_data_og$elev_factor == each_elev)
        phenometric_df$num_og_pos[counter] <- sum(cl_sp_data_og$elev_factor == each_elev & cl_sp_data_og$detect == 1)
        
        phenometric_df$num_used_records[counter] <- sum(elev_split$elev_factor == each_elev)
        phenometric_df$num_used_pos[counter] <- sum(elev_split$elev_factor == each_elev & elev_split$detect == 1)
        
        phenometric_df$season[counter] <- season_in
        phenometric_df$arr_or_dep[counter] <- arr_or_dep_in
        
        phenometric_df$num_diverge[counter] <- num_diverge
        phenometric_df$num_tree[counter] <- num_tree
        phenometric_df$num_BFMI[counter] <- num_BFMI
        phenometric_df$delta[counter] <- DELTA
        phenometric_df$tree_depth[counter] <- TREE_DEPTH
        phenometric_df$t_iter[counter] <- ITER
        phenometric_df$max_Rhat[counter] <- max_Rhat
        phenometric_df$min_neff[counter] <- min_neff
        
        #generate predict data
        rng_sp_cl_data <- range(unique(elev_split$julian_day[which(elev_split$julian_day >= window_start & 
                                                                     elev_split$julian_day <= window_end)]))
        predictDays <- rng_sp_cl_data[1]:rng_sp_cl_data[2]
        
        # Get the saved elevation factor, use that to make predictions
        elev_mode <- metric_targets$elev_band[each_metric_target]
        
        newdata <- data.frame(julian_day = predictDays, shr = 0)
        
        #predict response
        dfit <- rstanarm::posterior_linpred(fit2, newdata = newdata, transform = T)
        leading_halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))
        trailing_halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))
        max_fit <- rep(NA, ((ITER/2)*CHAINS))
        tlmax <- rep(NA, ((ITER/2)*CHAINS))
        
        #Metrics for detecting primarily migratory cells
        num_eligible_max <- rep(NA, ((ITER/2)*CHAINS))
        avg_over_max <- rep(NA, ((ITER/2)*CHAINS))
        
        lmin_before_max <- rep(NA, ((ITER/2)*CHAINS))
        lmin_after_max <- rep(NA, ((ITER/2)*CHAINS))
        min_depth_over_max <- rep(NA, ((ITER/2)*CHAINS))
        
        
        ### ARRIVAL SPECIFIC METRICS ####
        if(arr_or_dep_in == "arrival"){
          #day at which probability of occurrence is half local maximum value
          for (L in 1:((ITER/2)*CHAINS)){
            rowL <- as.vector(dfit[L,])
            #first detection within prediction day range
            fd <- min(elev_split$julian_day[which(elev_split$detect == 1 & 
                                                    elev_split$julian_day >= rng_sp_cl_data[1] & 
                                                    elev_split$julian_day <= rng_sp_cl_data[2])])
            #local maximum(s)
            #from: stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
            lmax_idx <- which(diff(sign(diff(rowL))) == -2) + 1
            lmax <- predictDays[lmax_idx]
            #first local max to come after first detection
            
            if(length(lmax_idx)>1){
              max_peak_height <- max(rowL[lmax_idx])
            } else {
              max_peak_height <- max(rowL)
            }
            
            #First local max (that fits requirements)
            flm <- which(rowL[lmax_idx] >= mean(rowL) & 
                           rowL[lmax_idx] >= ((0.15)*max_peak_height) & 
                           lmax > fd)
            
            num_eligible_max[L] <- length(flm)
            if (length(flm) > 0)
            {
              #first local max to come after first detection
              lmax2_idx <- lmax_idx[min(flm)]
              lmax2 <- lmax[min(flm)]
              tlmax[L] <- TRUE
            } else {
              #no local max
              lmax2_idx <- which.max(rowL)
              lmax2 <- predictDays[which.max(rowL)]
              tlmax[L] <- FALSE
            }
            #store local max
            max_fit[L] <- lmax2
            
            #Save the depth of the minimum, if one exists
            lmin_idx_post_max <- min(c(which.min(rowL[lmax2_idx:length(rowL)]), 
                                       which(diff(sign(diff(rowL[lmax2_idx:length(rowL)]))) == 2) + 1)) + lmax2_idx - 1
            
            lmin_after_max[L] <- predictDays[lmin_idx_post_max] - lmax2
            min_depth_over_max[L] <- rowL[lmin_idx_post_max]/rowL[lmax2_idx]
            
            #local mins before max (global and local mins)
            lmin_idx <- c(which.min(rowL[1:lmax2_idx]), 
                          which(diff(sign(diff(rowL[1:lmax2_idx]))) == 2) + 1)
            lmin <- predictDays[lmin_idx]
            
            #local min nearest to local max
            lmin2_idx <- lmin_idx[which.min(lmax2 - lmin)]
            lmin2 <- predictDays[lmin2_idx]
            
            #value at local max - value at min (typically 0)
            leading_dmm <- rowL[lmax2_idx] - rowL[lmin2_idx]
            #all positions less than or equal to half diff between max and min + value min
            leading_tlm <- which(rowL <= ((leading_dmm/2) + rowL[lmin2_idx]))
            #which of these come before max and after or at min
            
            leading_vgm <- leading_tlm[which(leading_tlm < lmax2_idx & leading_tlm >= lmin2_idx)]
            
            #insert halfmax (first day for situations where max is a julian_day = first day)
            if (length(leading_vgm) > 0)
            {
              leading_halfmax_fit[L] <- predictDays[max(leading_vgm)]
            } else {
              leading_halfmax_fit[L] <- predictDays[1]
            }
            
            #Now get trailing hm as well
            lmin_idx_post_max
            
            #value diff from last maximum to absolute minimum
            trailing_dmm <- rowL[lmax2_idx] - rowL[lmin_idx_post_max]
            
            #all positions less than or equal to half diff between max and min + value min
            trailing_tlm <- which(rowL <= ((trailing_dmm/2) + rowL[lmin_idx_post_max]))
            
            #which of these come before min but after max
            trailing_vgm <- trailing_tlm[which(trailing_tlm < lmin_idx_post_max & trailing_tlm >= lmax2_idx)]
            
            #insert trailing_halfmax (first day for situations where max is a julian_day = 1)
            if (length(trailing_vgm) > 0)
            {
              trailing_halfmax_fit[L] <- predictDays[min(trailing_vgm)]
            } else {
              trailing_halfmax_fit[L] <- predictDays[1]
            }
            
          }
        } else if (arr_or_dep_in == "departure") {
          ### DEPARTURE SPECIFIC METRICS ####
          #day at which probability of occurrence is half local maximum value
          for (L in 1:((ITER/2)*CHAINS)){
            rowL <- as.vector(dfit[L,])
            #first detection
            #fd <- min(cell_data$julian_day[which(cell_data$detect == 1)])
            #local maximum(s)
            #from: stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
            lmax_idx <- which(diff(sign(diff(rowL))) == -2) + 1
            lmax <- predictDays[lmax_idx]
            
            abs_max <- max(which(rowL==max(rowL)))
            abs_max_val <-rowL[abs_max]
            
            #Last local max before absolute minimum
            eligible_max <- rowL[lmax_idx] >= (0.15 * abs_max_val) | rowL[lmax_idx] >= mean(rowL)
            num_eligible_max[L] <- length(eligible_max)
            
            #if (length(lmax_idx[which(rowL[lmax_idx] >= (0.15 * abs_max_val))]) > 0){
            if (length(lmax_idx[which(eligible_max)]) > 0){
              #last_eligible_max <- max(lmax_idx[which(rowL[lmax_idx] >= (0.15 * abs_max_val))])
              last_eligible_max <- max(lmax_idx[which(eligible_max)])
              after_last_max_idx <- rowL[last_eligible_max:length(rowL)]
            } else {
              last_eligible_max <- abs_max
              after_last_max_idx <- rowL[last_eligible_max:length(rowL)]
            }
            
            
            abs_min <- min(which(after_last_max_idx == min(after_last_max_idx)) + last_eligible_max-1)
            
            #Find last max before absolute minimum
            #Last max has to be at least X% of absolute max
            #last_max_idx <- which(lmax_idx < abs_min & rowL[lmax_idx] >= (0.15 * abs_max_val))
            last_max_idx <- which(lmax_idx < abs_min & rowL[lmax_idx] >= (0.15 * abs_max_val) | rowL[lmax_idx] >= mean(rowL))
            
            if (length(last_max_idx) > 0)
            {
              #last local max to come before absolute minimum
              last_max <- max(last_max_idx) 
              lmax2_idx <- lmax_idx[last_max]
              ##lmax2_idx <- lmax_idx[min(last_max)]
              lmax2 <- lmax[last_max]
              tlmax[L] <- TRUE
            } else {
              #no local max
              lmax2_idx <- which.max(rowL)
              lmax2 <- predictDays[which.max(rowL)]
              tlmax[L] <- FALSE
            }
            #store local max
            if(is.na(lmax2)){
              print(last_max)
              print(which.max(rowL))
              print(rowL)
            }
            max_fit[L] <- lmax2 
            
            #Save the average probability of detection before the max as a % of the max OVER THE LAST 90 DAYS
            if (lmax2_idx >= 91){
              day_window <- rowL[(lmax2_idx-90):lmax2_idx]
              lmin_before_max[L] <- predictDays[which.min(day_window)] + 90
            } else{
              #Unless lmax2_idx comes before 91 (likely bad fit), in which case start at 1
              day_window <- rowL[1:lmax2_idx]
              lmin_before_max[L] <- predictDays[which.min(day_window)]
            }
            
            avg_over_max[L] <- mean(day_window)/rowL[lmax2_idx]
            min_depth_over_max[L] <- day_window[which.min(day_window)]/rowL[lmax2_idx]
            
            #local mins before max (global and local mins)
            leading_lmin_idx <- c(which.min(rowL[1:lmax2_idx]), 
                                  which(diff(sign(diff(rowL[1:lmax2_idx]))) == 2) + 1)
            leading_lmin <- predictDays[leading_lmin_idx]
            
            #local min nearest to local max
            leading_lmin2_idx <- leading_lmin_idx[which.min(lmax2_idx - leading_lmin_idx)]
            leading_lmin2 <- predictDays[leading_lmin2_idx]
            
            #value at local max - value at min
            leading_dmm <- rowL[lmax2_idx] - rowL[leading_lmin2_idx]
            #all positions less than or equal to half diff between max and min + value min
            leading_tlm <- which(rowL <= ((leading_dmm/2) + rowL[leading_lmin2_idx]))
            
            leading_vgm <- leading_tlm[which(leading_tlm < lmax2_idx & leading_tlm >= leading_lmin2_idx)]
            
            if (length(leading_vgm) > 0)
            {
              leading_halfmax_fit[L] <- predictDays[max(leading_vgm)]
            } else {
              leading_halfmax_fit[L] <- predictDays[1]
            }
            
            
            #Global min
            lmin_idx <- abs_min
            
            trailing_lmin_idx <- c(which.min(rowL[lmax2_idx:length(rowL)]), 
                                   which(diff(sign(diff(rowL[lmax2_idx:length(rowL)]))) == 2) + 1)
            
            #lmin <- predictDays[lmin_idx]
            
            #value diff from last maximum to absolute minimum
            dmm <- rowL[lmax2_idx] - rowL[abs_min]
            
            #all positions less than or equal to half diff between max and min + value min
            tlm <- which(rowL <= ((dmm/2) + rowL[lmin_idx]))
            
            #which of these come before min but after max
            vgm <- tlm[which(tlm < lmin_idx & tlm >= lmax2_idx)]
            
            #insert trailing_halfmax (first day for situations where max is a julian_day = 1)
            if (length(vgm) > 0)
            {
              trailing_halfmax_fit[L] <- predictDays[min(vgm)]
            } else {
              trailing_halfmax_fit[L] <- predictDays[1]
            }
          }
        }
        
        ### Save results ####
        #number of iterations that had local max
        phenometric_df$plmax[counter] <- round(sum(tlmax)/((ITER/2)*CHAINS), 3)
        
        #model fit
        mn_dfit <- apply(dfit, 2, mean)
        LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
        UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
        
        #check whether local max exists for mean model fit
        mlmax <- sum(diff(sign(diff(mn_dfit))) == -2)
        if (mlmax > 0){
          phenometric_df$mlmax[counter] <- TRUE
        } else {
          phenometric_df$mlmax[counter] <- FALSE
        }
        
        #estimated leading halfmax
        mn_leading_hm <- mean(leading_halfmax_fit)
        LCI_leading_hm <- quantile(leading_halfmax_fit, probs = 0.025)
        UCI_leading_hm <- quantile(leading_halfmax_fit, probs = 0.975)
        
        #estimated trailing halfmax
        mn_trailing_hm <- mean(trailing_halfmax_fit)
        LCI_trailing_hm <- quantile(trailing_halfmax_fit, probs = 0.025)
        UCI_trailing_hm <- quantile(trailing_halfmax_fit, probs = 0.975)
        
        #estimated difference in halfmax estimates
        halfmax_diff <- trailing_halfmax_fit - leading_halfmax_fit
        mn_diff_hm <- mean(halfmax_diff)
        LCI_diff_hm <- quantile(halfmax_diff, probs = 0.025)
        UCI_diff_hm <- quantile(halfmax_diff, probs = 0.975)
        
        #estimated max
        mn_max <- mean(max_fit)
        LCI_max <- quantile(max_fit, probs = 0.025)
        UCI_max <- quantile(max_fit, probs = 0.975)
        
        
        phenometric_df[counter,"median_eligible_max"] <- median(num_eligible_max)
        phenometric_df[counter,"other_eligible_max_perc"] <- length(which(num_eligible_max == median(num_eligible_max)))/L
        
        phenometric_df[counter,"avg_over_max"] <- mean(avg_over_max)
        
        phenometric_df[counter,"lmin_before_max"] <- mean(lmin_before_max)
        phenometric_df[counter,"lmin_after_max"] <- mean(lmin_after_max)
        phenometric_df[counter,"min_depth_over_max"] <- mean(min_depth_over_max)
        
        phenometric_df[counter,"leading_hm_mean"] <- mean(leading_halfmax_fit)
        phenometric_df[counter,"leading_hm_sd"] <- sd(leading_halfmax_fit)
        
        phenometric_df[counter,"trailing_hm_mean"] <- mean(trailing_halfmax_fit)
        phenometric_df[counter,"trailing_hm_sd"] <- sd(trailing_halfmax_fit)
        
        phenometric_df[counter,"arr_diff_hm_mean"] <- mean(halfmax_diff)
        phenometric_df[counter,"arr_diff_hm_sd"] <- sd(halfmax_diff)
        
        phenometric_df[counter,"arr_max_mean"] <- mean(max_fit)
        phenometric_df[counter,"arr_max_sd"] <- sd(max_fit)
        
        ###Plot model fit and data #### 
        
        # Uncomment to save PDF summary 
        
        pdf(paste0("data/output/spec_hm_pdfs/",season_in,"_",arr_or_dep_in,"/",n_spec, "_",
                   year_in, "_",cell_in,"_",season_in ,"_",arr_or_dep_in,"_",elev_mode,".pdf"))
        plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 2,
             ylim = c(0, max(UCI_dfit)),
             main = paste0(n_spec, ' - ', year_in, ' - ', cell_in,' - ',as.character(elev_mode)),
             xlab = "Julian Day", ylab = "Probability of occurrence")
        lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 2)
        lines(predictDays, mn_dfit, lwd = 2)
        dd <- elev_split$detect[elev_split$elev_factor == elev_mode]
        dd[which(dd == 1)] <- max(UCI_dfit)
        points(elev_split$julian_day[elev_split$elev_factor == elev_mode], dd, col = rgb(0,0,0,0.25))

        abline(v = mn_trailing_hm, col = rgb(0,0,1,0.5), lwd = 2)
        abline(v = LCI_trailing_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
        abline(v = UCI_trailing_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
        
        abline(v = mn_max, col = rgb(0,1,0,0.5), lwd = 2)
        abline(v = LCI_max, col = rgb(0,1,0,0.5), lwd = 2, lty = 2)
        abline(v = UCI_max, col = rgb(0,1,0,0.5), lwd = 2, lty = 2)
        
        legend("bottomleft",
               legend = c("Model fit", "CI fit", "Trailing HM", "CI HM"),
               col = c('black', 'red', rgb(0,0,1,0.5), rgb(0,0,1,0.5)),
               lty = c(1,2,1,2), lwd = c(2,2,2,2), cex = 1.3)
        
        dev.off()
      }
      
    }
  }
}
Sys.time()


### Subplot 2 - Data density map ####

#Model fit and data
mdl_data <- readRDS("data/input/mdl_data/full_data_ss.rds")
indexing_tab <- readRDS("data/input/mdl_data/indexing_tab_data_ss.rds")
post_ind <- indexing_tab
site_enviro_metrics <- mdl_data %>% 
  select(sp_site,sd_temp_60,sd_total_precip60) %>% 
  unique()
post_ind$species <- unique(mdl_data$species)[post_ind$sp]

grid_elev_bands <- st_read("data/input/spatial/grid_elev_bands.shp")

# Get number of datapoint in phenology dataset
n_by_grid <- mdl_data %>% 
  group_by(cell,elev) %>%
  summarise(count = n())

# Merge with elev, cell
grid_elev_bands$seqnum <- 1:nrow(grid_elev_bands)
grid_elev_bands <- merge(grid_elev_bands,n_by_grid,by=c("elev","cell"))

#hexgrid6 <- dggridR::dgconstruct(res = 6) 
#grid  <- dgcellstogrid(hexgrid6,grid_elev_bands$cell,return_sf = FALSE)

# Convert multipolygon to points, keeping attributes
grid_elev_bands_cols <- st_drop_geometry(grid_elev_bands)
grid_elev_bands_cols$L3 <- 1:nrow(grid_elev_bands_cols)

grid  <- grid_elev_bands %>%
  st_coordinates()
grid <- merge(grid,grid_elev_bands_cols,by="L3")

colnames(grid) <- c("group","long","lat","L1","L2","elev","cell","seqnum","count")

#grid  <- merge(grid,samp_size,by=c("cell","elev"))
#grid <- st_drop_geometry(grid,)

library(rnaturalearth)
library(rnaturalearthdata)

water_col = "grey40"
country_col = "grey20"
midpoint_col = "grey90"
country_outline_col = "grey50"
water_col = "grey40"

proj_info <- "+proj=natearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"

# Get world countries data
world <- ne_countries(scale = "medium", returnclass = "sf")
world_rp <- st_transform(world, crs = proj_info)
grid_elev_bands_rp <- st_transform(grid_elev_bands, 
                                           crs = proj_info)


grid_elev_bands_rp$log10_count <- log10(grid_elev_bands_rp$count)

#bird_cols <- bird_colors("Cassin's Finch",expand_palette = 10)[5:10]
bird_cols <- bird_colors("Variable Sunbird",reverse = F,expand_palette = 100)[1:80]

# Plot with geom_sf (simplest approach)
ggplot(grid_elev_bands_rp) +
  geom_sf(data = world_rp, fill = country_col, color = country_outline_col, size = 0.1) +
  geom_sf(aes(fill = count), color = "transparent", size = 0.2) +
  scale_fill_bird(bird_cols = bird_cols) +
  scale_y_log10()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        legend.key.size = unit(1, 'cm'),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position = c(0.2,0.35),
        legend.background = element_rect(fill=alpha(water_col))) +
  theme(panel.background = element_rect(fill=alpha(water_col)))

### Projection figure ####

#Simple conceptual figure of projections into the future of departure

plot(NULL,xlim=c(2002,2100),ylim=c(-5,5),ylab="Departure Date",xlab="Year",
     yaxt="n",cex.lab=1.6,cex.axis=1.6)
#abline(h=0,lty=2)

dep_trend <- c(-0.025,0,0.025)
ints <- c(-1,0,1)
cols <- bird_colors("Lazuli Bunting",rev=TRUE)
#Add 10 lines
for (nn in 1:3){
  dep_trend_in <- dep_trend[nn]
  sim_known_trend <- ints[nn] + dep_trend_in*1:23 + rnorm(23,0,0.25)
  points(2002:2024,sim_known_trend,type="l",col = cols[nn],lwd=6)
  
  sim_future_range_rec <- c()
  for (mm in 1:100){
    proj_unc <- 1:76 * 0.05
    sim_future_range <- ints[nn] + dep_trend_in*24:99 + rnorm(76,0,0.25+proj_unc)
    sim_future_range_rec <- rbind(sim_future_range_rec,sim_future_range)
  }
  mean_proj <- c(sim_known_trend[23],apply(sim_future_range_rec, 2, mean))
  upper_proj <- c(sim_known_trend[23],apply(sim_future_range_rec, 2, quantile,.975))
  lower_proj <- c(sim_known_trend[23],apply(sim_future_range_rec, 2, quantile,.025))
  
  points(2024:2100,mean_proj,type="l",lty=3,col=alpha(cols[nn],.95),lwd=3)
  #polygon(c(2024:2100,rev(2024:2100)),c(lower_proj,rev(upper_proj)),
  #        col=alpha(cols[nn],.5),border=NA)
}




