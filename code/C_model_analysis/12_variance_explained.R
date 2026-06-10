# Code to look at variance explained at the species- site- and species-site-level

library(MCMCvis)
library(dplyr)
library(ggplot2)

#Read in model
mdl_in <- readRDS("data/output/departure_model/departure_climate_fit_ss2.rds")
ind_data <- readRDS("data/input/mdl_data/indexing_tab_data_ss.rds")
full_data <- readRDS("data/input/mdl_data/full_data_ss.rds")
list_data <- readRDS("data/input/mdl_data/list_format_data_ss.rds")

sp_ind <- full_data$sp
sp_site_ind <- full_data$sp_site
site_ind <- full_data$site_ind

temp_anom <- full_data$temp60_anom
precip_anom <- full_data$precip60_anom

yy_trues <- MCMCchains(mdl_in,params="yy_true")
alphas <- MCMCchains(mdl_in,params="alpha")
beta1s <- MCMCchains(mdl_in,params="beta1")
beta2s <- MCMCchains(mdl_in,params="beta2")

mu_beta1s <- MCMCchains(mdl_in,params="mu_beta1")
mu_beta2s <- MCMCchains(mdl_in,params="mu_beta2")

mdl_in <- NULL #Clear memory

var_exp_temp_rec <- matrix(NA,nrow=max(sp_site_ind),ncol=20000)
var_exp_precip_rec <- matrix(NA,nrow=max(sp_site_ind),ncol=20000)
var_exp_env_rec <- matrix(NA,nrow=max(sp_site_ind),ncol=20000)
total_variance <- matrix(NA,nrow=max(sp_site_ind),ncol=20000)

for(sp_site_target in 1:max(sp_site_ind)){
  print(sp_site_target)
  sp_site_target_ind <- which(sp_site_ind == sp_site_target)
  
  anom <- t(yy_trues[,sp_site_target_ind]) - alphas[,sp_site_target]
  
  total_var <- apply(anom, 2, var)
  
  b1_mat <- matrix(beta1s[,sp_site_target], nrow = length(sp_site_target_ind), ncol = length(beta1s[,sp_site_target]), byrow = TRUE)
  b2_mat <- matrix(beta2s[,sp_site_target], nrow = length(sp_site_target_ind), ncol = length(beta2s[,sp_site_target]), byrow = TRUE)
  
  grep_temp = b1_mat * temp_anom[sp_site_target_ind] + b2_mat * precip_anom[sp_site_target_ind]
  grep_precip = b2_mat * precip_anom[sp_site_target_ind]
  grep_env = grep_temp + grep_precip
  
  e_temp = anom - grep_temp
  e_precip = anom - grep_precip
  e_env = anom - grep_env
  
  grep_temp_sp_cl <- apply(grep_temp, 2, var)
  e_temp_sp_cl <- apply(e_temp, 2, var)
  
  grep_precip_sp_cl <- apply(grep_precip, 2, var)
  e_precip_sp_cl <- apply(e_precip, 2, var)
  
  grep_env_sp_cl <- apply(grep_env, 2, var)
  e_env_sp_cl <- apply(e_env, 2, var)
  
  var_exp_temp_rec[sp_site_target,] <- (grep_temp_sp_cl)/(grep_temp_sp_cl+e_temp_sp_cl)
  var_exp_precip_rec[sp_site_target,] <- (grep_precip_sp_cl)/(grep_precip_sp_cl+e_precip_sp_cl)
  var_exp_env_rec[sp_site_target,] <- (grep_env_sp_cl)/(grep_env_sp_cl+e_env_sp_cl)
  
  total_variance[sp_site_target,] <- total_var
}

### Now calculate by species ####
sp_level_est_temp <- c()
sp_level_est_precip <- c()
sp_level_est_env <- c()
sp_level_var_exp <- c()

for (each_sp in 1:max(ind_data$sp)){
  
  sp_ind <- which(ind_data$sp == each_sp)
  
  iter_est_temp <- colMeans(var_exp_temp_rec[sp_ind,])
  iter_est_precip <- colMeans(var_exp_precip_rec[sp_ind,])
  iter_est_env <- colMeans(var_exp_env_rec[sp_ind,])
  iter_est_var <- colMeans(total_variance[sp_ind,])
  
  sp_level_est_temp <- rbind(sp_level_est_temp,c(mean(iter_est_temp),quantile(iter_est_temp,.055),quantile(iter_est_temp,.945)))
  sp_level_est_precip <- rbind(sp_level_est_precip,c(mean(iter_est_precip),quantile(iter_est_precip,.055),quantile(iter_est_precip,.945)))
  sp_level_est_env <- rbind(sp_level_est_env,c(mean(iter_est_env),quantile(iter_est_env,.055),quantile(iter_est_env,.945)))
  
  sp_level_var_exp <- rbind(sp_level_var_exp,c(mean(iter_est_var),quantile(iter_est_var,.055),quantile(iter_est_var,.945)))
}

sp_level_est_temp <- as.data.frame(sp_level_est_temp)
sp_level_est_temp$species <- unique(full_data$species)
colnames(sp_level_est_temp)[1:3] <- c("mean_temp","lower_5p5_temp","upper_94p5_temp")

sp_level_est_precip <- as.data.frame(sp_level_est_precip)
sp_level_est_precip$species <- unique(full_data$species)
colnames(sp_level_est_precip)[1:3] <- c("mean_precip","lower_5p5_precip","upper_94p5_precip")

sp_level_est_env <- as.data.frame(sp_level_est_env)
sp_level_est_env$species <- unique(full_data$species)
colnames(sp_level_est_env)[1:3] <- c("mean_env","lower_5p5_env","upper_94p5_env")

sp_level_var_exp <- as.data.frame(sp_level_var_exp)
sp_level_var_exp$species <- unique(full_data$species)
colnames(sp_level_var_exp)[1:3] <- c("mean_var","lower_5p5_var","upper_94p5_var")

saveRDS(sp_level_est_temp,"data/output/variance_explained/var_exp_temp.rds")
saveRDS(sp_level_est_precip,"data/output/variance_explained/var_exp_precip.rds")
saveRDS(sp_level_est_env,"data/output/variance_explained/var_exp_env.rds")

saveRDS(sp_level_var_exp,"data/output/variance_explained/total_var.rds")
