library(dplyr)
library(birdcolors)
library(ggplot2)
library(rstanarm)
library(MCMCvis)

# Only need to run this section once! This will save data for coloring the plots.

# mdl_fit <- readRDS("data/output/departure_model/departure_climate_fit_ss2.rds")
# 
# mu_beta1s <- MCMCsummary(mdl_fit,params="mu_beta1")
# mu_beta2s <- MCMCsummary(mdl_fit,params="mu_beta2")
# mdl_fit <- NULL
# gc()

#saveRDS(mu_beta1s,"data/output/variance_explained/mu_beta1s.rds")
#saveRDS(mu_beta2s,"data/output/variance_explained/mu_beta2s.rds")

mu_beta1s <- readRDS("data/output/variance_explained/mu_beta1s.rds")
var_exp_env <- readRDS("data/output/variance_explained/var_exp_env.rds")
var_total <- readRDS("data/output/variance_explained/total_var.rds")

best_corr <- readRDS("data/output/variance_explained/site_data_best_corr.rds")

sp_summary_temp <- readRDS("data/output/variance_explained/sp_summary_temp.rds")
colnames(sp_summary_temp)[2:5] <- c("mean_max_r_sq_t","mean_95_r_sq_t","mean_mean_r_sq_t","count_areas_t")

sp_summary_precip <- readRDS("data/output/variance_explained/sp_summary_precip.rds")
colnames(sp_summary_precip)[2:5] <- c("mean_max_r_sq_p","mean_95_r_sq_p","mean_mean_r_sq_p","count_areas_p")

var_exp_comb <- merge(var_exp_env,var_total,by="species")
var_exp_comb$mu_beta1s <- mu_beta1s$mean
var_exp_comb <- merge(var_exp_comb,sp_summary_temp,by="species")
var_exp_comb <- merge(var_exp_comb,sp_summary_precip,by="species")

var_exp_comb_analysis <- var_exp_comb[complete.cases(var_exp_comb),]
var_exp_comb_analysis$best_case_max <- best_corr$mean_best_R2
#Remove species with recent lumps/splits that make range maps outdated
sp_remove <- c("Siberian Pipit","Eastern Cattle-Egret","Siberian Stonechat",
               "Amur Stonechat","European Stonechat")

var_exp_comb_analysis <- var_exp_comb_analysis %>%
  filter(!species %in% sp_remove)

### Combined (best case scenario temp and precip) ####

var_env_lm <- stan_lm(formula = mean_env~best_case_max, data=var_exp_comb_analysis,prior = NULL)
MCMCsummary(var_env_lm,pg0 = TRUE,probs = c(0.055,.5,0.945))
var_env_lm_chains <- MCMCchains(var_env_lm)

#saveRDS(var_env_lm,"data/output/mig_dist/var_env_lm.rds")

xx <- seq(min(var_exp_comb_analysis$mean_max_r_sq_t,na.rm=TRUE),max(var_exp_comb_analysis$mean_max_r_sq_t,na.rm=TRUE),by=.01)

e_rec <- matrix(NA,ncol=length(xx),nrow = 4000)
for(nn in 1:4000){
  e_rec[nn,] <- var_env_lm_chains[nn,1] + xx * var_env_lm_chains[nn,2]
}
exp_mean <- colMeans(e_rec)
exp_lower89 <- apply(e_rec,2,quantile,.055)
exp_upper89 <- apply(e_rec,2,quantile,.945)


ggplot()+
  geom_point(data=var_exp_comb_analysis,aes(best_case_max,mean_env,col=mu_beta1s)) +
  scale_color_bird(bird_colors("Bluethroat",rev=F)[c(2:4)],midpoint = 0) +
  geom_ribbon(aes(x = xx,ymin=exp_lower89,ymax=exp_upper89),fill="grey",alpha=0.8) +
  geom_line(aes(x = xx,y=exp_mean),col= "black",size=2,alpha=0.9) +  
  xlab("Climatic Connectivity") +
  ylab("Variance Explained") +
  labs(col="Temperature Sensitivity (days/SD)") +
  theme_classic(base_size = 16)

MCMCsummary(var_env_lm,params="R2",probs = c(0.055,.5,0.945))
#saveRDS(var_env_lm,"data/output/clim_con/var_env_lm.rds") # To save, uncomment this line


var_total_lm <- stan_lm(formula = mean_var~best_case_max, data=var_exp_comb_analysis,prior = NULL)
MCMCsummary(var_total_lm,pg0 = TRUE,probs = c(0.055,.5,0.945))
var_total_lm_chains <- MCMCchains(var_total_lm)

xx <- seq(min(var_exp_comb_analysis$mean_max_r_sq_p,na.rm=TRUE),
          max(var_exp_comb_analysis$mean_max_r_sq_p,na.rm=TRUE),by=.01)

e_rec <- matrix(NA,ncol=length(xx),nrow = 4000)
for(nn in 1:4000){
  e_rec[nn,] <- var_total_lm_chains[nn,1] + xx * var_total_lm_chains[nn,2]
}
exp_mean <- colMeans(e_rec)
exp_lower89 <- apply(e_rec,2,quantile,.055)
exp_upper89 <- apply(e_rec,2,quantile,.945)

ggplot()+
  geom_point(data=var_exp_comb_analysis,aes(best_case_max,mean_var,col=mu_beta1s)) +
  scale_color_bird(bird_colors("Bluethroat",rev=F)[c(2:4)],midpoint = 0) +
  geom_ribbon(aes(x = xx,ymin=exp_lower89,ymax=exp_upper89),fill="grey",alpha=0.8) +
  geom_line(aes(x = xx,y=exp_mean),col= "black",size=2,alpha=0.9) +  
  labs(col="Temperature Sensitivity (days/SD)") +
  xlab("Climatic Connectivity") +
  ylab("Variance") +
  theme_classic(base_size = 16)

MCMCsummary(var_total_lm,params="R2",probs = c(0.055,.5,0.945))
#saveRDS(var_total_lm,"data/output/clim_con/var_total_lm.rds") # To save, uncomment this line


### To read in saved model data ####
var_env_lm <- readRDS("data/output/clim_con/var_env_lm.rds")
MCMCsummary(var_env_lm,pg0 = T,probs = c(0.055,.5,0.945))

var_total_lm <- readRDS("data/output/clim_con/var_total_lm.rds")
MCMCsummary(var_total_lm,pg0 = T,probs = c(0.055,.5,0.945))
