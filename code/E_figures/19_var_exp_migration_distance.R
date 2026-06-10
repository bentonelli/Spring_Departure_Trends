library(dplyr)
library(birdcolors)
library(ggplot2)
library(rstanarm)
library(MCMCvis)

mdl_fit <- readRDS("data/output/departure_model/departure_climate_fit_ss2.rds")

mu_beta1s <- MCMCsummary(mdl_fit,params="mu_beta1")

mdl_fit <- NULL
gc()

mig_dist <- readRDS("data/output/mig_dist/sp_mig_data_all.rds")
var_exp_env <- readRDS("data/output/variance_explained/var_exp_env.rds")
var_total <- readRDS("data/output/variance_explained/total_var.rds")

var_exp_comb <- merge(var_exp_env,mig_dist,by="species")
var_exp_comb <- merge(var_exp_comb,var_total,by="species")

var_exp_comb_analysis <- var_exp_comb[complete.cases(var_exp_comb),]

var_env_lm <- stan_lm(formula = mean_env~migration_dist_km, data=var_exp_comb_analysis,prior = NULL)
MCMCsummary(var_env_lm,pg0 = TRUE,probs = c(0.055,.5,0.945))
var_env_lm_chains <- MCMCchains(var_env_lm)

saveRDS(var_env_lm,"data/output/mig_dist/var_env_lm.rds")

xx <- seq(min(var_exp_comb_analysis$migration_dist_km,na.rm=TRUE),max(var_exp_comb_analysis$migration_dist_km,na.rm=TRUE),by=10)

e_rec <- matrix(NA,ncol=length(xx),nrow = 4000)
for(nn in 1:4000){
  e_rec[nn,] <- var_env_lm_chains[nn,1] + xx * var_env_lm_chains[nn,2]
}
exp_mean <- colMeans(e_rec)
exp_lower89 <- apply(e_rec,2,quantile,.055)
exp_upper89 <- apply(e_rec,2,quantile,.945)


ggplot()+
  geom_point(data=var_exp_comb,aes(migration_dist_km,mean_env,col=mu_beta1s$mean)) +
  scale_color_bird(bird_colors("Bluethroat",rev=TRUE)[c(2:4)],midpoint = 0) +
  geom_ribbon(aes(x = xx,ymin=exp_lower89,ymax=exp_upper89),fill="grey",alpha=0.8) +
  geom_line(aes(x = xx,y=exp_mean),col= "black",size=2,alpha=0.9) +  
  xlab("Migration Distance (km)") +
  ylab("Variance Explained") +
  labs(col="Temperature Sensitivity (days/SD)") +
  theme_classic(base_size = 16)

MCMCsummary(var_env_lm,params="R2",probs = c(0.055,.5,0.945))

### Total Variance ####
var_total_lm <- stan_lm(formula = mean_var~migration_dist_km, data=var_exp_comb_analysis,prior = NULL)
MCMCsummary(var_total_lm,pg0 = TRUE,probs = c(0.055,.5,0.945))
var_total_lm_chains <- MCMCchains(var_total_lm)

e_rec <- matrix(NA,ncol=length(xx),nrow = 4000)
for(nn in 1:4000){
  e_rec[nn,] <- var_total_lm_chains[nn,1] + xx * var_total_lm_chains[nn,2]
}
exp_mean <- colMeans(e_rec)
exp_lower89 <- apply(e_rec,2,quantile,.055)
exp_upper89 <- apply(e_rec,2,quantile,.945)

ggplot()+
  geom_point(data=var_exp_comb,aes(migration_dist_km,mean_var,col=mu_beta1s$mean)) +
  scale_color_bird(bird_colors("Bluethroat",rev=TRUE)[c(2:4)],midpoint = 0) +
  geom_ribbon(aes(x = xx,ymin=exp_lower89,ymax=exp_upper89),fill="grey",alpha=0.8) +
  geom_line(aes(x = xx,y=exp_mean),col= "black",size=2,alpha=0.9) +  
  labs(col="Temperature Sensitivity (days/SD)") +
  xlab("Migration Distance (km)") +
  ylab("Variance") +
  theme_classic(base_size = 16)

MCMCsummary(var_total_lm,params="R2",probs = c(0.055,.5,0.945))
saveRDS(var_total_lm,"data/output/mig_dist/var_total_lm.rds")
