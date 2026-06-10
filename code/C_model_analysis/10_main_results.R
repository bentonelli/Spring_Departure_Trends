library(MCMCvis)
library(ggplot2)
library(shinystan)
library(dplyr)
library(birdcolors)
library(shredder)

#Model fit and data
mdl_fit <- readRDS("data/output/departure_model/departure_climate_fit_ss2.rds")

#Save diagnostics
#MCMCdiag(mdl_fit,file_name = "data/output/departure_model/mdl_diagnostics.txt")

mdl_fit_filt <- rstan::extract(mdl_fit,pars=c("mu_nu1","mu_nu2","zeta1","zeta2","mu_eta","xi"),permuted=TRUE)
mat <- as.matrix(as.data.frame(mdl_fit_filt))
MCMCsummary(mat,pg0 = TRUE,probs = c(0.055,.5,0.945))

mas <- rstan::extract(mdl_fit,pars=c("mu_alpha"),permuted=TRUE)
mas_df <- as.matrix(as.data.frame(mas))
mas_sum <- MCMCsummary(mas_df,pg0 = TRUE,probs = c(0.055,.5,0.945))

mb1s <- rstan::extract(mdl_fit,pars=c("mu_beta1"),permuted=TRUE)
mb1s_df <- as.matrix(as.data.frame(mb1s))
mb1s_sum <- MCMCsummary(mb1s_df,pg0 = TRUE,probs = c(0.055,.5,0.945))

mb2s <- rstan::extract(mdl_fit,pars=c("mu_beta2"),permuted=TRUE)
mb2s_df <- as.matrix(as.data.frame(mb2s))
mb2s_sum <- MCMCsummary(mb2s_df,pg0 = TRUE,probs = c(0.055,.5,0.945))

mdl_data <- readRDS("data/input/mdl_data/full_data_ss.rds")
indexing_tab <- readRDS("data/input/mdl_data/indexing_tab_data_ss.rds")

mas_sum$sp <- unique(mdl_data$species)

mb1s_sum$sp <- unique(mdl_data$species)

mb2s_sum$sp <- unique(mdl_data$species)

# Relative shifts
sum(mb1s_sum$mean < -1 | mb1s_sum$mean > 1)/nrow(mb1s_sum)
sum(mb2s_sum$mean < -1 | mb2s_sum$mean > 1)/nrow(mb2s_sum)

mdl_fit <- NULL
gc()


  

