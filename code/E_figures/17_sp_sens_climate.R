library(dplyr)
library(MCMCvis)
library(ggplot2)
library(birdcolors)
library(ggbeeswarm)
library(clootl)
library(phytools)
library(readr)

mdl_fit <- readRDS("data/output/departure_model/departure_climate_fit_ss2.rds")
mdl_data <- readRDS("data/input/mdl_data/full_data_ss.rds")
indexing_tab <- readRDS("data/input/mdl_data/indexing_tab_data_ss.rds")
list_data <- readRDS("data/input/mdl_data/list_format_data_ss.rds")

sp_mean_dep <- MCMCsummary(mdl_fit,"mu_alpha")
mu_beta1s <- MCMCsummary(mdl_fit,"mu_beta1")
mu_beta2s <- MCMCsummary(mdl_fit,"mu_beta2")
theta1_out <- MCMCsummary(mdl_fit,"theta1")
zeta1_out <- MCMCsummary(mdl_fit,"zeta1")
theta1_chains <- MCMCchains(mdl_fit,params = "theta1")
zeta1_chains <- MCMCchains(mdl_fit,params = "zeta1")
theta2_chains <- MCMCchains(mdl_fit,params = "theta2")
zeta2_chains <- MCMCchains(mdl_fit,params = "zeta2")

mdl_fit <- NULL

sp_info <- indexing_tab %>% 
  select(sp,mean_lon,mean_lat,mean_temp_niche,mean_precip_niche) %>%
  unique()

avg_untransformed_temp <- mean(sp_info$mean_temp_niche)
avg_untransformed_precip <- mean(sp_info$mean_precip_niche)
avg_untransformed_temp_C <- avg_untransformed_temp - 273.15

sp_info$mu_alpha <- sp_mean_dep$mean
sp_info$mu_beta1s <- mu_beta1s$mean
sp_info$mu_beta2s <- mu_beta2s$mean

sp_info$mean_temp_niche <- list_data$MT
sp_info$mean_precip_niche <- list_data$MP
sp_info$mean_precip_niche_mm_p_day <- 1000*(list_data$MP/1000 + avg_untransformed_precip)

sp_info$pred_mu_beta1 <- theta1_out$mean + zeta1_out$mean * sp_info$mean_temp_niche
sp_info$resid_mu_beta1 <- sp_info$mu_beta1s - sp_info$pred_mu_beta1

sp_info$species <- unique(mdl_data$species)
#
col_div <- bird_colors("Northern Parula",reverse = TRUE)

point_size <- 2.5
alpha_points <- .75

line_width <- 2
line_alpha <- 0.75
line_col <- "grey30"

ribbon_col <- "grey60"


### Temperature Plot ####

xx <- min(sp_info$mean_temp_niche):max(sp_info$mean_temp_niche)

e_rec <- matrix(NA,ncol=length(xx),nrow = 5000)
for(nn in 1:5000){
  e_rec[nn,] <- theta1_chains[nn] + xx * zeta1_chains[nn]
}
exp_mean <- colMeans(e_rec)
exp_lower89 <- apply(e_rec,2,quantile,.055)
exp_upper89 <- apply(e_rec,2,quantile,.945)

ggplot() +
  geom_point(sp_info,mapping=aes(mean_temp_niche+avg_untransformed_temp_C, mu_beta1s,
                                 col=abs(mean_lat)),alpha=alpha_points,size=point_size) +
  geom_ribbon(aes(x = xx+avg_untransformed_temp_C,ymin=exp_lower89,ymax=exp_upper89),
              fill=ribbon_col,alpha=line_alpha) +
  
  labs(col="") +
  geom_line(aes(x = xx+avg_untransformed_temp_C,y=exp_mean),col= line_col,size=line_width,alpha=line_alpha) +
  geom_hline(yintercept = 0,col=line_col,lty=2) +
  scale_color_bird(col_div) + 
  ylim(-4.5,2) +
  xlim(-12,32) +
  xlab("Average Temperature (C)") +
  ylab("Temperature Sensitivity (days/+1SD)") +
  theme_classic() +
  theme_classic(base_size = 20)+
  theme(legend.key.size = unit(0.75, 'cm'),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.position = c(0.85,0.25),
        legend.background = element_rect(fill=alpha("transparent"))) +
  geom_vline(xintercept = sp_info$mean_temp_niche[182]+avg_untransformed_temp_C) +
  geom_hline(yintercept = sp_info$mu_beta1s[182])

# Snow Bunting = 304
# Dunlin = 107
# Indian Paradise-Flycatcher = 182

### Precipitation plot ####

xx <- seq(min(sp_info$mean_precip_niche),max(sp_info$mean_precip_niche),by=.1)

e_rec <- matrix(NA,ncol=length(xx),nrow = 5000)
for(nn in 1:5000){
  e_rec[nn,] <- theta2_chains[nn] + xx * zeta2_chains[nn]
}
exp_mean <- colMeans(e_rec)
exp_lower89 <- apply(e_rec,2,quantile,.055)
exp_upper89 <- apply(e_rec,2,quantile,.945)

xx_trans <- 1000*(xx/1000 + avg_untransformed_precip)

ggplot() +
  geom_point(sp_info,mapping=aes(mean_precip_niche_mm_p_day,mu_beta2s,
             col=abs(mean_lat)),alpha=alpha_points,size=point_size) +
  geom_ribbon(aes(x = xx_trans,ymin=exp_lower89,ymax=exp_upper89),fill=ribbon_col,alpha=line_alpha) +
  labs(col="") +
  geom_line(aes(x = xx_trans,y=exp_mean),col= line_col,size=line_width,alpha=line_alpha) +
  geom_hline(yintercept = 0,col=line_col,lty=2) +
  scale_color_bird(col_div) + 
  xlab("Average Precipitation (mm/day)") +
  ylab("Precipitation Sensitivity (days/+1SD)") +
  theme_classic(base_size = 20)+
  theme(legend.key.size = unit(0.75, 'cm'),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.position = c(0.85,0.25),
        legend.background = element_rect(fill=alpha("transparent")))
                
