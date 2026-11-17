# Script to read in departure predictions, create ensemble forecasts

library(dplyr)
library(ggplot2)
library(birdcolors)
library(zoo)
library(ggbeeswarm)
library(tidyr)


# Read in model data
mdl_data <- readRDS("data/input/mdl_data/full_data_ss.rds")
indexing_tab <- readRDS("data/input/mdl_data/indexing_tab_data_ss.rds")

sp_it <- mdl_data %>% 
  select(c("species","mean_dep","elev",
           "lat","lon",
           "species","mean_dep",
           "mean_temp60","sd_temp_60",
           "mean_total_precip60","sd_total_precip60")) %>%
  unique()

# Read in modeled projections
mdl1_proj <- readRDS("data/output/depature_projections/est_dep_CESM2_SSP1_2p6.rds")
mdl2_proj <- readRDS("data/output/depature_projections/est_dep_CNRM-ESM2-1_SSP1_2p6.rds")
mdl3_proj <- readRDS("data/output/depature_projections/est_dep_MIROC6_SSP1_2p6.rds")

ensemble_avg_proj_2p6 <- as.matrix((mdl1_proj+mdl2_proj+mdl3_proj)/3)

mdl1_proj <- readRDS("data/output/depature_projections/est_dep_CESM2_SSP2_4p5.rds")
mdl2_proj <- readRDS("data/output/depature_projections/est_dep_CNRM-ESM2-1_SSP2_4p5.rds")
mdl3_proj <- readRDS("data/output/depature_projections/est_dep_MIROC6_SSP2_4p5.rds")

ensemble_avg_proj_4p5 <- as.matrix((mdl1_proj+mdl2_proj+mdl3_proj)/3)

mdl1_proj <- readRDS("data/output/depature_projections/est_dep_CESM2_SSP3_7p0.rds")
mdl2_proj <- readRDS("data/output/depature_projections/est_dep_CNRM-ESM2-1_SSP3_7p0.rds")
mdl3_proj <- readRDS("data/output/depature_projections/est_dep_MIROC6_SSP3_7p0.rds")

ensemble_avg_proj_7p0 <- as.matrix((mdl1_proj+mdl2_proj+mdl3_proj)/3)

ens_list <- list(ensemble_avg_proj_2p6,ensemble_avg_proj_4p5,ensemble_avg_proj_7p0)

### Global options ###
xlim_all <- c(2010,2105)
export_size_x <- 14
export_size_y <- 6
base_size_all <- 26


### Plot 1, a few species ####

bc_in <- bird_colors(palette_name = "Northern Parula",expand_palette = 536,rev=F)
sp_plot_col <- c(bc_in[137],bc_in[413],bc_in[479])
sp_list <- c("Snow Bunting","Dunlin","Indian Paradise-Flycatcher")

sp_targ_long_all <- c()
new_max_pop <- 0
for (ll in 1:3){
  sp_it_filt <- which(sp_it$species == sp_list[ll])
  sp_targ_change <- ens_list[[2]][sp_it_filt,]
  sp_targ_change_df <- as.data.frame(sp_targ_change)
  colnames(sp_targ_change_df) <- 1:ncol(sp_targ_change_df)
  sp_targ_change_df$scale = "Population"
  sp_targ_change_df$pop_num = 1:nrow(sp_targ_change_df)+new_max_pop
  sp_targ_change_df$lat <- mean(sp_it$lat[sp_it_filt])
  
  new_max_pop <- max(sp_targ_change_df$pop_num)
  
  sp_targ_change_df <- rbind(sp_targ_change_df,c(colMeans(sp_targ_change,na.rm=TRUE),
                                                 scale="Species",1000+ll,lat=mean(sp_it$lat[sp_it_filt])))
  
  sp_targ_long <- pivot_longer(sp_targ_change_df,cols=1:86,
                               names_to = "year",values_to = "est_dep")
  sp_targ_long$species <- sp_list[ll]
  sp_targ_long$pop_num <- as.numeric(sp_targ_long$pop_num)
  sp_targ_long$year <- as.numeric(sp_targ_long$year)
  sp_targ_long$est_dep <- as.numeric(sp_targ_long$est_dep)
  
  sp_targ_long_all <- rbind(sp_targ_long_all,sp_targ_long)
}

head(sp_targ_long_all)

pp1a <- ggplot(sp_targ_long_all) +
  geom_line(aes(x=2015+year,y=est_dep,group=pop_num,
                col=factor(species,levels=c("Snow Bunting","Dunlin","Indian Paradise-Flycatcher")),
                alpha = factor(scale),
                linewidth = factor(scale)))+
  ylab("Departure Date Change (days)") +
  xlab("Year") +
  scale_x_continuous(
    limits = xlim_all,
    breaks = seq(2015, 2100, by = 5),
    expand=c(0,0)) +
  scale_linewidth_manual(values=c(.5,2.25),guide="none") +
  scale_alpha_manual(values=c(.25,.9),guide="none") + 
  scale_color_manual(values = sp_plot_col)+
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  labs(col="") +
  theme_bw(base_size = base_size_all) +
  theme(legend.key.size = unit(1, 'cm'),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.position = c(0.175,0.25),
        legend.background = element_rect(fill=alpha("transparent")),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  coord_cartesian(xlim = xlim_all, ylim = c(-11, 7))
ggsave("plot1a.pdf", plot = pp1a, width = export_size_x, height = export_size_y, units = "in")   

pp1b <- ggplot(sp_targ_long_all) +
  geom_line(aes(x=2015+year,y=est_dep,group=pop_num,
                col=factor(species,levels=c("Snow Bunting","Dunlin","Indian Paradise-Flycatcher")),
                alpha = factor(scale),
                linewidth = factor(scale)))+
  ylab("Departure Date Change (days)") +
  xlab("Year") +
  scale_x_continuous(
    limits = xlim_all,
    breaks = seq(xlim_all[1],xlim_all[2],by=5),
    labels = function(x) paste0(x-4, "-", x),
    expand=c(0,0)
  ) +
  scale_linewidth_manual(values=c(.5,2.25),guide="none") +
  scale_alpha_manual(values=c(.25,.9),guide="none") + 
  scale_color_manual(values = sp_plot_col)+
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  labs(col="") +
  theme_bw(base_size = base_size_all) +
  theme(legend.key.size = unit(1, 'cm'),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.position = c(0.175,0.25),
        legend.background = element_rect(fill=alpha("transparent")),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  coord_cartesian(xlim = xlim_all, ylim = c(-11, 7))
ggsave("plot1b.pdf", plot = pp1b, width = export_size_x, height = export_size_y, units = "in")   

### Plot 2 - Ensemble forecast across species ####
ensemble_avg_proj_sp <- ens_list[[2]] %>% 
  as.data.frame() %>% 
  mutate(sp = sp_it$species)

sp_dep_proj <- c()
for(each_sp in unique(ensemble_avg_proj_sp$sp)){
  sp_dt <- filter(ensemble_avg_proj_sp,sp == each_sp)
  sp_add <- c(colMeans(sp_dt[1:(ncol(sp_dt)-1)]))
  sp_dep_proj <- rbind(sp_dep_proj,sp_add)
}

colnames(sp_dep_proj) <- 2015:2100
rownames(sp_dep_proj) <- unique(ensemble_avg_proj_sp$sp)

rollmean_dep_proj <- apply(sp_dep_proj,1,zoo::rollmean,k=5,fill=NA,na.pad=TRUE,align = "right")
rownames(rollmean_dep_proj) <- 2015:2100

rollmean_dep_proj_filt <- rollmean_dep_proj %>%
  as.data.frame() %>%
  filter(rownames(rollmean_dep_proj) %in% seq(2025,2100,by=5))

rollmean_dep_proj_filt_rn <- t(rollmean_dep_proj_filt) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "species")

long_dep <- pivot_longer(rollmean_dep_proj_filt_rn,2:ncol(rollmean_dep_proj_filt_rn),
                         values_to = "dep",names_to = "year")

sp_characteristics <- indexing_tab %>% 
  select(sp,mean_lat,mean_temp_niche,mean_precip_niche) %>%
  unique()

sp_characteristics$species <- mdl_data$species %>% unique()
sp_characteristics$mean_temp_niche <- sp_characteristics$mean_temp_niche - 273.15

long_dep <- merge(long_dep,sp_characteristics,by="species")

long_dep$year <- as.numeric(long_dep$year)

# For species level summary:
eoc_change <- long_dep %>% filter(year %in% c(2095,2100)) %>% group_by(species) %>% summarise(mean_dep_eoc = mean(dep))

bird_col <- bird_colors("Northern Parula",rev=TRUE)

quantile_segments <- long_dep %>%
  group_by(year) %>%
  summarise(q10 = quantile(dep,.1),
            q90 = quantile(dep,.9))

bcs_quantiles <- bird_colors(palette_name = "Cassin's Finch",expand_palette = 4)[2:4]#c("grey10","grey30","grey40")#bird_colors("Blue Pitta",rev=FALSE)[2:4]

pp2 <- ggplot()+
  geom_hline(yintercept = 0,linetype="dashed",col="darkgrey")+
  geom_quasirandom(data=long_dep,aes(x=(year),y = dep,
                                     color=abs(mean_lat)),
                   size=0.1,alpha=0.75,pch=19) +
  scale_color_bird(bird_col) +
  labs(col="Latitude (abs)") +
  geom_segment(data=quantile_segments,aes(x=year-2,xend=year+2,y=q10),lwd=1.25,col=alpha(bcs_quantiles[2],.9))+
  geom_segment(data=quantile_segments,aes(x=year-2,xend=year+2,y=q90),lwd=1.25,col=alpha(bcs_quantiles[2],.9))+
  ylab("Departure Date Change (days)")+
  xlab("Year") +
  scale_x_continuous(
    limits = xlim_all,
    breaks = seq(xlim_all[1],xlim_all[2],by=5),
    labels = function(x) paste0(x-4, "-", x),
    expand=c(0,0)
  ) +
  theme_bw(base_size = base_size_all) +
  theme(legend.key.size = unit(0.75, 'cm'),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position = c(0.075,0.275),
        legend.background = element_rect(fill=alpha("transparent")),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  coord_cartesian(xlim = xlim_all, ylim = c(-11, 7))
ggsave("plot2.pdf", plot = pp2, width = export_size_x, height = export_size_y, units = "in")   

### Plot 3, SD change across three SSP scenarios ####
quantile_segments <- list()
for(ssp in 1:length(ens_list)){
  ensemble_avg_proj_sp <- ens_list[[ssp]] %>% 
    as.data.frame() %>% 
    mutate(sp = sp_it$species)
  
  sp_dep_proj <- c()
  for(each_sp in unique(ensemble_avg_proj_sp$sp)){
    sp_dt <- filter(ensemble_avg_proj_sp,sp == each_sp)
    sp_add <- c(colMeans(sp_dt[1:(ncol(sp_dt)-1)]))
    sp_dep_proj <- rbind(sp_dep_proj,sp_add)
  }
  
  colnames(sp_dep_proj) <- 2015:2100
  rownames(sp_dep_proj) <- unique(ensemble_avg_proj_sp$sp)
  
  rollmean_dep_proj <- apply(sp_dep_proj,1,zoo::rollmean,k=5,fill=NA,na.pad=TRUE,align = "right")
  rownames(rollmean_dep_proj) <- 2015:2100
  
  rollmean_dep_proj_filt <- rollmean_dep_proj %>%
    as.data.frame() %>%
    filter(rownames(rollmean_dep_proj) %in% seq(2025,2100,by=5))
  
  rollmean_dep_proj_filt_rn <- t(rollmean_dep_proj_filt) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var = "species")
  
  long_dep <- pivot_longer(rollmean_dep_proj_filt_rn,2:ncol(rollmean_dep_proj_filt_rn),
                           values_to = "dep",names_to = "year")
  
  sp_characteristics <- indexing_tab %>% 
    select(sp,mean_lat,mean_temp_niche,mean_precip_niche) %>%
    unique()
  
  sp_characteristics$species <- mdl_data$species %>% unique()
  sp_characteristics$mean_temp_niche <- sp_characteristics$mean_temp_niche - 273.15
  
  long_dep <- merge(long_dep,sp_characteristics,by="species")
  
  long_dep$year <- as.numeric(long_dep$year)
  
  quantile_segments[[ssp]] <- long_dep %>%
    group_by(year) %>%
    summarise(q10 = quantile(dep,.1),
              q90 = quantile(dep,.9))
}

quantile_segments[[1]]$scenario <- "SSP1-2.6"
quantile_segments[[2]]$scenario <- "SSP2-4.5"
quantile_segments[[3]]$scenario <- "SSP3-7.0"
quantile_segments_all <- do.call("rbind", quantile_segments)

pp3 <- ggplot()+
  geom_hline(yintercept = 0,linetype="dashed",col="darkgrey")+
  geom_segment(data = quantile_segments_all,aes(x=year-2,xend=year+2,y=q10,color = as.factor(scenario)),lwd=1.25)+
  geom_segment(data = quantile_segments_all,aes(x=year-2,xend=year+2,y=q90,color = as.factor(scenario)),lwd=1.25)+
  scale_color_manual(values=bcs_quantiles)+
  labs(col = "Scenario") +
  ylab("Departure Date Change (days)")+
  xlab("Year") + 
  scale_x_continuous(
    limits = xlim_all,
    breaks = seq(xlim_all[1],xlim_all[2],by=5),
    labels = function(x) paste0(x-4, "-", x),
    expand=c(0,0)
  ) +
  theme_bw(base_size = base_size_all) +
  theme(legend.key.size = unit(0.75, 'cm'),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position = c(0.075,0.2),
        legend.background = element_rect(fill=alpha("transparent")),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = xlim_all, ylim = c(-11, 7))
ggsave("plot3.pdf", plot = pp3, width = export_size_x, height = export_size_y, units = "in")   

