library(dplyr)
library(ggplot2)
library(dggridR)
library(sf)
library(terra)
library(MCMCvis)
library(mapdata)
library(birdcolors)

#Model fit and data
mdl_fit <- readRDS("data/output/departure_model/departure_climate_fit_ss2.rds")
mdl_data <- readRDS("data/input/mdl_data/full_data_ss.rds")
indexing_tab <- readRDS("data/input/mdl_data/indexing_tab_data_ss.rds")

site_enviro_metrics <- mdl_data %>% 
  select(sp_site,sd_temp_60,sd_total_precip60) %>% 
  unique()

#Maps
alphas <- MCMCsummary(mdl_fit,params=c("alpha"),pg0 = TRUE)
betas1 <- MCMCsummary(mdl_fit,params=c("beta1"),pg0 = TRUE)
betas2 <- MCMCsummary(mdl_fit,params=c("beta2"))

post_ind <- indexing_tab
post_ind$alphas <- alphas$mean
post_ind$betas_mean_temp <- betas1$mean
post_ind$betas_mean_precip <- betas2$mean

post_ind <- merge(post_ind,site_enviro_metrics,by="sp_site")

#post_ind$shift_per_degree <- post_ind$betas_mean_temp * (1/post_ind$sd_temp_60)
post_ind$species <- unique(mdl_data$species)[post_ind$sp]

col_fact_temp <- rep(2,nrow(post_ind))
col_fact_temp[which(post_ind$betas_mean_temp > 0.5)] <- 3
col_fact_temp[which(post_ind$betas_mean_temp < -0.5)] <- 1

col_fact_precip <- rep(2,nrow(post_ind))
col_fact_precip[which(post_ind$betas_mean_precip > 0.5)] <- 3
col_fact_precip[which(post_ind$betas_mean_precip < -0.5)] <- 1

col_opt_temp <- bird_colors(palette_name = "Bluethroat")[2:4]
col_opt_precip <- bird_colors(palette_name = "Ruby-throated Hummingbird",rev=T)[2:4]
col_opt_precip[2] <- col_opt_temp[2]

post_ind$col_fact_temp <- as.factor(col_fact_temp)
post_ind$col_fact_precip <- as.factor(col_fact_precip)

### GENERATE X-Y LOCATIONS ####
grid_elev_bands <- st_read("data/input/spatial/grid_elev_bands.shp")

# Get unique sites in phenology dataset
unq_cl_elev <- mdl_data %>% 
  select(cell,elev,species) %>%
  unique() %>%
  group_by(cell,elev) %>%
  summarise(count = n())
  
# Merge with elev, cell
grid_elev_bands <- merge(grid_elev_bands,unq_cl_elev,by=c("elev","cell"))

#Get number of populations for each cell

samp_point_test <- st_sample(grid_elev_bands,grid_elev_bands$count,type="random",exact=T)

# This could be done in a much simpler way...
scaffold <- c()
for (nn in 1:nrow(grid_elev_bands)){
  for(mm in 1:as.numeric(grid_elev_bands$count[nn])){
    scaffold <- rbind(scaffold,rep(as.numeric(grid_elev_bands[nn,])[1:2]))
  }
}

scaffold <- as.data.frame(scaffold)
colnames(scaffold) <- c("elev","cell")

#plot(test_poly)
plot(samp_point_test,cex=.1,pch=19)

sum(grid_elev_bands$count)
length(samp_point_test)

point_coord <- as.data.frame(st_coordinates(samp_point_test))

point_coord$elev <- scaffold$elev
point_coord$cell <- scaffold$cell

cl_elev_sp <- post_ind %>% 
  arrange(elev,cell,species)

point_coord <- point_coord %>% arrange(elev,cell)

cl_elev_sp$X <- point_coord$X
cl_elev_sp$Y <- point_coord$Y

water_col = "grey40"
country_col = "grey20"
midpoint_col = "grey90"
#grid_outline_col = "grey10"
country_outline_col = "grey50"
countries <- map_data("world") 
states <- map_data("state") 
lakes <- map_data("lakes") 
see_thru <- .8
point_size <- 2
pch_points <- 19

panel_size <- 6

#Global
map_res_x <- c(-179,179)
map_res_y <- c(-70,70)

p_departure <- ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group),alpha=1,fill=country_col) +
  geom_path(data=countries, aes(x=long, y=lat, group=group),alpha=1,cex=.1,color=country_outline_col) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),alpha=1,fill=water_col) +
  geom_path(data=lakes, aes(x=long, y=lat, group=group),alpha=1,cex=.05,color=country_outline_col) +
  geom_point(data=cl_elev_sp,aes(x=X,y=Y,col=alphas+mean(mdl_data$trailing_hm_mean)),alpha=see_thru,size=0.5,pch=pch_points) + 
  scale_color_bird(bird_cols = bird_colors("Eastern Bluebird",rev=F)) +
  labs(col="Mean Departure") +
  coord_sf(xlim=map_res_x,ylim=map_res_y) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(1.25, 'cm'),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12),
        legend.position = c(.2,.4),
        legend.background = element_rect(fill=alpha(water_col))) +
  theme(panel.background = element_rect(fill = alpha(water_col)))

p_departure
ggsave("~/Desktop/global_average_departure.pdf",p_departure,width = 24,height = 8)

p_temp <- ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group),alpha=1,fill=country_col) +
  geom_path(data=countries, aes(x=long, y=lat, group=group),alpha=1,cex=.1,color=country_outline_col) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),alpha=1,fill=water_col) +
  geom_path(data=lakes, aes(x=long, y=lat, group=group),alpha=1,cex=.05,color=country_outline_col) +
  geom_point(data=cl_elev_sp,aes(x=X,y=Y,col=col_fact_temp),alpha=see_thru,size=point_size*.25,pch=pch_points) + 
  scale_color_manual(values=col_opt_temp) + 
  coord_sf(xlim=map_res_x,ylim=map_res_y) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(1.25, 'cm'),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12),
        legend.position = "none",
        legend.background = element_rect(fill=alpha(water_col))) +
  theme(panel.background = element_rect(fill = alpha(water_col)))

p_temp
ggsave("~/Desktop/global_temp.pdf",p_temp,width = 24,height = 8)

#Global Precip
p_precip <- ggplot() +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group),alpha=1,fill=country_col) +
  geom_path(data=countries, aes(x=long, y=lat, group=group),alpha=1,cex=.1,color=country_outline_col) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),alpha=1,fill=water_col) +
  geom_path(data=lakes, aes(x=long, y=lat, group=group),alpha=1,cex=.05,color=country_outline_col) +
  geom_point(data=cl_elev_sp,aes(x=X,y=Y,col=col_fact_precip),alpha=see_thru,size=point_size*.25,pch=pch_points) + 
  scale_color_manual(values=col_opt_precip) + 
  coord_sf(xlim=map_res_x,ylim=map_res_y) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        legend.key.size = unit(1.25, 'cm'),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12),
        legend.position = "none",
        legend.background = element_rect(fill=alpha(water_col))) +
  theme(panel.background = element_rect(fill = alpha(water_col)))
p_precip
ggsave("~/Desktop/global_precip.pdf",p_precip,width = 24,height = 8)

# Temp Panels
#NA Temp 70 by 65
map_res_x <- c(-125,-55)
map_res_y <- c(-5,60)
p_temp_mod <- p_temp
p_temp_mod$layers[[5]]$aes_params$size <- .25
p_temp_na <- p_temp_mod + coord_fixed(ratio=1,xlim=map_res_x,ylim=map_res_y)
ggsave("~/Desktop/temp_na.pdf",p_temp_na,width = panel_size,height = panel_size)

#Europe Temp 70 by 65
map_res_x <- c(-10,25)
map_res_y <- c(30,62.5)
p_temp_mod <- p_temp
p_temp_mod$layers[[5]]$aes_params$size <- point_size
p_temp_eu <- p_temp_mod + coord_fixed(ratio=1,xlim=map_res_x,ylim=map_res_y)
ggsave("~/Desktop/temp_eu.pdf",p_temp_eu,width = panel_size,height = panel_size)

#S/SE Asia Temp 70 x 60
map_res_x <- c(70,105)
map_res_y <- c(2,34.5)
p_temp_mod <- p_temp
p_temp_mod$layers[[5]]$aes_params$size <- point_size
p_temp_sse <- p_temp_mod + coord_fixed(ratio=1,xlim=map_res_x,ylim=map_res_y)
ggsave("~/Desktop/temp_sse.pdf",p_temp_sse,width = panel_size,height = panel_size)

#Precip Panels
#NA Temp 70 by 65
map_res_x <- c(-125,-55)
map_res_y <- c(-5,60)
p_precip_mod <- p_precip
p_precip_mod$layers[[5]]$aes_params$size <- .25
p_precip_na <- p_precip_mod + coord_fixed(ratio=1,xlim=map_res_x,ylim=map_res_y)
ggsave("~/Desktop/precip_na.pdf",p_precip_na,width = panel_size,height = panel_size)

#Europe precip 70 by 65
map_res_x <- c(-10,25)
map_res_y <- c(30,62.5)
p_precip_mod <- p_precip
p_precip_mod$layers[[5]]$aes_params$size <- point_size
p_precip_eu <- p_precip_mod + coord_fixed(ratio=1,xlim=map_res_x,ylim=map_res_y)
ggsave("~/Desktop/precip_eu.pdf",p_precip_eu,width = panel_size,height = panel_size)

#S/SE Asia precip 70 x 60
map_res_x <- c(70,105)
map_res_y <- c(2,34.5)
p_precip_mod <- p_precip
p_precip_mod$layers[[5]]$aes_params$size <- point_size
p_precip_sse <- p_precip_mod + coord_fixed(ratio=1,xlim=map_res_x,ylim=map_res_y)
ggsave("~/Desktop/precip_sse.pdf",p_precip_sse,width = panel_size,height = panel_size)



