library(MCMCvis)
library(ggplot2)

#Model fit and data
mdl_fit <- readRDS("data/output/departure_model/departure_climate_fit_ss2.rds")
mdl_data <- readRDS("data/output/departure_model/full_data_ss.rds")
indexing_tab <- readRDS("data/output/departure_model/indexing_tab_data_ss.rds")

yy_trues_est <- MCMCchains(mdl_fit,params=c("yy_true"))
nrow(yy_trues_est)
ncol(yy_trues_est)

yy_measured <- mdl_data$trailing_hm_mean
adj_trues <- mean(yy_measured)

pdf("PPC.pdf",width = 6, height = 6, units = "in")
plot(density(yy_measured,bw=2),main="PPC",lwd=5,ylim=c(0,0.015),xlab="Day")
for(nn in seq(1,5000,by=10)){
  points(density(yy_trues_est[nn,]+adj_trues,bw=2),type="l",col=alpha("firebrick3",.01),lwd=.5)
}
dev.off()