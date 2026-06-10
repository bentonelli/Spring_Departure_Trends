library(birdcolors)
library(MCMCvis)

#Model fit and data
mdl_fit <- readRDS("data/output/departure_model/departure_climate_fit_ss2.rds")

bc_temp <- bird_colors("Lazuli Bunting")
bc_precip <- bird_colors("Costa's Hummingbird")

MCMCplot(mdl_fit,params=c("theta1","theta2","zeta1","zeta2","mu_nu1","mu_nu2","mu_eta","xi"),
         col=c(bc_temp[3],bc_precip[1],bc_temp[3],bc_precip[1],bc_temp[3],bc_precip[1],"black","black"),
         labels=c(expression(theta[1]),expression(theta[2]),
                  expression(zeta[1]),expression(zeta[2]),
                  expression(mu[nu[1]]),expression(mu[nu[2]]),
                  expression(mu[eta]),expression(xi)))

         