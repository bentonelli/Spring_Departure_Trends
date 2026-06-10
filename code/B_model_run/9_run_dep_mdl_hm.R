library(rstan)

setwd("Departure/")

data_mdl <- readRDS("list_format_data_ss.rds")

var_fit <- stan(
  file = "dep_clim.stan",  # Stan program
  data = data_mdl,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 3000,          # number of warmup iterations per chain
  iter = 8000,            # total number of iterations per chain
  cores = 4, 
  thin = 1, 
  pars = c("beta1_raw","beta2_raw","b1s","b2s",'z_b1s',"z_b2s","yy_raw"),
  include = FALSE,
  save_warmup = FALSE,
  control = list(adapt_delta = .8)
)

saveRDS(var_fit,"departure_climate_fit_ss1.rds")

