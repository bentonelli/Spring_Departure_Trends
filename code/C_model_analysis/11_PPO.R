library(MCMCvis)
library(ggplot2)
library(stats)

#Model fit and data
var_fit <- readRDS("data/output/departure_model/departure_climate_fit_ss2.rds")

var_fit2 <- var_fit %>% 
  shredder::stan_select("omega","sigma_mu_alpha","theta1","theta2","xi","alpha_sigma_alpha","beta_sigma_alpha",
                          "mu_eta","sigma_eta","zeta1","zeta2",
                          "sigma_beta1","sigma_beta2",
                          "mu_nu1","mu_nu2","mu_sigma","sigma_sigma")

var_fit <- NULL 

set.seed(1234)
iter_length <- 20000

sigma_b1s <- rgamma(iter_length,3,2);        
sigma_b2s <- rgamma(iter_length,3,2);        

omega <- rnorm(iter_length,0,10);
sigma_mu_alpha <- rnorm(iter_length,30,30);

theta1 <- rnorm(iter_length,0,10);
theta2 <- rnorm(iter_length,0,10);

xi <- rnorm(iter_length,0,10);

alpha_sigma_alpha <- rexp(iter_length,.1);
beta_sigma_alpha <- rexp(iter_length,.1);

mu_eta <- rnorm(iter_length,0,10);
sigma_eta <- rgamma(iter_length,2,1);

zeta1 <- rnorm(iter_length,0,10);
zeta2 <- rnorm(iter_length,0,10);

sigma_beta1 <- rnorm(iter_length,0.75,0.25);
sigma_beta2 <- rnorm(iter_length,0.75,0.25);

mu_nu1 <- rnorm(iter_length,0,5);
mu_nu2 <- rnorm(iter_length,0,5);

mu_sigma <- rnorm(iter_length,0,3);
sigma_sigma <- rgamma(iter_length,2,2);

PR <- cbind(omega,sigma_mu_alpha,theta1,theta2,xi,alpha_sigma_alpha,beta_sigma_alpha,
            mu_eta,sigma_eta,zeta1,zeta2,sigma_beta1,sigma_beta2,
            mu_nu1,mu_nu2,mu_sigma,sigma_sigma)

MCMCtrace(var_fit2,
          params=c("omega","sigma_mu_alpha","theta1","theta2","xi","alpha_sigma_alpha","beta_sigma_alpha",
                   "mu_eta","sigma_eta","zeta1","zeta2",
                   "sigma_beta1","sigma_beta2",
                   "mu_nu1","mu_nu2","mu_sigma","sigma_sigma"),
          priors = PR,PPO_out = TRUE,post_zm = FALSE)
