
data {
  int<lower=0> N;
  vector[N] yy;
  vector[N] yy_unc;
  vector[N] temp;
  vector[N] precip;

  int Nsp;
  int<lower=0,upper=Nsp> sp [N];
  
  int Nsp_c;
  int<lower=0,upper=Nsp_c> sp_c [N];
  int Nsite;
  
  int<lower=0,upper=Nsp> sp_match [Nsp_c];
  int<lower=0,upper=Nsp_c> sp_cl_match [Nsp_c];
  int<lower=0,upper=Nsite> site_match [Nsp_c];
  
  vector[Nsp_c] lat;

  vector[Nsp_c] TN;
  vector[Nsp_c] PN;
  
  vector[Nsp] MT;
  vector[Nsp] MP;
  
  vector[Nsp] ML;
}

parameters {
  
  real omega;
  real<lower=0> sigma_mu_alpha;
  
  real theta1;
  real theta2;

  vector[Nsp] mu_alpha;
  
  real<lower=0>alpha_sigma_alpha;
  real<lower=0>beta_sigma_alpha;
  vector<lower=0>[Nsp] sigma_alpha;
  
  real mu_eta;
  real<lower=0> sigma_eta;
  vector[Nsp] eta;
  
  real mu_nu1;
  real mu_nu2;
  
  real zeta1;
  real zeta2;
  
  real xi;
  
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_beta2;
  
  vector[Nsp_c] alpha;
  
  real mu_sigma;
  real<lower=0> sigma_sigma;
  vector<lower=0>[Nsp] sigma;

  //Non centered stuff
  vector[N] yy_raw;
  
  vector[Nsp_c] beta1_raw;
  vector[Nsp_c] beta2_raw;
  
  vector<lower=0>[2] sigma_b1s;
  cholesky_factor_corr[2] L_Rho_b1s;             
  matrix[Nsp, 2] z_b1s;                          
  
  vector<lower=0>[2] sigma_b2s;
  cholesky_factor_corr[2] L_Rho_b2s;             
  matrix[Nsp, 2] z_b2s;                          
  
}

transformed parameters {
  vector[N] yy_true;
  
  matrix[Nsp, 2] b1s;                 
  matrix[2, 2] Rho_b1s;               // corr matrix
  
  matrix[Nsp, 2] b2s;                 
  matrix[2, 2] Rho_b2s;               // corr matrix
  
  vector[Nsp] mu_beta1;
  vector[Nsp] mu_beta2;
  
  vector[Nsp] nu1;
  vector[Nsp] nu2;
  
  vector[Nsp_c] beta1;
  vector[Nsp_c] beta2;
  
  b1s = z_b1s * diag_pre_multiply(sigma_b1s, L_Rho_b1s)';
  Rho_b1s = multiply_lower_tri_self_transpose(L_Rho_b1s);
  
  b2s = z_b2s * diag_pre_multiply(sigma_b2s, L_Rho_b2s)';
  Rho_b2s = multiply_lower_tri_self_transpose(L_Rho_b2s);
  
  mu_beta1 = theta1 + zeta1 * MT + b1s[,1];
  nu1 = mu_nu1 + b1s[,2];
  
  mu_beta2 = theta2 + zeta2 * MP + b2s[,1];
  nu2 = mu_nu2 + b2s[,2];

  beta1 = mu_beta1[sp_match] + nu1[sp_match] .* TN + sigma_beta1 * beta1_raw;
  beta2 = mu_beta2[sp_match] + nu2[sp_match] .* PN + sigma_beta2 * beta2_raw;
  
  yy_true = alpha[sp_c] + beta1[sp_c] .* temp + beta2[sp_c] .* precip + sigma[sp] .* yy_raw;
}

model {
  sigma_b1s ~ gamma(2,2);        
  sigma_b2s ~ gamma(2,2);        
  
  omega ~ normal(0,10);
  sigma_mu_alpha ~ normal(30,30);
  
  theta1 ~ normal(0,10);
  theta2 ~ normal(0,10);

  beta1_raw ~ normal(0,1);
  beta2_raw ~ normal(0,1);
  
  xi ~ normal(0,10);
  
  mu_alpha ~ normal(omega + xi * ML,sigma_mu_alpha);
  
  alpha_sigma_alpha ~ exponential(.1);
  beta_sigma_alpha ~ exponential(1);
  sigma_alpha ~ gamma(alpha_sigma_alpha,beta_sigma_alpha);

  mu_eta ~ normal(0,10);
  sigma_eta ~ gamma(2,1);
  eta ~ normal(mu_eta,sigma_eta);
  
  zeta1 ~ normal(0,10);
  zeta2 ~ normal(0,10);
  
  sigma_beta1 ~ normal(0.75,0.25);
  sigma_beta2 ~ normal(0.75,0.25);
  
  mu_nu1 ~ normal(0,5);
  mu_nu2 ~ normal(0,5);
  
  alpha ~ normal(mu_alpha[sp_match] + eta[sp_match] .* lat,sigma_alpha[sp_match]);
  
  mu_sigma ~ normal(0,3);
  sigma_sigma ~ gamma(2,2);
  sigma ~ lognormal(mu_sigma,sigma_sigma);
  
  
  yy_raw ~ normal(0,1);
  yy ~ normal(yy_true,yy_unc);
  
  // MVN
  to_vector(z_b1s) ~ std_normal();
  L_Rho_b1s ~ lkj_corr_cholesky(1);
  
  to_vector(z_b2s) ~ std_normal();
  L_Rho_b2s ~ lkj_corr_cholesky(1);
  
}
