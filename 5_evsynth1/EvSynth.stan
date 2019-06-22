data {
  int<lower=1> S;
  int<lower=1> H;
  int<lower=0> r0[S];
  int<lower=0> r1[S];
  int<lower=1> n0[S];
  int<lower=1> n1[S];
  int<lower=0> x[H];
  int<lower=1> m[H];
  real mu_inf;
  real<lower=0> sigma_inf;
  real mu_l;
  real<lower=0> sigma_l;
  real mean_alpha;
  real<lower=0> sd_alpha;
  real mean_mu_delta;
  real<lower=0> sd_mu_delta;
  real mean_mu_gamma;  
  real<lower=0> sd_mu_gamma;
}

parameters {
  real gamma[H];
  real mu_delta;
  real mu_gamma;
  real<lower=0> sigma_delta;
  real<lower=0> sigma_gamma;
  real phi;
  real c_inf;
  real<lower=0> l;
  real delta[S];
  real alpha[S];
}

transformed parameters {
  real<lower=0,upper=1> beta[H];
  real<lower=0,upper=1> pi0[S];
  real<lower=0,upper=1> pi1[S];
  for (h in 1:H) {
    beta[h] = inv_logit(gamma[h]);  
  }
  for (s in 1:S) {
    pi0[s] = inv_logit(alpha[s]);
    pi1[s] = inv_logit(alpha[s]+delta[s]);
  }
}

model {
  // priors
  mu_delta ~ normal(mean_mu_delta, sd_mu_delta);
  mu_gamma ~ normal(mean_mu_gamma, sd_mu_gamma);
  sigma_delta ~ uniform(0,10);
  sigma_gamma ~ uniform(0,10);
  phi ~ normal(0,1000000); 
  // costs of influenza
  c_inf ~ normal(mu_inf, sigma_inf);
  // Length of time to recovery when infected by influenza
  l ~ lognormal(mu_l,sigma_l);
  // Evidence synthesis on incidence of influenza in "healthy" adults in placebo group  
  x ~ binomial(m, beta);  
  gamma ~ normal(mu_gamma, sigma_gamma);
  // Evidence synthesis for efffectiveness of NIs
  delta ~ normal(mu_delta, sigma_delta);
  alpha ~ normal(mean_alpha, sd_alpha);  
  r0 ~ binomial(n0,pi0);
  r1 ~ binomial(n1,pi1);
}

generated quantities {
  real<lower=0> rho;
  real<lower=0,upper=1> p1;
  real<lower=0,upper=1> p2;
  // odds Ratio of influenza under treatment with NIs
  rho = exp(mu_delta); 
  // estimated probability of influenza in "healthy adults" under t=0
  p1 = inv_logit(mu_gamma);
  // estimated probability of influenza in "healthy adults" under t=1
  p2 = inv_logit(mu_gamma + mu_delta);   
}
