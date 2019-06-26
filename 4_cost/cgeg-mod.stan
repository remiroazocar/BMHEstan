data {
  int<lower=1> N1;
  int<lower=1> N2;
  real<lower=0> cost1[N1];
  real<lower=0> cost2[N2];
  real eff1[N1];
  real eff2[N2];
}

parameters {
  real<lower=0> shape_c1;
  real<lower=0> mu_c1;
  real beta_c1;
  real<lower=0> shape_e1;
  real mu_e1;
  real<lower=0> shape_c2;
  real<lower=0> mu_c2;
  real beta_c2;
  real<lower=0> shape_e2;
  real mu_e2;
}

transformed parameters {
  real phi_c1[N1];
  real phi_e1[N1];
  real<lower=0> rate_c1[N1];
  real<lower=0> rate_e1[N1];
  real phi_c2[N2];
  real phi_e2[N2];
  real<lower=0> rate_c2[N2];
  real<lower=0> rate_e2[N2];  
  for (i in 1:N1) {
    phi_c1[i] = mu_c1 + beta_c1*(eff1[i]-mu_e1);
    phi_e1[i] = mu_e1;
    rate_c1[i] = shape_c1/phi_c1[i];
    rate_e1[i] = shape_e1/phi_e1[i];
  }
  for (i in 1:N2) {
    phi_c2[i] = mu_c2 + beta_c2*(eff2[i]-mu_e2);
    phi_e2[i] = mu_e2;
    rate_c2[i] = shape_c2/phi_c2[i];
    rate_e2[i] = shape_e2/phi_e2[i];
  }
}


model {
  // arm 1 prior distributions
  shape_c1 ~ uniform(0,10);
  mu_c1 ~ uniform(0,1000);
  beta_c1 ~ uniform(-10,10);
  shape_e1 ~ uniform(0,10);
  mu_e1 ~ uniform(0,1000);
  cost1 ~ gamma(shape_c1, rate_c1);
  eff1 ~ gamma(shape_e1, rate_e1);
  // arm 2 prior distributions
  shape_c2 ~ uniform(0,10);
  mu_c2 ~ uniform(0,1000);
  beta_c2 ~ uniform(-10,10);
  shape_e2 ~ uniform(0,10);
  mu_e2 ~ uniform(0,1000);
  cost2 ~ gamma(shape_c2, rate_c2);
  eff2 ~ gamma(shape_e2, rate_e2);
}

generated quantities {
  real<lower=0> mu_c[2];
  real mu_e[2];
  real delta_c;
  real delta_e;
  real<lower=0> K_space; // CEAC curves: space between willingness-to-pay tick marks
  real K[11];
  real INB[11];
  real CEAC[11];  
  // HE outcomes
  mu_c[1] = mu_c1;
  mu_c[2] = mu_c2;
  mu_e[1] = -mu_e1;
  mu_e[2] = -mu_e2;
  delta_c = mu_c2-mu_c1;
  delta_e = -(mu_e2-mu_e1); // as days in hospital is a bad thing
  // CEAC curves
  K_space = 0.1; // 100 pounds
  for (j in 1:11) {
    K[j] = (j-1)*K_space;
    INB[j] = K[j]*delta_e - delta_c;
    CEAC[j] = step(INB[j]);
  }  
}
