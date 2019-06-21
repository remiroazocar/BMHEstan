data {
  real<lower=0> a;
  real<lower=0> b;
  int<lower=1> m;
  int<lower=1> n;
  int<lower=0> ncrit;
}

parameters { 
  real theta; 
}

model { 
  theta ~ beta(a,b); // prior distribution
}

generated quantities {
  int y;
  int ypred;
  real Pcrit;
  y = binomial_rng(m, theta); // sampling dist.
  ypred = binomial_rng(n, theta); // predictive dist.
  Pcrit = step(ypred - ncrit + 0.5); // =1 if ypred>=ncrit, 0 otherwise  
}
