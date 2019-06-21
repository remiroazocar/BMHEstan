data {
  real<lower=0> a;
  real<lower=0> b;
  int<lower=0> y;
  int<lower=1> m;
  int<lower=1> n;
  int<lower=0> ncrit;
}

parameters { 
  real theta; 
}

model { 
  theta ~ beta(a,b); // prior distribution
  y ~ binomial(m, theta); // sampling dist.
}

generated quantities {
  int ypred;
  real Pcrit;
  ypred = binomial_rng(n, theta); // predictive dist.
  Pcrit = step(ypred - ncrit + 0.5); // =1 if ypred>=ncrit, 0 otherwise  
}
