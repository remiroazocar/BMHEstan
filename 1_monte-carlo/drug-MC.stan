data {
  int<lower=1> N;
}

parameters { 
  real theta; 
}

model { 
  // can also pass beta distribution parameters as data if desired
  theta ~ beta(9.2,13.8); // prior distribution
}

generated quantities {
  int Y;
  real Pcrit;
  Y = binomial_rng(N, theta); // sampling dist.
  Pcrit = step(Y-14.5); // =1 if y >= 15, 0 otherwise  
}
