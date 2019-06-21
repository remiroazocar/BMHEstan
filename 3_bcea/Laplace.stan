data {
  int<lower=0> y;
  int<lower=1> n;
}

parameters { 
  real theta; 
}

model { 
  theta ~ uniform(0,1); 
  y ~ binomial(n, theta);
}

generated quantities {
}
