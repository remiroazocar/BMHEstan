data {
  int<lower=1> N;
  real<lower=0,upper=1> theta;
}

parameters { // no parameters
}

model { // no parameters
}

generated quantities {
  int Y;
  real P8;
  Y = binomial_rng(N, theta); 
  P8 = step(Y-7.5); // = 1 if Y is 8 or more
}
