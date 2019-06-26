data {
  int<lower=1> N1;
  int<lower=1> N2;
  real<lower=0> cost1[N1];
  real<lower=0> cost2[N2];
}

parameters {
  real mu[2];
  real ls[2]; // log-st.dev
}

transformed parameters {
  real delta_c;
  real<lower=0> sigma[2];
  real<lower=0> ss[2];
  delta_c = mu[2] - mu[1];
  // compute standard deviation from its log
  sigma[1] = 10 ^ ls[1];
  sigma[2] = 10 ^ ls[2];
  // compute variance from standard deviation
  ss[1] = sigma[1]*sigma[1];
  ss[2] = sigma[2]*sigma[2];
}

model {
  // fit Normal model to observed costs from arms 1 and 2
  for (i in 1:N1) {
    cost1[i] ~ normal(mu[1], sigma[1]);    
  }
  for (i in 1:N2) {
    cost2[i] ~ normal(mu[2], sigma[2]);   
  }
  // priors
  mu ~ uniform(0,1000000);
  ls ~ uniform(-5,5);
}

generated quantities {
  real d1[N1];
  real d2[N2];
  real dev[2];
  real total_dev;
  real pred_cost1[N1];
  real pred_cost2[N2];
  real<lower=0, upper=1> p_cost1[N1];
  real<lower=0, upper=1> p_cost2[N2];
  real total_cost;
  real cum[50];
  real cum_step[50];
  real<lower=0,upper=1> por;
  // Deviance
  for (i in 1:N1) {
    d1[i] = log(exp(-pow(cost1[i]-mu[1],2)/(2*ss[1]))/(sqrt(2*3.1415927*ss[1])));
  } 
  for (i in 1:N2) {
    d2[i] = log(exp(-pow(cost2[i]-mu[2],2)/(2*ss[2]))/(sqrt(2*3.1415927*ss[2])));
  }
  dev[1] = -2*sum(d1);
  dev[2] = -2*sum(d2);
  total_dev = sum(dev);
  for(i in 1:N1){
    // predict costs for a set of future patients
    pred_cost1[i] = normal_rng(mu[1],sigma[1]);
    // predictive p-values to measure "extremeness" of observed costs
    p_cost1[i] = step(pred_cost1[i]-cost1[i]);
  }
  for(i in 1:N2){
    pred_cost2[i] = normal_rng(mu[2],sigma[2]);
    p_cost2[i] = step(pred_cost2[i]-cost2[i]); 
  }
  // predicted total cost for 10 patients in group 2
  total_cost = sum(pred_cost2[1:10]);
  // cumulate cost of treating patients with treatment 2
  cum[1] = pred_cost2[1];
  for(i in 2:50) {
    cum[i] = cum[i-1]+pred_cost2[i];  // cumulative cost for i patients
  }
  // cum_step = i if cost for i patients < Budget, 0 otherwise, e.g. if cum[5] < 200
  // but cum[6] > 200, cum_step = 1,2,3,4,5,0,0,0,....
  for (i in 1:50) {
    cum_step[i]=i*step(200-cum[i]);
  }
  // Probability of ruin given a budget of 300
  por = step(total_cost-300);
}
