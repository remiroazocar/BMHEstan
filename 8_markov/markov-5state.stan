data {
  // changed MxN to NxM in order to facilitate Stan operations with vectors
  int<lower=0> r_sfc[5,4];
  int<lower=0> n_sfc[4];
  vector[4] weekly_cost_sfc;
  int<lower=0> r_fp[5,4];
  int<lower=0> n_fp[4];
  vector[4] weekly_cost_fp;
  matrix<lower=0>[5,4] prior_sfc;
  matrix<lower=0>[5,4] prior_fp;
  vector<lower=0,upper=1>[5] p_fixed;
  row_vector<lower=0,upper=1>[5] s_start;
}

parameters {
  // cannot assigned fixed values to model parameter
  // split p_sfc, p_fp into non-absorbing/absorbing components and assign fixed val in transformed params.
  simplex[5] p_sfc_nonabs[4];
  simplex[5] p_fp_nonabs[4];
}

transformed parameters {
  vector<lower=0,upper=1>[5] p_sfc_abs;
  vector<lower=0,upper=1>[5] p_fp_abs;
  // Fixed transition probabilites for the absorbing state
  p_sfc_abs = p_fixed;
  p_fp_abs = p_fixed;
}

// vectorise
// probabilities parameter is not a valid simplex
model {
  for (i in 1:4) {
    // multinomial distribution for r for the i=4 non-absorbing states  
    r_sfc[,i] ~ multinomial(p_sfc_nonabs[i]);
    r_fp[,i] ~ multinomial(p_fp_nonabs[i]);
    // Dirichlet prior distributions for the transition probabilities
    p_sfc_nonabs[i] ~ dirichlet(prior_sfc[,i]);
    p_fp_nonabs[i] ~ dirichlet(prior_fp[,i]);
  }
}

generated quantities {
  matrix[5,4] p_sfc_nonabs_mat; // convert simplex to matrices for matrix operations
  matrix[5,4] p_fp_nonabs_mat;
  matrix[5,5] p_sfc;
  matrix[5,5] p_fp;
  matrix<lower=0,upper=1>[13,5] s_sfc;
  matrix<lower=0,upper=1>[13,5] s_fp;
  vector[13] weekly_cost_tf_fp;
  vector[13] cost_sfc;
  vector[13] cost_fp;
  real total_cost_sfc; 
  real total_cost_fp;
  real total_cost_diff;
  real stw_sfc;
  real stw_fp;
  real stw_diff;
  real<lower=0> K_space;
  real K[11];
  real INB[11];
  real Q[11];  
  real ce[4];
  // Starting state
  s_sfc[1,] = s_start;
  s_fp[1,] = s_start;
  // concatenate non-absorbing and absorbing probabilities
  for (m in 1:size(p_sfc_nonabs)) {
    p_sfc_nonabs_mat[,m]=p_sfc_nonabs[m];
    p_fp_nonabs_mat[,m]=p_fp_nonabs[m];
  }
  p_sfc = append_col(p_sfc_nonabs_mat, p_sfc_abs);
  p_fp = append_col(p_fp_nonabs_mat, p_fp_abs);
  // Markov model
  for (t in 2:13) {
    s_sfc[t,] = s_sfc[t-1,]*p_sfc'; // state proportions
    s_fp[t,] = s_fp[t-1,]*p_fp';
    weekly_cost_tf_fp[t] = (s_fp[t,1:4]*weekly_cost_fp)/sum(s_fp[t,1:4]); // cost of TF at time t
  }
  cost_sfc = s_sfc[,1:4]*weekly_cost_sfc+s_sfc[,5].*weekly_cost_tf_fp; // cost at time t
  cost_fp = s_fp[,1:4]*weekly_cost_fp+s_fp[,5].*weekly_cost_tf_fp;
  // Costs and utilities
  total_cost_sfc = sum(cost_sfc[2:13]);
  total_cost_fp = sum(cost_fp[2:13]);
  total_cost_diff = total_cost_sfc-total_cost_fp;
  stw_sfc = sum(s_sfc[2:13,1]);
  stw_fp = sum(s_fp[2:13,1]);
  stw_diff = stw_sfc-stw_fp; 
  // INB and CEAC
  K_space = 5;
  for (j in 1:11) {
    K[j] = (j-1)*K_space;
    INB[j] = K[j]*stw_diff-total_cost_diff;
    Q[j] = step(INB[j]);
  }  
  // ce[1]=c1, ce[2]=c2, ce[3]=e1, ce[4]=e2
  ce[1] = total_cost_fp;
  ce[2] = total_cost_sfc;
  ce[3] = stw_fp;
  ce[4] = stw_sfc;
}  