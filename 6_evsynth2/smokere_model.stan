// Based on program 1(c) of NICE DSU TSD 2

data {
  int ns; 
  int nt; 
  int na[ns];
  int max_na;
  int r[ns,max_na];
  int n[ns,max_na];
  int t[ns,max_na];
}

parameters {
  real mu[ns];
  real d_other[nt-1]; // treatments beside reference, which is fixed to zero
  real delta_other[ns,max_na-1];
  real<lower=0> stdev; // between trial standard deviation
  real alpha;
  real L;  
}

transformed parameters {
  real p[ns,max_na];
  real w[ns,max_na]; 
  real sw[ns,max_na-1]; // no entry for reference arm in these
  real md[ns,max_na-1]; 
  real taud[ns,max_na-1];
  real sdd[ns,max_na-1];
  real<lower=0> tau; // between-trial precision
	real d_1;
	real d[nt];
	real delta_1[ns];
	real delta[ns,max_na];
	tau = pow(stdev,-2); // between-trial precision = (1/between-trial variance)
	d_1 = 0; // treatment effect is zero for reference treatment
	d[1] = d_1;
	d[2:nt] = d_other;
	delta_1 = rep_array(0,ns); // treatment effect is zero for control arm
	delta[,1] = delta_1; 
	delta[,2:max_na] = delta_other;
	for (i in 1:ns) { // LOOP TROUGH STUDIES
	  w[i,1] = 0; // adjustment for multi-arm trials is zero for control arm
	  for (k in 1:na[i]) { // LOOP TROUGH ARMS 
	    p[i,k] = inv_logit(mu[i] + delta[i,k]);
	  }
	  for (k in 2:na[i]) { // LOOP TROUGH ARMS 
	    w[i,k] = (delta[i,k] - d[t[i,k]] + d[t[i,1]]); // adjustment for multi-arm RCTs
	    sw[i,k-1] = sum(w[i,1:(k-1)])/(k-1); // cumulative adjustment for multi-arm trials
	    md[i,k-1] = d[t[i,k]] - d[t[i,1]] + sw[i,k-1]; // mean of LOR distributions (with multi-arm trial correction)
	    taud[i,k-1] = tau*2*(k-1)/k; // precision of LOR distributions (with multi-arm trial correction)
	    sdd[i,k-1] = pow(taud[i,k-1],-0.5);
	  }
	}
}

model {
  // Binomial likelihood, logit link
  // Random effects model for multi-arm trials
  // priors changed from super-vague to weakly informative for slightly faster convergence
	mu ~ cauchy(0,5); // priors for all trial baselines 
  d_other ~ cauchy(0,5); // priors for treatment effects
  stdev ~ inv_gamma(0.1, 0.1); // prior for between-trial SD.
	for(i in 1:ns) { // LOOP THROUGH STUDIES 
		for (k in 1:na[i]) { // LOOP TROUGH ARMS 
			r[i,k] ~ binomial(n[i,k], p[i,k]); // binomial likelihood
		} 
		for (k in 1:(na[i]-1)) {
			delta_other[i,k] ~ normal(md[i,k], sdd[i,k]); // trial-specific LOR distributions
		}
	}  
  // Log odds of quitting successfully under no intervention (from published data)
	alpha ~ normal(-2.6, 0.38); // precision=6.925
	// Life years gained by quitting
	L ~ normal(15,4); // precision=0.0625	
}

generated quantities {
	real o_r[nt,nt];
	real pq[nt];
	// odds ratios for all treatment comparisons
	for (c in 1:(nt-1)) {
	  o_r[c,c] = 1;
	  for (k in (c+1):nt)  {
		  o_r[c,k] = exp(d[c]-d[k]);
		  o_r[k,c] = 1/o_r[c,k];
		}
	}	
	o_r[nt,nt] = 1;
	// Absolute probability of quitting successfully under each intervention
	for (i in 1:nt) {
		pq[i] = inv_logit(alpha + d[i]);
	}
}
