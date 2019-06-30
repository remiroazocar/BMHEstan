data {
  int ns; 
  int nt; 
  int na[ns];
  int r[ns,max(na)];
  int n[ns,max(na)];
  int t[ns,max(na)];
}

parameters {
  real mu[ns];
  real d_other[nt-1]; // treatments beside reference, which is fixed to zero
  real alpha;
  real L;
}

transformed parameters {
	real delta[ns,max(na)];
  real p[ns,max(na)];
	real d_1;
	real d[nt];
	d_1 = 0; // treatment effect is zero for reference treatment
	d[1] = d_1;
	d[2:nt] = d_other;
	for (i in 1:ns) { // LOOP TROUGH STUDIES
	  for (k in 1:(na[i])) { // LOOP TROUGH ARMS 
	    delta[i,k] = d[t[i,k]] - d[t[i,1]]; // lin. predictor model
	    p[i,k] = inv_logit(mu[i] + delta[i,k]);
	  }   
	}    
}

model {
  // Binomial likelihood, logit link
  // Fixed effect model
	mu ~ normal(0,1000000); // vague priors for all trial baselines 
  d_other ~ normal(0,1000000); //  vague priors for treatment effects
	for(i in 1:ns) { // LOOP THROUGH STUDIES 
		for (k in 1:na[i]) { // LOOP TROUGH ARMS 
			r[i,k] ~ binomial(n[i,k], p[i,k]); // binomial likelihood
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
	for (c in 1:nt) {
	    o_r[c,c] = 1;
	}
	for (c in 1:(nt-1)) {
	  for (k in (c+1):nt)  {
		  o_r[c,k] = exp(d[c]-d[k]);
		  o_r[k,c] = 1/o_r[c,k];
		}
	}	
	// Absolute probability of quitting successfully under each intervention
	for (i in 1:nt) {
		pq[i] = inv_logit(alpha + d[i]);
	}
}
