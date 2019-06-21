rm(list=ls()) # clean working environment

## set working directory
# setwd("1_monte-carlo")

# load rstan
library("rstan")
# parallelise if multiple cores/available RAM
options(mc.cores = parallel::detectCores())
# automatically save compiled stan program
rstan_options(auto_write = TRUE)

data <- list(N=20) # to be passed on to Stan

# fit Stan model. adapt_delta increased from 0.8 to 0.99 to reduce divergent
# transitions
drugMC.fit <- stan(file = 'drug-MC.stan', iter=10000, chains=1, data=data, 
                   control=list(adapt_delta=0.99), seed=2019)

# useful commands for quick overview
print(drugMC.fit)  # a rough summary
plot(drugMC.fit)   # a visual representation

# more plotting examples

# trace plot
traceplot(drugMC.fit, pars=c("theta"))

# package "bayesplot" produces nice plots of MCMC draws
library("bayesplot")
mcmc_hist(as.array(drugMC.fit), pars=c("theta"), binwidth=0.01)
mcmc_dens(as.array(drugMC.fit), pars=c("theta"))
