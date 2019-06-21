## SETUP
rm(list=ls()) 
# setwd("2_mcmc")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# data in separate list to be passed on to Stan (usually better)
data <- list(
            a = 9.2, b = 13.8, # prior parameters
            m = 20, # number of trials
            n = 40, # future number of trials
            ncrit = 25) # critical value of future successes

# initial parameter values to be passed on to Stan
init <- list(list(theta=0.1))

drugMCMC.fit <- stan(file = 'drug-MCMC.stan', iter=10000, chains=1, 
                     data=data, init=init, control=list(adapt_delta=0.99), 
                     seed=2019)

# more useful plots. these are diagnostic plots
library("bayesplot")

# rhat
rhats <- rhat(drugMCMC.fit) 
mcmc_rhat(rhats)

# effective sample size
ratios <- neff_ratio(drugMCMC.fit)
mcmc_neff(ratios)

# autocorrelation
draws <- as.array(drugMCMC.fit)
mcmc_acf(draws, pars = c("theta"), lags = 10)
