## SETUP
rm(list=ls()) 
# setwd("3_bcea")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load data
y <- 241945
n <- 493527
data <- list(y=y,n=n)

filein <- "Laplace.stan" # model file location
# parameter initialisation values (2 chains)
inits_det <- list(list(theta=.1),list(theta=.9)) # deterministic
inits_ran <- function(){list(theta=runif(1))} # random

laplace.fit <- stan(file=filein, iter=10000, warmup=4500, chains=2,
                    data=data, init=inits_det, control=list(adapt_delta=0.99),
                    seed=2019, thin=1)

# summary results
print(laplace.fit)

# histogram
library("bayesplot")
mcmc_hist(as.array(laplace.fit), pars=c("theta"), binwidth=0.0005)
