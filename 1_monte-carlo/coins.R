rm(list=ls()) # clean working environment

## set working directory
# setwd("1_monte-carlo")

# load rstan
library("rstan")
# parallelise if multiple cores/available RAM
options(mc.cores = parallel::detectCores())
# automatically save compiled stan program
rstan_options(auto_write = TRUE)

data <- list(N=10, theta=0.5) # to be passed on to Stan

# fit Stan model. set Fixed param as model has no parameters 
coins.fit <- stan(file = 'coins.stan', iter=1000, chains=1, data=data, 
                  warmup=0, algorithm="Fixed_param", seed=2019)

# useful commands for quick overview
print(coins.fit, digits=3)  # a rough summary
plot(coins.fit)   # a visual representation

chain <- 1
as.array(coins.fit)[1:1000,chain,]  # array: element, chain, column (theta/deviance)

# Collect posterior samples across all chains:
Y <- extract(coins.fit)$Y
P8 <- extract(coins.fit)$P8

# example: plot histogram for Y

# some settings to make histogram look good (integer values)
breaks=seq(min(Y)-0.5, max(Y)+0.5, by=1)
hist(Y,breaks=breaks, freq=FALSE, xlim=c(0,10), 
     xlab="Y", ylab="P(Y)", main="Y posterior density")
