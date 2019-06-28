# Setup
rm(list=ls()) 
library("rstan")
# setwd("8_markov")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Specify Stan options
model.file <- "markov-5state.stan"
# Data specified in markov-5-state-data.txt
r_sfc <- matrix(c(210,60,0,1,1, 88,641,0,4,13, 0,0,0,0,0, 1,0,0,0,1), 
                ncol=4, byrow=FALSE)
n_sfc <- c(272,746,0,2)
weekly_cost_sfc <- c(7.96,7.96,1821.17,100.79)
r_fp <- matrix(c(66,32,0,0,2, 42,752,0,5,20, 0,0,0,0,0, 0,4,0,1,0), 
               ncol=4, byrow=FALSE)
n_fp <- c(100,819,0,5)
weekly_cost_fp <- c(2.38,2.38,1851.58,95.21)
p_fixed <- c(0,0,0,0,1)
s_start <- c(1,0,0,0,0)
prior_sfc <- matrix(rep(1,20), ncol=4, byrow=FALSE)
prior_fp <- matrix(rep(1,20), ncol=4, byrow=FALSE)
data <- list(r_sfc=r_sfc,n_sfc=n_sfc,weekly_cost_sfc=weekly_cost_sfc,
             r_fp=r_fp,n_fp=n_fp,weekly_cost_fp=weekly_cost_fp,
             p_fixed=p_fixed,s_start=s_start,prior_sfc=prior_sfc,prior_fp=prior_fp)
params <- c("total_cost_sfc","total_cost_fp","stw_sfc","stw_fp","stw_diff",
            "INB","p_sfc","p_fp","Q","ce")
n.iter <- 20000
n.burnin <- 1000
n.chains <- 1
n.thin <- 1

## Run Stan model
markov5.fit <- stan(file=model.file, iter=n.iter, warmup=n.burnin, chains=n.chains,
                    data=data, control=list(adapt_delta=0.995), seed=2019, 
                    thin=n.thin, pars=params)

## Cost effectiveness analysis
library(BCEA)
bc <- bcea(e = cbind("FP"=extract(markov5.fit)$stw_fp, 
                     "SFC"=extract(markov5.fit)$stw_sfc),
           c = cbind("FP"=extract(markov5.fit)$total_cost_fp, 
                     "SFC"=extract(markov5.fit)$total_cost_sfc),
           wtp = seq(0, 150, 5))
contour2(bc, wtp=50, graph="ggplot2")
ceac.plot(bc)

## Useful package to analyse Stan output, transforms output into tidy data frame, and produces pdf
## with ggplot2 plots
# library("ggmcmc")
# tidy <- ggmcmc::ggs(markov5.fit)
# ggmcmc::ggmcmc(tidy)
