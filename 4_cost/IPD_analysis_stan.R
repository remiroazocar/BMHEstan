## SETUP
rm(list=ls()) 
# setwd("4_cost")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load cost data
cost.data=source("cost-data.txt")$value

filein <- "normal-mod.stan"
# Defines list of parameters to be tracked
params <- c("mu","ss","ls","delta_c","dev","total_dev")	

cost.fit <- stan(file=filein, iter=2000, warmup=1000, chains=2,
                 data=cost.data, control=list(adapt_delta=0.99),
                 seed=2019, thin=1, pars=params)

# print and manipulate output (stan simulation draws stored in cost.fit)
print(cost.fit, digits=3)
mu.draws <- extract(cost.fit)$mu 

plot(
  mu.draws[,1],									# What to plot on the x-axis?
  mu.draws[,2],									# What to plot on the y-axis?
  xlab="Population average cost in arm 1 (x £1000)",	# Label for the x-axis
  ylab="Population average cost in arm 2 (x £1000)",	# Label for the y-axis
  pch=20,												# Type of symbol used in the plot 
  # (see http://www.endmemo.com/program/R/pchsymbols.php)
  cex=.8,												# Amount by which plotting text and symbols should be scaled 
  # relative to the default (here 80%)
  main="Joint distribution of mean costs"				# Title for the plot
)

# Rescales the difference in cost by £1000
delta.c=extract(cost.fit)$delta_c*1000
hist(
  delta.c,							# variable to use for the histogram 
  xlab="Cost differential (£)",		# Label for the x-axis
  main="Posterior distribution"		# Title of the graph
)
# Adds a vertical line at 0
abline(
  v=0,		# Plots a vertical line at 0
  lwd=2,		# line width (=2pt)
  lty=2		# line type (=2 = a line)
)

# Computes the posterior probability that delta.c is positive (increase in average costs) as the 
# number of times that delta.c is positive divided by the total number of simulations 
nsims <- length(extract(cost.fit)$delta_c)
sum(delta.c>0)/nsims

# Compute the posterior means for the costs
mean(mu.draws[,1])
mean(mu.draws[,2])
mean(delta.c)

## Exercise 2.
# The deviance was calculated in the Stan model and saved as total_dev
mean(extract(cost.fit)$total_dev)
# Unlike for BUGS, the DIC is not implemented within Stan

## Exercise 4.
# Load data for costs & utilities into R from the txt file
cost.utility=source("cost-util-data.txt")$value

# Runs the Stan model from R
model.file="cgeg-mod.stan"	 # Specifies the new file with the cost-effectiveness model

# Here we load initial values for the 3 chains in the model
inits1=source("cgeg-inits1.txt")$value			# Loads the initial values for the 1st chain
inits2=source("cgeg-inits2.txt")$value			# Loads the initial values for the 2nd chain
inits3=source("cgeg-inits3.txt")$value			# Loads the initial values for the 3rd chain
# match names in Stan (underscores instead of periods)
names(inits1) <- names(inits2) <- names(inits3) <- c("mu_c1", "mu_c2", "mu_e1",
                                                     "mu_e2", "beta_c1", "beta_c2",
                                                     "shape_c1", "shape_c2", "shape_e1",
                                                     "shape_e2")

# Initial values must be provided as a list with each name matching with a variable name in the model
inits1
# inits must then be a list of lists where each list element is the initial values for EACH chain
inits=list(inits1,inits2,inits3) # Combines them into a single list --- can be lazy and only use 2...
params=c("mu_c","mu_e","delta_c","delta_e","INB","CEAC")

# Finally runs the model by calling Stan in the background
n.burnin=1000
n.iter=4000			# this adds 3000 simulations to the 1000 of burnin
n.chains=3

cost.fit.2 <- stan(file=model.file, iter=n.iter, warmup=n.burnin, chains=n.chains,
                   data=cost.utility, control=list(adapt_delta=0.995),
                   seed=2019, thin=1, init=inits, pars=params)

# Now can print & manipulate the output (stored in the object 'cost.fit.2')
print(cost.fit.2, digits=3)

# Plots the cost-effectiveness plane
plot(extract(cost.fit.2)$delta_e, extract(cost.fit.2)$delta_c, pch=20, cex=.8,
     xlab="Effectiveness differential",
     ylab="Cost differential", main="Cost-effectiveness plane")

# Calculates the summary measures for the cost-effectiveness analysis
mean(extract(cost.fit.2)$delta_e)
mean(extract(cost.fit.2)$delta_c)

# Calculates credible intervals for delta_e and delta_c
# Can change the level to compare with the slides using prob=
quantile(extract(cost.fit.2)$delta_e,prob=c(0.025,0.975))
quantile(extract(cost.fit.2)$delta_c,prob=c(0.025,0.975))

## Exercise 5.
# Figures out which of the indices is associated with willingness to pay = £500
K=numeric()				# defines an empty numeric vector, named 'K'
K.space=0.1				# in line with the Stan model
for (j in 1:11) {		# Loop to compute the K vector for all j values between 1 and 11
  K[j]=(j-1)*K.space
}
idx=which(K==0.5)		# Now finds the index for which K is equal to 0.5 
# NB: When dealing with logical operations (eg if, while, absolute equality), R wants the "double ="

# Now can compute the statistics for INB at K=£500
mean(extract(cost.fit.2)$INB[,idx])	# NB: the object INB is a matrix with as many rows
# as you've saved simulations and 11 columns so you have to pick all the rows in column 'idx'
# 2.5% quantile of the posterior distribution = lower end of the 95% interval
quantile(extract(cost.fit.2)$INB[,idx],.025)	
# 97.5% quantile of the posterior distribution = upper end of the 95% interval
quantile(extract(cost.fit.2)$INB[,idx],.975)	

# Now can identify the probability of cost-effectiveness for the new treatment (arm 2)
# First a graphical representation --- "visual inspection"
hist(extract(cost.fit.2)$INB[,idx],xlab="Incremental Net Benefit for k=£500",
     main="Posterior distribution")
abline(v=0,lwd=2,lty=2)
# Then can actually compute the probability as the proportion of simulations that are positive
n.sims <- length(extract(cost.fit.2)$INB)
prob.ce=sum(extract(cost.fit.2)$INB[,idx]>0)/n.sims
# Notice the commands used to calculate this probability
# n.sims is the total number of simulations  taken from the Bayesian model
# extract(cost.fit.2)$INB[,idx]>0 gives TRUE if greater than 0 and FALSE if not
# Numerically this is encoded as 1 for TRUE and 0 for FALSE
# Therefore, sum(...) adds up a 1 for each TRUE, giving the number of positive INB values.

# Finally plot the CEAC
# First extracts the mean values of the CEAC for each willingness to pay
CEAC=numeric()
for (i in 1:length(K)) { # length(K) is the length of the vector K = number of willingness to pay points
  CEAC[i]=mean(extract(cost.fit.2)$CEAC[,i])	# takes the mean of each column for the matrix CEAC
}
plot(K,CEAC,xlab="Willingness to pay (x £1000)",ylab="Cost-effectiveness acceptability curve", main="",
     t="l"	# This instructs R to plot using a solid curve (t="type", "l"=line)
)
