# Setup
rm(list=ls()) 

# Load packages
library("rstan")
library("coda")
library("BCEA")

# setwd("6_evsynth2")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

### Smoking cessation network meta-analysis data in format obtained
### from Lu & Ades tutorial "Introduction to Mixed Treatment
### Comparisons"

smoke <- read.table("smoke_data_orig.txt", header=TRUE, nrow=24)
tnames <- c("A: None","B: Self-help","C: Individual","D: Group")
names(smoke) <- gsub("\\.", "", names(smoke))
ns <- nrow(smoke)
nt <- 4
r <- smoke[,c("r1","r2","r3")]
n <- smoke[,c("n1","n2","n3")]
t <- smoke[,c("t1","t2","t3")]
## Order by comparison
comp <- paste(t[,1], t[,2], ifelse(is.na(t[,3]),"",t[,3]), sep="")
r <- as.matrix(r)[order(comp),]
n <- as.matrix(n)[order(comp),]
t <- as.matrix(t)[order(comp),]
na <- smoke$na[order(comp)]
n[is.na(r)] <- NA

### Format data for table in the slides
rc <- nc <- matrix(nrow=ns, ncol=nt)
for (i in 1:2) {   # ith col should be treatment i
  rc[cbind(1:ns, t[,i])] <- r[,i]
  nc[cbind(1:ns, t[,i])] <- n[,i]
}
threes <- which(na==3)
rc[cbind(threes, t[threes,3])] <- r[threes,3]
nc[cbind(threes, t[threes,3])] <- n[threes,3]
dc <- matrix(paste(rc, nc, sep="/"), nrow=ns, ncol=nt)
dc[dc=="NA/NA"] <- ""
write.table(dc[order(comp),], sep=" & ", quote=FALSE, eol="\\\\\n")

# Prepare data for Stan: Stan does not support NAs so replace by dummy value e.g. -1
# (we do not loop over NAs in Stan model so these do not affect inference)

r[is.na(r)] <- -1
n[is.na(n)] <- -1
t[is.na(t)] <- -1
smoke.list <- list(r=r, n=n, t=t, na=na, ns=ns, nt=nt)

###############    FIT MODELS IN Stan via RStan    ###############

### Initial values
inits <- list(list(mu=rep(0,24), d_other=c(0,0,0)),
              list(mu=rep(-1,24), d_other=c(1,1,1)))

### FIXED EFFECTS MODEL

### Run Stan model 
res <- stan(file="smokefix_model.stan", data=smoke.list, init=inits,
            pars=c("d", "o_r", "L", "pq"), chains=2, warmup=1000, iter=20000, 
            seed=2019, thin=1, control=list(adapt_delta=0.995)) 
print(res, digits=3)

# Evidence of convergence:
# i) trace plot ("hairy caterpillars" indicate convergence)
rstan::traceplot(res, pars=c("d[2]","d[3]","d[4]"))
# ii) Rhat < 1.1 (Brooks, Gelman and Rubin diagnostic) is more evidence of adequate convergence
library("bayesplot")
rhats <- rhat(res, pars=c("d[2]","d[3]","d[4]"))
mcmc_rhat(rhats) 
print(summary(res, pars=c("d[2]","d[3]","d[4]"))$summary) # see Rhat and n_eff
# iii) How many further iterations do we need after convergence?  Consider
# the "effective sample size".  Rule of thumb is you need n.eff>4000
# to ensure 95\% credible limits have 94.5-95.5% coverage (true in
# this case).  Or look at Monte Carlo error, or just check results
# don't change when you run more.
ratios <- neff_ratio(res, pars=c("d[2]","d[3]","d[4]"))
mcmc_neff(ratios)

st <- summary(res, pars=c("o_r"))$summary  # Gives mean and 95% CIs for odds ratios, for plotting

## Organise results for plotting
or <- st[grep("o_r", rownames(st)),c("2.5%", "50%", "97.5%")]
colnames(or) <- c("l95","est","u95")
or <- as.data.frame(or)
or$act <- as.numeric(sub("o_r\\[([0-9]+),([0-9]+)\\]", "\\1", rownames(or)))
or$com <- as.numeric(sub("o_r\\[([0-9]+),([0-9]+)\\]", "\\2", rownames(or)))
or <- or[or$act>or$com,]
or <- or[order(or$com,or$act),]

### RANDOM EFFECTS MODEL

smoke.list <- list(r=r, n=n, t=t, max_na=max(na),na=na, ns=ns, nt=nt)
inits <- list(list(mu=rep(0,24), d_other=c(0,0,0), stdev=1),
              list(mu=rep(-1,24), d_other=c(1,1,1), stdev=2))

# Note: the current implementation of this model is slow (takes 5 minutes to run). 
# The maximum allowed tree depth has been increased from 10 to 10.5 to avoid hitting
# the maximum tree depth. Hitting the max. tree depth is an efficiency concern:
# if it is reached, the algorithm terminates prematurely to avoid excessively long
# execution time
res.re <- stan(file="smokere_model.stan", data=smoke.list, init=inits,
               chains=2, warmup=1000, iter=2000, pars=c("o_r", "d", "stdev", "pq", "L"),
               seed=2019, thin=1, control=list(adapt_delta=0.99,
                                               max_treedepth=10.5))

# Evidence of convergence: check convergence of random effects SD in particular!
# i) trace plot ("hairy caterpillars" indicate convergence)
rstan::traceplot(res.re, pars=c("d[2]","d[3]","d[4]", "stdev"))
# ii) Rhat < 1.1 (Brooks, Gelman and Rubin diagnostic) is more evidence of adequate convergence
library("bayesplot")
rhats <- rhat(res.re, pars=c("d[2]","d[3]","d[4]","stdev"))
mcmc_rhat(rhats) 
print(summary(res.re, pars=c("d[2]","d[3]","d[4]","stdev"))$summary) # see Rhat and n_eff
# iii) effective sample size
ratios <- neff_ratio(res.re, pars=c("d[2]","d[3]","d[4]","stdev"))
mcmc_neff(ratios)

st.re <- summary(res.re, pars=c("o_r","stdev"))$summary  # Gives mean and 95% CIs for odds ratios, for plotting

## Organise results for plotting
or.re <- st.re[grep("o_r", rownames(st.re)),c("2.5%", "50%", "97.5%")]
colnames(or.re) <- c("l95","est","u95")
or.re <- as.data.frame(or.re)
or.re$act <- as.numeric(sub("o_r\\[([0-9]+),([0-9]+)\\]", "\\1", rownames(or.re)))
or.re$com <- as.numeric(sub("o_r\\[([0-9]+),([0-9]+)\\]", "\\2", rownames(or.re)))
or.re <- or.re[or.re$act>or.re$com,]
or.re <- or.re[order(or.re$com,or.re$act),]

## Calculate direct evidence - classical ORs and CIs.
## These exist for all comparisons here

fe <- dati <- ntrials <- numeric()
datij <- array(dim=c(4,4,4)); dimnames(datij) <- list(c("r1","r2","n1","n2"),NULL,NULL)
for (i in 1:3) { # comparator
  for (j in (i+1):4){ # active
    pw <- !is.na(rc[,i]) & !is.na(rc[,j])
    ntrials <- c(ntrials, sum(pw))
    r1 <- sum(rc[pw,i]); n1 <- sum(nc[pw,i])
    r2 <- sum(rc[pw,j]); n2 <- sum(nc[pw,j])
    lori <- log(r2) - log(n2-r2) - log(r1) + log(n1-r1)
    sei <- sqrt(1/r2 - 1/(n2-r2) + 1/r1 - 1/(n1-r1))
    fe <- rbind(fe, c(j, i, exp(lori+qnorm(c(0.025, 0.5, 0.975))*sei)))
    colnames(fe) <- c("act","com","l95","est","u95")
    datij[,i,j] <- c(r1,r2,n1,n2)
    dati <- rbind(dati, c(r1,r2,n1,n2))
  }
}
fe

## PLOT:  direct against fixed effects mixed

# par(lwd=2, mar=c(2.1, 0.1, 0.1, 2.1))
plot(or[,"est"], 1:6, type="n", xlim=c(-2, 5), ylim=c(0,6.5), ylab="", bty="n", xlab="Odds ratio", axes=FALSE)
points(or[,"est"], 1:6, pch=19, col="red")
axis(1, at=0:5)
abline(v=c(1), col="lightgray", lwd=1)
segments(or[,"l95"], 1:6, or[,"u95"], 1:6, col="red")
points(fe[,"est"], 1:6-0.2, pch=19)
segments(fe[,"l95"], 1:6-0.2, fe[,"u95"], 1:6-0.2)
text(-0.2,6.5,"Pooled N", pos=4)
text(-2, 1:6, paste(tnames[fe[,"act"]],tnames[fe[,"com"]],sep=" / "), pos=4, cex=0.8)
N <- dati[,3]+dati[,4]
text(-0.1, 1:6, N, pos=4)
legend("bottomright", col=c("red","black"), c("Fixed effects mixed","Fixed effects direct"), lwd=c(2,2,2), bty="n")

#dev.off()

## PLOT:  direct against fixed effects against random effects

plot(or[,"est"], 1:6, type="n", xlim=c(-2, 5), ylim=c(0,6.5), ylab="", bty="n", xlab="Odds ratio", axes=FALSE)
points(or[,"est"], 1:6, pch=19, col="red")
axis(1, at=0:5)
abline(v=c(1), col="lightgray", lwd=1)
segments(or[,"l95"], 1:6, or[,"u95"], 1:6, col="red")
points(fe[,"est"], 1:6-0.2, pch=19)
segments(fe[,"l95"], 1:6-0.2, fe[,"u95"], 1:6-0.2)
points(or.re[,"est"], 1:6+0.2, pch=19, col="blue")
segments(or.re[,"l95"], 1:6+0.2, or.re[,"u95"], 1:6+0.2, col="blue")
arrows(or.re[,"l95"], 3.2, 5.3, 3.2, col="blue", length=0.1)
text(-0.2,6.5,"Pooled N", pos=4)
text(-2, 1:6, paste(tnames[fe[,"act"]],tnames[fe[,"com"]],sep=" / "), pos=4, cex=0.8)
N <- dati[,3]+dati[,4]
text(-0.1, 1:6, N, pos=4)
legend("bottomright", col=c("blue","red","black"), c("Random effects mixed","Fixed effects mixed","Fixed effects direct"), lwd=c(2,2,2), bty="n")
text(x=or.re[,"est"], y=1:6+0.29,  col="blue", cex=0.7,
     labels=paste(round(or.re[,"est"],2), " (",round(or.re[,"l95"],2), ", ", round(or.re[,"u95"],2), ")", sep=""))   
text(x=2.2, y=6.29, col="blue", cex=0.7, labels="Posterior median (95% credible interval)", pos=4)

## Plot heterogeneity between trials in direct ORs
hres <- numeric()
for (i in 1:3) { # comparator
  for (j in (i+1):4){ # active
    pw <- !is.na(rc[,i]) & !is.na(rc[,j])
    r1 <- rc[pw,i]; n1 <- nc[pw,i]
    r2 <- rc[pw,j]; n2 <- nc[pw,j]
    lori <- log(r2) - log(n2-r2) - log(r1) + log(n1-r1)
    sei <- sqrt(1/r2 - 1/(n2-r2) + 1/r1 - 1/(n1-r1))
    hresi <- exp(cbind(lori, lori+qnorm(0.025)*sei, lori+qnorm(0.975)*sei))
    hres <- rbind(hres, cbind(j, i, hresi))
  }
}
colnames(hres)[3:5] <- c("or","l95","u95")
hres[c(7,20),"or"] <- -100
hres[c(7,20),"l95"] <- 1
hres[c(7,20),"u95"] <- 1

par(bg="white", xlog=TRUE)
comps <- paste(tnames[hres[,"j"]],tnames[hres[,"i"]],sep=" / ")
ypos <- 1:nrow(hres) + rep(1:6*2, table(comps))
plot(0, ylim=range(ypos), xlim=c(0.1, 10), axes=FALSE, bty="n", ylab="", xlab="Odds ratio (log scale axis)", type="n", log="x")
apos <- cumsum(table(comps)+2)+1
del <- 0.1
lrec <- 0.1; urec <- 10
rect(lrec, 0+del, urec, apos[1], col="lightgray", border=NA)
rect(lrec, apos[2]+del, urec, apos[3], col="lightgray", border=NA)
rect(lrec, apos[4]+del, urec, apos[5], col="lightgray", border=NA)
rect(lrec, apos[1]+del, urec, apos[2], col="lightblue", border=NA)
rect(lrec, apos[3]+del, urec, apos[4], col="lightblue", border=NA)
rect(lrec, apos[5]+del, urec, apos[6], col="lightblue", border=NA)
points(hres[,"or"], ypos, pch=19)
off.inds <- c(7,8,9,20)
hres[off.inds,"u95"] <- 10.5
segments(hres[,"l95"], ypos, hres[,"u95"], ypos)
rinf <- hres[off.inds,]
arrows(rinf[,"l95"], ypos[off.inds], rinf[,"u95"], ypos[off.inds], length=0.1)
axis(1, at=c(0.2, 0.5, 0:5, 10))
abline(v=1, col="purple")
amid <- c(0, apos[1:5]) + diff(c(0,apos))/2
text(0.1, amid, unique(comps), pos=4)

### Cost-effectiveness analysis (random effects model)
unit.cost <- c(0,200,6000,600)
ints <- c("No contact","Self help","Individual counselling","Group counselling")
e <- c <- matrix(NA,length(extract(res.re)$d[,1]),4)
L <- extract(res.re)$L # MCMC sample from distribution of life-years gained by quitting
pq <- extract(res.re)$pq # ...and from distributions of probability of quitting for each of 4 interventions

for (t in 1:4) {
  e[,t] <- L*pq[,t]
  c[,t] <- unit.cost[t]
}
colnames(e) <- colnames(c) <- ints
round(apply(e, 2, quantile, c(0.025, 0.5, 0.975)), 1) # results on slide
m <- bcea(e,c,interventions=ints)
summary(m)
plot(m)
