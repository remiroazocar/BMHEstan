# Bayesian methods in health economics in Stan

This repository contains implementations of the [BMHE Summer School](http://www.statistica.it/gianluca/teaching/summer-school/) practicals in Stan. 
Stan is a relatively new "probabilistic programming" language implementing MCMC algorithms for Bayesian inference. This repository will be of interest 
to the more avid students who want to conduct their analyses in Stan. 

## Stan

While BUGS/JAGS are currently the more popular languages for Bayesian analysis in health economics, Stan is a good alternative. Firstly, Stan uses the 
latest MCMC methodology (Hamiltonian Monte Carlo and the no-U-turn sampler). Secondly, it provides informative error messages, is very well documented
and has an active user base and development group. Finally, it can be integrated with R via RStan to produce great graphics. Notwithstanding, Stan does have its cons.
While it provides faster convergence (lesser iterations) for more complex models, sampling can be slow for simpler models. In addition, it has a steeper learning curve than
BUGS/JAGS and does not support missing data. 

Stan can called from within R (see instructions [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)). Stan code is more structured
than BUGS code and is set out in program blocks, each of which has a specific purpose:

* `data`: defines the data read into the model 
* `parameters`: defines the parameters sampled in the MCMC chains
* `transformed parameters`: defines parameter transformations
* `model`: defines model likelihoods
* `generated quantities`: defines functions of the parameters and the data

In R, we use the Rstan function `stan` to draw posterior samples. This function wraps the following three steps: 1.) translation of
the model in Stan code to C++ code, 2.) compilation of the C++ code, 3.) sampling. Stan takes some time to compile the model, which 
increases overall run time, but when the same model is fit again (possibly with new data/setting), recompilation can be avoided (as
the compiled object is saved). 

## Topics

* Introduction to health economic evaluations (3)
* Introduction to Bayesian inference (1)
* Introduction to Markov Chain Monte Carlo in Stan (2)
* Cost-effectiveness analysis with individual-level data (4)
* Aggregated-level data and hierarchical models (5)
* Evidence synthesis and network meta-analysis (6)
* Markov models (8)
