## R Script for Fitting Generalized Exponential Random Graph Models to Migration Data Set
## Written by James D. Wilson and Matthew Denny 2/10/15

# modifications to run on windows desktop
# setwd("J:/DB/Dropbox/Cramner_ERGM_group/Denny_GERGM_Testing")
# PKG_CPPFLAGS = "-std=c++0x"
# Sys.setenv(PKG_CPPFLAGS = PKG_CPPFLAGS)

#run on Cluster
#source('~/Dropbox/Cramner_ERGM_group/Denny_GERGM_Testing/GERGM_Testing/Migration_Reproducible_Script.R')

setwd("/home/Denny/Dropbox/Cramner_ERGM_group/Denny_GERGM_Testing")
PKG_CPPFLAGS = "-std=c++0x"
Sys.setenv(PKG_CPPFLAGS = PKG_CPPFLAGS)

## Source script files 

library(RcppArmadillo)
Rcpp::sourceCpp('MH_Sampler_Normal_Together.cpp')
source('GERGM.R', echo=F)
Report <- function(x){print(x)}

## Read in data

load("./Data/Migration_Data.Rdata")

# Notes: This .Rdata file contains the following data structures:
# data - a list variable of length 2 containing 
#        a) Graph - the differential graph between 2007 and 2006 U.S. Migration: 51 x 51 adjacency matrix
#        b) Z - transformation variables described in Wilson et al. 2015. This is a 51 x 51 x 10 array
# transform.data - Z (see (b) above)

## Construct the formula object for the gergm() function

net <- data$Graph

formula.obj = net ~ in2star + out2star + ctriads + ttriads + recip

## Set seed
seed = 123
CUSTOM = F

## Run Gibbs Model Fit
Gibbs.fit <- gergm(formula.obj, directed = TRUE, seed = seed,
                   transform.data = transform.data,
                   method = "Gibbs", max.num.iterations = 10, 
                   mc.num.iterations = 100, nsim = 1000, 
                   MCMC.burnin = 500, tolerance = 1e-04, shape.parameter = 1, 
                   together = 1, gain.factor = 0)

## Run MH Model Fit (April 5th, 2015)
MH.fit <- gergm(formula.obj, 
				directed = TRUE, 
				seed = seed,
                transform.data = transform.data,
                method = "Metropolis", 
				max.num.iterations = 10, 
                mc.num.iterations = 100, 
				nsim = 1000000, 
                MCMC.burnin = 10000, 
				tolerance = 1e-04, 
				shape.parameter = 0.01, 
                together = 1, 
				thin = 0.001, 
				gain.factor = 0
				)

## Run MPLE fit (April 8th, 2015) TODO: fix the theta comparison problem
MPLE.fit <- gergm(formula.obj, directed = TRUE, MPLE.only = TRUE, seed = seed,
                  max.num.iterations = 10, transform.data = transform.data,
                  mc.num.iterations = 100, tolerance = 1e-04, shape.parameter = 0.01, 
                  together = 1, thin = 0.001, gain.factor = 0)

## Save the results
save(Gibbs.fit, file = "Migration_Replication_Gibbsfit_Final.RData")
save(MH.fit, file = "Migration_Replication_MHfit_Final.RData")
save(MPLE.fit, file = "Migration_Replication_MPLEfit_Final.RData")


###Do this when I get the results back!
# ## Simulate 1000 networks and compare results with observed statistics
# Gibbs.sims <- simulate.gergm(MH.fit, 1000, seed = seed, method = "Gibbs", MCMC.burnin = 1000)
# 
# #We are keeping every 1000 networks in the MH to avoid autocorrelation
# MH.sims <- simulate.gergm(Gibbs.fit, 1000, seed = seed, method = "Metropolis", MCMC.burnin = 20000, thin = 0.01, shape.parameter = 0.01)
# ## Plot and compare network statistics
# #gof.plot(Gibbs.sims, MH.fit) 
# 
# ##Compare statistics (from MH parameter estimates )
# library(caroline)
# gof.compare.plot(Gibbs.sims, MH.sims, MH.fit)
# 
# save(MH.sims, Gibbs.sims, file = "MH_Gibbs_Simulations_from_MH_results.RData")


