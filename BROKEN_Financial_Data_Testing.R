## validation testing for GERGM function

# source('~/Dropbox/Cramner_ERGM_group/Denny_GERGM_Testing/Financial_Data_Testing.R')
#scl enable devtoolset-1.1 bash
#set wd 
rm(list = ls())
setwd("~/Dropbox/Cramner_ERGM_group/Denny_GERGM_Testing")
#setwd("J:/DB/Dropbox/Cramner_ERGM_group/Denny_GERGM_Testing")
# modifications to run on windows desktop
#setwd("J:/DB/Dropbox/Cramner_ERGM_group/Denny_GERGM_Testing")
PKG_CPPFLAGS = "-std=c++0x"
Sys.setenv(PKG_CPPFLAGS = PKG_CPPFLAGS)
# source script files 
library(RcppArmadillo)
#Rcpp::sourceCpp('MH_Sampler.cpp')
#Rcpp::sourceCpp('MH_Sampler_Together.cpp')
Rcpp::sourceCpp('MH_Sampler_Normal_Together.cpp')
source('GERGM.R', echo=F)
Report <- function(x){print(x)}

## Read in data.

## Country to Country Lending Data ##
 load("./Data/Transformed_Country_Lending_Data_1980-2005.Rdata")
# # transformed_adjacency_matrix_list
# # pareto_adjacency_matrix_list
# # log_adjacency_matrix_list

UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)
UMASS_RED <- rgb(153,0,51,255,maxColorValue = 255)
UMASS_GREEN <- rgb(0,102,102,255,maxColorValue = 255)
UMASS_YELLOW <- rgb(255,255,102,255,maxColorValue = 255)
UMASS_ORANGE <- rgb(255,204,51,255,maxColorValue = 255)

net <- log_adjacency_matrix_list[[1]]

transform.data <- matrix(1,18,18)
#formula.obj = net ~ in2star + out2star + ctriads + ttriads + recip
formula.obj = net ~  ttriads +recip

cat("Formula = ",toString(formula.obj),"\n")


## 3) gergm
shape = 0.1
sample_every = 2000
seed = 123
CUSTOM = F
MYTHETAS = list(par = c(0,0,0))
MH.fit <- gergm(formula.obj, directed = TRUE, seed = seed,
                transform.data = NULL,
                method = "Metropolis", max.num.iterations = 10, 
                mc.num.iterations = 100, nsim = 2000000, 
                MCMC.burnin = 200000, tolerance = 0.05, shape.parameter = 0.1, 
                together = 1, thin = 1/sample_every, weights = c(0.1,0.1),gain.factor = 0)

print(MH.fit)

MH.sims <- simulate.gergm(MH.fit, 2000000, seed = seed, method = "Metropolis")

save(list = c("MH.sims", "MH.fit"), file = "MHFitrecip.Rdata")

stats <- MH.fit@stats[1,]
par(mfrow = c(1,1))
boxplot(MH.sims$Statistics, medcol=UMASS_RED)
boxplot(rbind(stats,stats),add =T, medcol=UMASS_BLUE, names = F)

