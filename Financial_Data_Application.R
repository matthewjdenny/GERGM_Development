## Trials for Financial data

# source('~/Desktop/GERGM_Development/Financial_Data_Application.R')
#set year

YEAR = 1
#use this to open a shell in CentOS 6 that can be used to run the scripts
#scl enable devtoolset-1.1 bash

PKG_CPPFLAGS = "-std=c++0x"
Sys.setenv(PKG_CPPFLAGS = PKG_CPPFLAGS)

## This is working directly out of the GERGM folder in Github
setwd("~/Desktop/GERGM_Development")
library(RcppArmadillo)
Rcpp::sourceCpp('Scripts/MH_Sampler_Normal_Together.cpp')
source('Scripts/GERGM.R', echo=F)
Report <- function(x){print(x)}

load("Data/Transformed_Country_Lending_Data_1980-2005.Rdata")

RUN_MODEL = T
EVALUATE_OUTPUT = F


if(RUN_MODEL){
  # First model: transform.data only takes the mean of the realized values
  # The ERGM part of the model is a transitive.triads + reciprocity model
  net <- log_adjacency_matrix_list[[YEAR]]
  
  formula.obj = net ~  in2star + out2star+  ttriads
  
  shape = 0.05
  sample_every = 2000
  seed = 123
  CUSTOM = F
  
  MH.fit <- gergm(formula.obj, directed = TRUE, seed = seed,
                  transform.data = NULL,
                  method = "Metropolis", max.num.iterations = 10, 
                  mc.num.iterations = 100, nsim = 2000000, 
                  MCMC.burnin = 200000, tolerance = 0.001, shape.parameter = shape, 
                  together = 1, thin = 1/sample_every, weights = c(0.03,0.03,0.03), gain.factor = 0.5)
  
  #Simulate MH to see if we have any goodness of fit?
  MH.sims <- simulate.gergm(MH.fit, 2000000, seed = seed, method = "Metropolis", MCMC.burnin = 10000, thin = 1/sample_every, together = 1)
  
  
  True.stats <- MH.fit@stats[1,]
  temp.stats <- MH.sims$Statistics
  #Re weighting the recip and ttriads to correct level
  #temp.stats$recip = temp.stats$recip ^ (1/0.03)
  temp.stats$ttriads = temp.stats$ttriads ^ (1/0.03)
  temp.stats$in2stars = temp.stats$in2stars ^ (1/0.03)
  temp.stats$out2stars = temp.stats$out2stars ^ (1/0.03)

  ## plot goodness of fit
  indx = order(True.stats)
  
  
  #setwd("~/Desktop/GERGM_Development/Output")
  pdf(file = paste("MH_GOF_",YEAR,".pdf",sep = ""), width = 7, height = 5)
  par(oma=c(1,1,1,1))
  par(mar=c(4,4.5,2,1))
  spacing = seq(1,401, by = 15)
  par(mfrow = c(1,1))
  boxplot(log(temp.stats)[indx], ylab = "Log value", cex.axis = 1.5, cex.lab = 1.5, main = paste("Metropolis-Hastings", YEAR), cex.main = 1.5)
  lines(log(True.stats)[indx], type = "b", lwd = 3, pch = 5, lty = 8, col = "blue")
  legend("bottomright", "Observed", lwd = 3, pch = 5, lty = 8, col = "blue", cex = 1.25)
  dev.off()
  save(MH.fit, MH.sims, file = paste("Financial_Data_RTIO_",YEAR,".RData",sep = ""))
  #setwd("~/Desktop/GERGM_Development")
}

if(EVALUATE_OUTPUT){
  #### take a look at our saved output
  setwd("~/Desktop/GERGM_Development")
  load("~/Desktop/GERGM_Development/Financial_model_1.RData")
  
  #look at the MH parameter estimates
  MH.fit
  
  temp <- MH.fit@theta.coef
  z.stats <- as.numeric(temp[1,]/temp[2,])
  2*pnorm(-abs(z.stats[1]))
  ## 0.01065612
  2*pnorm(-abs(z.stats[2]))
  ## 0.214362
}






