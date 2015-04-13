## Trials for Financial data
## This is working directly out of the GERGM folder in Github
library(RcppArmadillo)
Rcpp::sourceCpp('Scripts/MH_Sampler_Normal_Together.cpp')
source('Scripts/GERGM.R', echo=F)
Report <- function(x){print(x)}

load("Data/Transformed_Country_Lending_Data_1980-2005.Rdata")

# First model: transform.data only takes the mean of the realized values
# The ERGM part of the model is a transitive.triads + reciprocity model
net <- log_adjacency_matrix_list[[1]]

formula.obj = net ~  recip + ttriads

shape = 0.05
sample_every = 2000
seed = 123
CUSTOM = F

MH.fit <- gergm(formula.obj, directed = TRUE, seed = seed,
                transform.data = NULL,
                method = "Metropolis", max.num.iterations = 10, 
                mc.num.iterations = 100, nsim = 2000000, 
                MCMC.burnin = 200000, tolerance = 0.001, shape.parameter = shape, 
                together = 1, thin = 1/sample_every, weights = c(0.1, 0.1), gain.factor = 0.10)

MPLE.fit <- gergm(formula.obj, directed = TRUE, seed = seed,
                  transform.data = NULL, MPLE.only = TRUE,
                  method = "Metropolis", max.num.iterations = 10, 
                  mc.num.iterations = 100, nsim = 2000000, 
                  MCMC.burnin = 200000, tolerance = 0.001, shape.parameter = shape, 
                  together = 1, thin = 1/sample_every, weights = c(0.1, 0.1), gain.factor = 0.10)

#Simulate MH to see if we have any goodness of fit?
MH.sims <- simulate.gergm(MH.fit, 2000000, seed = seed, method = "Metropolis", MCMC.burnin = 10000, thin = 1/sample_every, together = 1)
MPLE.sims <- simulate.gergm(MPLE.fit, 2000000, seed = seed, method = "Metropolis", MCMC.burnin = 10000, thin = 1/sample_every, together = 1)

True.stats <- MH.fit@stats[1,]
temp.stats <- MH.sims$Statistics
temp2.stats <- MPLE.sims$Statistics
#Re weighting the recip and ttriads to correct level
temp.stats$recip = temp.stats$recip ^ (1/0.1)
temp.stats$ttriads = temp.stats$ttriads ^ (1/0.1)

temp2.stats$recip = temp2.stats$recip ^ (1/0.1)
temp2.stats$ttriads = temp2.stats$ttriads ^ (1/0.1)

## plot goodness of fit
indx = order(True.stats)
par(oma=c(1,1,1,1))
par(mar=c(4,4.5,2,1))
spacing = seq(1,401, by = 15)
par(mfrow = c(1,1))
boxplot(log(temp.stats)[indx], ylab = "Log value", cex.axis = 1.5, cex.lab = 1.5, main = "Metropolis-Hastings", cex.main = 1.5)
lines(log(True.stats)[indx], type = "b", lwd = 3, pch = 5, lty = 8, col = "blue")
legend("bottomright", "Observed", lwd = 3, pch = 5, lty = 8, col = "blue", cex = 1.25)

par(mfrow = c(1,1))
boxplot(log(temp2.stats)[indx], ylab = "Log value", cex.axis = 1.5, cex.lab = 1.5, main = "MPLE", cex.main = 1.5)
lines(log(True.stats)[indx], type = "b", lwd = 3, pch = 5, lty = 8, col = "blue")
#legend("bottomright", "Observed", lwd = 3, pch = 5, lty = 8, col = "blue", cex = 1.25)

## plot MCMC trace plot on edge density
plot(unique(temp.stats[,6]), type = "l")

save(MH.fit, MH.sims, MPLE.fit, MPLE.sims, file = "Financial_model_1.RData")
