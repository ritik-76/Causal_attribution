source("MCMC_out.R")

gamma.chain.mat.1 <- matrix(NA, nrow = 12000, ncol = 35)

for ( i in 1:nrow(gamma.chain.mat.1)){
  gamma.chain.mat.1[i,] <- out.mcmc.f$gamma.chain[[i]]
}

library(coda)

# Geweke.stat Diagnostics

geweke.diag(mcmc(gamma.chain.mat.1))

# Effective Sample Size
ess(mcmc(gamma.chain.mat.1))