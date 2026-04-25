
source("MCMC_out.R")

eta.chain.mat.1 <- matrix(NA,nrow = 12000, ncol = 1750)
gamma.chain.mat.1 <- matrix(NA, nrow = 12000, ncol = 35)

for ( i in 1:nrow(eta.chain.mat.1)){
  eta.chain.mat.1[i,] <- out.mcmc.f$eta.chain[[i]]
  gamma.chain.mat.1[i,] <- out.mcmc.f$gamma.chain[[i]]
}
burn <- 2000
n_iter <- 12000
eta.chain.mat <- eta.chain.mat.1[-(1:burn),]
gamma.chain.mat <- gamma.chain.mat.1[-(1:burn),]

eta.mean.post <- colSums(eta.chain.mat)/10000
gamma.mean.post <- colSums(gamma.chain.mat)/1000
gamma.sd.post <- c()
for(i in 1:ncol(gamma.chain.mat)){
  gamma.sd.post <- append(gamma.sd.post, sd(gamma.chain.mat[,i]))
}
colname_coeff <- c("alpha_0","alpha_1","beta_0","beta_1","sigma_star", "xi_star","lambda")
coeff_df <- data.frame(matrix(ncol = length(colname_coeff), nrow = 5))
colnames(coeff_df) <- colname_coeff
# posterior means of all hyperparameters in \gamma
coeff_df <- matrix(gamma.mean.post, ncol = 7, nrow = 5)
# posterior sds of all hyperparameters in \gamma
coeff_df_sd <- matrix(gamma.sd.post, ncol = 7, nrow = 5)