#load("Design_Matrix.RData")
#load("Estimated_eta_bootstrap.RData")
#load("Adjacency_Matrix.RData")
#load("Sigma_star.RData")
source("Update_MCMC.R")
source("Likelihood_estimation.R")
source("Model_components.R")


eta.hat <- as.vector(eta_hat_mat)
n <- 250
p <- 7
nu.Sigma <- 0.1
Si.Sigma <- 0.1 * diag(p)
eta.init <- eta.hat
tXX <- crossprod.spam(X, X)
tXX.inv <- chol2inv.spam(chol.spam(as.spam(tXX)))
tXeta <- crossprod.spam(X, eta.init)
gamma.init <- tXX.inv %*% tXeta

eta.init.1 <- eta.init - X %*% gamma.init
matrix.eta.init <- matrix(data = eta.init.1, nrow = p, ncol = n )
Sigma.init <- matrix.eta.init %*% Adj_mat %*% t(matrix.eta.init)/n
out.mcmc <- mcmc_update_causal(eta.init, eta.hat, gamma.init, Sigma.init, nu.Sigma, Si.Sigma, X, n_iter = 12000,burn = 2000, thin = 5)
save(out.mcmc, file="final_output_MCMC.RData")