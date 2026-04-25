############################ ETA UPDATE #####################################

update_eta <- function(n, p, gamma, eta.hat, Sigma.star, Adj_mat, Sigma, X){
  
  # Inverse( Sigma )
  Sigma.inv <- chol2inv(chol(Sigma))
  # (D_W - W) kroneckar prod Sigma_inv
  Sigma.eta.1 <- kronecker.spam( X = as.spam(Adj_mat), Y = Sigma.inv)
  # Sigma.star + (D_W - W) kroneckar prod Sigma_inv
  Sigma.eta.inv <- Sigma.star + Sigma.eta.1
  # prod( Sigma.star, eta.hat)
  mu.eta.1 <- Sigma.star %*% eta.hat
  mu.eta.3 <- Sigma.eta.1 %*% (X %*% gamma)
  mu.eta.2 <- mu.eta.1 + mu.eta.3
  #mu.eta <- (Sigma.eta^-1) %*% mu.eta.2
  mu.eta <- solve.spam(Sigma.eta.inv, mu.eta.2)
  # mu.eta + chol(Sigma.eta) %*% e,  where e ~ rnorm(1750)
  #return( mu.eta + t(chol(Sigma.eta)) %*% rnorm(n*p) )
  return(as.numeric(mu.eta + backsolve(chol(Sigma.eta.inv),rnorm(n*p))))
}


################################ GAMMA  UPDATE ###############################################
update_gamma <- function(eta, Sigma, X){
  len.gamma <- ncol(X)
  # Inverse( Sigma )
  Sigma.inv <- chol2inv(chol(Sigma))
  # (D_W - W) kroneckar prod Sigma_inv
  Sigma.eta <- kronecker.spam( X = Adj_mat, Y = Sigma.inv)
  
  Sigma.gamma.inv <- crossprod.spam(X, Sigma.eta %*% X) + 1/(100)^2 * diag(len.gamma)
  
  Sigma.gamma <- chol2inv(chol(Sigma.gamma.inv))
  
  mu.gamma.1 <- crossprod.spam(X, Sigma.eta %*% eta)
  
  mu.gamma <- Sigma.gamma %*% mu.gamma.1
  
  return(mu.gamma + t(chol(Sigma.gamma)) %*% rnorm(len.gamma))
}


############################  SIGMA UPDATE ####################################


update_Sigma <- function(n, p, eta, gamma, Adj_mat, nu.Sigma, Si.Sigma, X) {
  
  # v' = v+n-1
  nu.post.Sigma <- nu.Sigma+n-1
  # Si'= matrix(eta - X %*% gamma) %*% (D_W-W) %*% t(matrix (eta - X %*% gamma))
  eta.1 <- eta - X %*% gamma
  
  matrix.eta <- matrix(data = eta.1, nrow = p, ncol = n )
  
  
  #Si.post.Sigma.1 <- crossprod.spam(as.spam(Adj_mat), t(matrix.eta))
  
  #eta_mat <- t(matrix(eta, 7,250))
  
  Si.post.Sigma <- matrix.eta %*% Adj_mat %*% t(matrix.eta) + Si.Sigma
  
  return( rinvwishart( nu = nu.post.Sigma, S = Si.post.Sigma))
}

########################### FINAL MCMC UPADATE ################################

mcmc_update_causal <- function( eta.init, eta.hat, gamma.init, Sigma.init, nu.Sigma, Si.Sigma, X, n_iter = 60000,burn = 10000,thin = 5){
  
  n <- length(eta.init) / nrow(Sigma.init) #nsites
  p <- nrow(Sigma.init) # dim(parameter vector at each site)
  
  eta.chain <- list()
  Sigma.chain <- list()
  gamma.chain <- list()
  
  eta <- eta.init
  Sigma <- Sigma.init
  gamma <- gamma.init
  
  return.iter <- (burn+1):n_iter
  
  
  for ( i in 1:n_iter){
    for (t in 1:thin){
      eta <- update_eta(n, p, gamma, eta.hat, Sigma.star, Adj_mat, Sigma, X)
      gamma <- update_gamma(eta, Sigma, X)
      Sigma <- update_Sigma(n, p, eta, gamma, Adj_mat, nu.Sigma, Si.Sigma, X)
    }
    
    
    
    eta.chain[[i]] <- eta
    Sigma.chain[[i]] <- Sigma
    gamma.chain[[i]] <- gamma 
    
    cat("\t iter", i, "\n")
  }
  
  out <- list(eta.chain=eta.chain, gamma.chain=gamma.chain, Sigma.chain=Sigma.chain )
  return(out)
}