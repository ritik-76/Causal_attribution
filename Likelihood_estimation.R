
library(extRemes)
library(ismev)
library(ggplot2)
library(fields)
library(spam)
library(Matrix)
library(LaplacesDemon)
library(parallel)
library(doParallel)

load("annual_max_temp_USA_nat.RData")
load("annual_max_temp_USA_hist.RData")
t<-c(1850:2014)

################################ HR LIKELIHOOD #################################
loglik_hr <- function(par1, y1, y2) {
  # Parameter
  r <- exp(par1)
  eps <- 1e-12
  
  y1 <- pmax(y1,eps)
  y2 <- pmax(y2,eps)
  # Auxiliary terms
  a1 <- 1/r + 0.5 * r * log(y1 / y2)
  a2 <- 1/r + 0.5 * r * log(y2 / y1)
  phi1 <- dnorm(a1)
  phi2 <- dnorm(a2)
  Phi1 <- pnorm(a1)
  Phi2 <- pnorm(a2)
  
  
  # Exponent measure and derivatives
  V <- y1 * Phi1 + y2 * Phi2
  V1 <- Phi1 + (r/2)*phi1 - (r/2)*(y2/y1)*phi2
  V2 <- Phi2 + (r/2)*phi2 - (r/2)*(y1/y2)*phi1
  V12 <- -(r/2)*(phi1/y2 + phi2/y1) + (r^2/4)*(a1*phi1/y2 - a2*phi2/y1)
  
  # Avoid log of negative or zero arguments
  val <- V1 * V2 - V12
  val <- pmax(val, eps)
  
  
  # Log-likelihood contributions
  
  
  
  logf <- -V + log(val)
  return(-sum(logf))
  
}

############################### LIKELIHOOD ESTIMATION #############################

ncores <- detectCores() - 3
cl <- makeCluster(ncores)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1234)
B <- 10000
t <- c(1850:2014)
info.mat_list <- list()
eta_hat_mat <- matrix(NA, nrow = 7, ncol = 250)
theta <- matrix(data = NA, nrow = B, ncol = 7)


clusterExport(cl, c("annual_max_temp_1", "annual_max_temp_2","B","t", "loglik_hr",
                    "optim", "info.mat_list", "eta_hat_mat", "theta"))
clusterEvalQ(cl, library(evd))

bootstrap_par <- function(i){
  z1 <- as.numeric(annual_max_temp_1[i, paste0("annual", 1:165)])
  z2 <- as.numeric(annual_max_temp_2[i, paste0("annual", 1:165)])
  theta <- matrix(data = NA, nrow = B, ncol = 7)
  fit <- fgev(z1,nsloc = (t-1800)/50,
              method= "BFGS",
              control = list(maxit = 5000),
              std.err = FALSE)
  fit1 <- fgev(z2,nsloc = (t-1800)/50,
               method= "BFGS",
               control = list(maxit = 5000),
               std.err = FALSE)
  
  #print(as.numeric(fit$estimate))
  #print(as.numeric(fit1$estimate))
  
  mu1 <- fit$estimate[1] + (t-1800)/50 * fit$estimate[2]
  sigma1 <- fit$estimate[3]
  shape1 <- fit$estimate[4]
  k1 <- 1 + shape1 * (z1-mu1)/sigma1
  y1 <- k1 ^ (-1/shape1)
  
  
  mu2 <- fit1$estimate[1] + (t-1800)/50 * fit1$estimate[2]
  sigma2 <- fit1$estimate[3]
  shape2 <- fit1$estimate[4]
  k2 <- 1 + shape2 * (z2-mu2)/sigma2
  y2 <- k2 ^ (-1/shape2)
  
  
  
  logjc <- - log(sigma1) - log(sigma2) - ((1/shape1)+1)* log (k1) - ((1/shape2)+1)* log (k2)
  
  negbvlik_hr <- function(par1, x1, x2){
    loglik_hr(par1,x1,x2) - sum(logjc)
  }
  
  
  start1 <- log(0.55)
  
  fit2 <- optim(par = start1, fn = negbvlik_hr,
                x1 = y1, x2 = y2, method = "BFGS",hessian = TRUE,
                control = list(maxit = 500))
  
  r_hat <- exp(fit2$par)
  #print(r_hat)
  x <- cbind(z1,z2)
  
  shape1.start <- max(min(max(fit$estimate[4],fit1$estimate[4]), 0.15), -0.45)
  scale1.start <- max(max(fit$estimate[3],fit1$estimate[3]), 0.05)
  dep.start <- max(min(r_hat, 5), 0.2)
  
  fit.rh <- fbvevd(x, model = 'hr',
                   start = list(loc1 = fit$estimate[1],loc1trend = fit$estimate[2],
                                scale1 = scale1.start, 
                                shape1 = shape1.start,
                                loc2 = fit1$estimate[1],loc2trend = fit1$estimate[2],dep = dep.start),
                   nsloc1 = (t-1800)/50, nsloc2 = (t-1800)/50, 
                   cshape = TRUE, cscale = TRUE,
                   method = "BFGS",
                   std.err = FALSE, warn.inf = FALSE)
  
  mu_1 <- as.numeric(fit.rh$estimate[1])+(t-1800)/50 * as.numeric(fit.rh$estimate[2])
  sigma_1 <- rep(as.numeric(fit.rh$estimate[3]),165)
  shape_1 <- rep(as.numeric(fit.rh$estimate[4]),165)
  margins1 <- cbind(mu_1,sigma_1,shape_1)
  
  #margins1
  
  mu_2 <- as.numeric(fit.rh$estimate[5])+(t-1800)/50 * as.numeric(fit.rh$estimate[6])
  
  margins2 <- cbind(mu_2,sigma_1,shape_1)
  
  #margins2
  for ( b in 1:B){
    
    
    f <- rbvevd(165, dep = fit.rh$estimate[7], model = "hr", mar1 = margins1, mar2 = margins2)
    z1 <- f[,1]
    z2 <- f[,2]
    #print(cbind(z1,z2))
    x <- cbind(z1,z2)
    
    fit.rh <- fbvevd(x, model = 'hr',
                     start = list(loc1 = fit$estimate[1],loc1trend = fit$estimate[2],
                                  scale1 = scale1.start, 
                                  shape1 = shape1.start,
                                  loc2 = fit1$estimate[1],loc2trend = fit1$estimate[2],dep = dep.start),
                     nsloc1 = (t-1800)/50, nsloc2 = (t-1800)/50, 
                     cshape = TRUE, cscale = TRUE,
                     method = "BFGS",
                     control = list(maxit = 5000),
                     std.err = FALSE, warn.inf = FALSE)
    
    theta[b,] <- as.numeric(fit.rh$estimate)
    
  }
  
  Cov_mat <- cov(theta)
  print(diag(Cov_mat))
  
  eta_hat_mat[,i] <- colMeans(theta)
  
  info.mat_list[[i]] <- chol2inv(chol(Cov_mat))
  
  return(list(eta_hat_mat_i = eta_hat_mat[,i], 
              info.mat_list_i = info.mat_list[[i]]))
}

output_list <- parLapply(cl, 1:250, fun = bootstrap_par)

stopCluster(cl)


eta_hat_mat_par <- sapply(output_list, `[[`, "eta_hat_mat_i")
info.mat_list_par <- lapply(output_list, `[[`, "info.mat_list_i")
eta_hat_mat_par <- eta_hat_mat_par[c(1,2,5,6,3,4,7), ]
#eta_hat_mat_par[5,]

eta_hat_mat <- eta_hat_mat_par
info.mat_list <- info.mat_list_par

A1 <- rbind(c(1, (mean(t)-1800)/50), c(0, sd(t)/50))
A <- as.matrix(bdiag(A1,A1,diag(3)))

for ( i in 1:ncol(eta_hat_mat)){
  eta_hat_mat[,i] <- A %*% eta_hat_mat[,i]
}

save(info.mat_list, file = "Information_Matrix.RData")
save(eta_hat_mat, file = "Estimated_eta_bootstrap.RData")