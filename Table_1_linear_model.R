############################ REGRESSION COEFFICIENTS ##########################

source("Likelihood_estimation.R")


Design_Mat_list <- list()
locs <- as.matrix(annual_max_temp_1[,1:2])
locs <- as.data.frame(locs)
elevs <- read.csv("USA_1deg_elevation_stats.csv")
elevs <- as.data.frame(elevs)
mean.sea.dist <- as.data.frame(read.csv("USA_1deg_mean_sea_distance.csv"))

for ( i in 1:length(locs$lon)){
  Design_Mat_list[[i]] <- kronecker(diag(7), t(c(1,locs$lon[i],locs$lat[i], 
                                                 elevs$mean_elev[i],
                                                 mean.sea.dist$mean_sea_dist_m[i]/1000)))
  
}
X <- do.call(rbind, Design_Mat_list)


coeff.pred <- matrix(NA, nrow = 7, ncol = 6)
residual.mat <- matrix(NA, nrow=7, ncol=250)
for( k in 1:7){
  XX <- list()
  
  for ( i in 1:length(locs$lon)){
    XX[[i]] <-  t(c(eta_hat_mat[k,i],locs$lon[i],locs$lat[i],
                    elevs$mean_elev[i],elevs$sd_elev[i],
                    mean.sea.dist$mean_sea_dist_m[i]/1000))
  }
  covariate.fit <- do.call(rbind, XX)
  colnames(covariate.fit) <- c("GEV_coeff","X1","X2","X3","X4","X5")
  covariate.fit <- as.data.frame(covariate.fit)
  fit.res <- lm(GEV_coeff ~ X1+X2+X3+X4+X5, data = covariate.fit)
  print(summary(fit.res))
  covariate.fit <- as.matrix(covariate.fit)
  covariate.fit[,1] <- rep(1,250)
  residual.mat[k,] <- eta_hat_mat[k,] - 
    as.matrix(covariate.fit) %*%(as.numeric(coefficients(fit.res)))
  
}