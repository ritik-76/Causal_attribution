load("Estimated_eta_bootstrap.RData")
load("Information_Matrix.RData")

#source("Likelihood_estimation.R")
######################## TRANSFORMATION OF PARAMETER ###########################

t <- c(1850:2014)

A1 <- rbind(c(1, (mean(t)-1800)/50), c(0, sd(t)/50))
A <- as.matrix(bdiag(A1,A1,diag(3)))
A.inv <- solve(A)
for ( i in 1:ncol(eta_hat_mat)){
  eta_hat_mat[,i] <- A %*% eta_hat_mat[,i]
}
A.inv <- as.spam(A.inv)

eta_hat_mat[5,]<- log(eta_hat_mat[5,])
eta_hat_mat[6,]<- h1.si(eta_hat_mat[6,])
eta_hat_mat[7,]<-log(eta_hat_mat[7,])

########################### BLOCK DIAGONAL MATRIX ##############################
c.phi <- 0.005
b.phi <- 0.25 * (1/c.phi) * (1- (0.25)^ c.phi)
a.phi <- -0.25 - (b.phi * log( (0.25)^c.phi / (1- (0.25)^c.phi)))


g.phi.der <- function(phi){
  (1/(b.phi * c.phi)) * (g1.phi(phi)+0.5)^(c.phi + 1) * exp(-(phi - a.phi)/b.phi)
}
trans.inf.matrix <- list()
Sigma.star <- NULL
D <- diag(c(rep(1,4), exp(eta_hat_mat[5,1]), g.phi.der((eta_hat_mat[6,1])), 
            exp(eta_hat_mat[7,1])))
trans.inf.matrix[[1]] <- crossprod.spam(D, (info.mat_list[[1]]) %*% D )
trans.inf.matrix[[1]] <- crossprod.spam(A.inv, trans.inf.matrix[[1]] %*% A.inv)
Sigma.star <- trans.inf.matrix[[1]]
for ( j in 2:250){
  D <- diag(c(rep(1,4), exp(eta_hat_mat[5,j]), g.phi.der((eta_hat_mat[6,j])), 
              exp(eta_hat_mat[7,j])))
  trans.inf.matrix[[j]] <- crossprod.spam(D, (info.mat_list[[j]]) %*% D )
  trans.inf.matrix[[j]] <- crossprod.spam(A.inv, trans.inf.matrix[[j]] %*% A.inv)
  Sigma.star <- bdiag.spam(Sigma.star, trans.inf.matrix[[j]])
}


############################ DEFINING ADJACENCY MATRIX ########################
locs <- as.matrix(annual_max_temp_1[,1:2])
D <- rdist(locs)
tol <- 1.5
grid.space <- min(D[D > 0])
grid.space
W <- (abs(D - grid.space) < tol) * 1
diag(W) <- 0
Adj_mat <- diag(rowSums(W))-W
Adj_mat <- as.spam(Adj_mat)

####################### DESIGN MATRIX FORMATION ###############################

Design_Mat_list <- list()
locs <- as.data.frame(locs)
elevs <- read.csv("USA_1deg_elevation_stats.csv")
elevs <- as.data.frame(elevs)
mean.sea.dist <- as.data.frame(read.csv("USA_1deg_mean_sea_distance.csv"))

for ( i in 1:length(locs$lon)){
  Design_Mat_list[[i]] <- kronecker(diag(7), t(c(1,locs$lon[i],locs$lat[i], 
                                                 elevs$mean_elev[i],
                                                 elevs$sd_elev[i],
                                                 mean.sea.dist$mean_sea_dist_m[i])))
  
}
X <- do.call(rbind, Design_Mat_list)
dim(X)

save(X, file = "Design_Matrix.RData")
save(Sigma.star, file = "Sigma_star.RData")
save(Adj_mat, file="Adjacency_Matrix.RData")