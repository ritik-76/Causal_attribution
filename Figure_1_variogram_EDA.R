########################## VARIOGRAM PLOT #####################################
source("Likelihhod_estimation.R")


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
  #print(summary(fit.res))
  covariate.fit <- as.matrix(covariate.fit)
  covariate.fit[,1] <- rep(1,250)
  residual.mat[k,] <- eta_hat_mat[k,] - 
    as.matrix(covariate.fit) %*%(as.numeric(coefficients(fit.res)))
  
}

d_max <- 15
d <- seq(0,d_max, length = 20)
vg <- variog(coords = as.matrix(annual_max_temp_1[,(1:2)]),data = t(residual.mat),
             uvec = d)


for ( i in 1:ncol(vg$v)){
  vg$v[,i] <- vg$v[,i]/var(residual.mat[i,])
}

vg_mat <- vg$v[-1,]

df_vg <- data.frame(vg_mat)
colnames(df_vg) <- c("alpha0","alpha1","beta0","beta1","sigma","xi","lambda")
df_vg$x <- vg$u[-1]

params <- c("alpha0","alpha1","beta0","beta1","sigma","psi","lambda")

param_labels <- c(
  expression(alpha[0]),
  expression(alpha[1]),
  expression(beta[0]),
  expression(beta[1]),
  expression(sigma^"*"),  # if you want σ*
  expression(psi),
  expression(lambda^"*")
)

library(ggplot2)

var_plot <- ggplot(df_vg, aes(x)) +
  
  geom_line(aes(y = alpha0, color = "alpha0", group = "alpha0")) +
  geom_point(aes(y = alpha0, color = "alpha0", shape = "alpha0"), size = 2.4) +
  
  geom_line(aes(y = alpha1, color = "alpha1", group = "alpha1")) +
  geom_point(aes(y = alpha1, color = "alpha1", shape = "alpha1"), size = 2.4) +
  
  geom_line(aes(y = beta0,  color = "beta0",  group = "beta0")) +
  geom_point(aes(y = beta0,  color = "beta0",  shape = "beta0"), size = 2.4) +
  
  geom_line(aes(y = beta1,  color = "beta1",  group = "beta1")) +
  geom_point(aes(y = beta1,  color = "beta1",  shape = "beta1"), size = 2.4) +
  
  geom_line(aes(y = sigma,  color = "sigma",  group = "sigma")) +
  geom_point(aes(y = sigma,  color = "sigma",  shape = "sigma"), size = 2.4) +
  
  geom_line(aes(y = xi,     color = "psi",    group = "psi")) +
  geom_point(aes(y = xi,     color = "psi",    shape = "psi"), size = 2.4) +
  
  geom_line(aes(y = lambda, color = "lambda", group = "lambda")) +
  geom_point(aes(y = lambda, color = "lambda", shape = "lambda"), size = 2.4)+
  
  scale_color_manual(
    values = c(
      "alpha0"="#1b9e77","alpha1"="#d95f02",
      "beta0"="#7570b3","beta1"="#e7298a",
      "sigma"="#66a61e","psi"="#e6ab02",
      "lambda"="#a6761d"
    ),
    breaks = params,
    labels = param_labels
  ) +
  
  scale_shape_manual(
    values = c(
      "alpha0"=1,"alpha1"=0,
      "beta0"=2,"beta1"=5,
      "sigma"=6,"psi"=15,
      "lambda"=16
    ),
    breaks = params,
    labels = param_labels
  ) +
  theme_grey(base_size = 10) +
  labs(x = "Distance(deg)", y = "Variogram", color = "", shape = "") +
  theme(
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    
    legend.position = c(0.99, 0.01),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(color = "black"),
    legend.key = element_rect(fill = NA),
    
    legend.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18),
    legend.text = element_text(size = 14,face="bold")
  )
