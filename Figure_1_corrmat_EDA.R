########################## CORRELATION MATRIX PLOT ############################
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
  #print(summary(fit.res))
  covariate.fit <- as.matrix(covariate.fit)
  covariate.fit[,1] <- rep(1,250)
  residual.mat[k,] <- eta_hat_mat[k,] - 
    as.matrix(covariate.fit) %*%(as.numeric(coefficients(fit.res)))
  
}



library(ggplot2)
library(reshape2)

#set.seed(1)
mat <- cor(t(residual.mat))

df <- melt(mat)
colnames(df) <- c("x", "y", "value")


q <- quantile(df$value, 0.95)
df$plot_val <- pmin(df$value, q)

res_mat_cor <- ggplot(df, aes(factor(y), factor(x), fill = plot_val)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", value)),
            size = 5,
            fontface = "bold") +
  coord_fixed(ratio = 1) +
  scale_fill_distiller(palette="RdBu", type="div", direction=-1, guide = "colourbar",
                       labels = scales::label_wrap(50),
                       limits= c(-1,1)) +
  labs(x = NULL, y=NULL,fill=NULL) +
  theme_minimal()+ 
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.text=element_text(size=14), 
    legend.key.height = unit(1, "cm"), 
    legend.key.width = unit(1, "cm"), 
    legend.title = element_text(size=14)
  ) +
  scale_x_discrete(labels=c(expression(alpha[0]),expression(alpha[1]),expression(beta[0]),
                            expression(beta[1]),expression(sigma^"*"),expression(psi),
                            expression(lambda^"*")), position = "top")+
  scale_y_discrete(labels=c(expression(lambda^"*"),expression(psi),expression(sigma^"*"),
                            expression(beta[1]),expression(beta[0]),expression(alpha[1]),
                            expression(alpha[0])), limits = rev)