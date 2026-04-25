
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

library(ggmap)
# register_google("API_Key")
map.heatmap <- function (lat, lon, data, mainTitle = NULL, legendTitle = NULL, 
                         xlim = NULL, ylim = NULL, zlim = NULL){
  
  # Set limits for x, y, z if not specified as parameters
  if(is.null(xlim)){xlim <- range(lon, na.rm = TRUE)}
  if(is.null(ylim)){ylim <- range(lat, na.rm = TRUE)}
  if(is.null(zlim)){zlim <- range(data, na.rm = TRUE)}
  
  usabox <- make_bbox(lon = c(min(lon)-1, max(lon)+1), lat = c(max(lat)+1, min(lat)-1), f =.1)
  us <- get_map(location = usabox, zoom = 2, maptype = "terrain", source = "google")
  dataframe <- data.frame(x = lon, y = lat, data = data)
  
  p <- ggmap(us)
  p <- p + geom_tile(dataframe, mapping = aes(x = x, y = y, fill = data), 
                     alpha = 1) + # 
    labs(title=mainTitle, x="Longitude", y="Latitude") + 
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=16),
          plot.title = element_text(size=20, hjust = 0.5),
          legend.text = element_text(size = 16), 
          legend.title = element_text(size = 18), 
          legend.key.size = unit(0.3, "in")) + 
    scale_fill_gradientn(name = legendTitle,
                         colours = c("#2c7bb6", "#ffffbf", "#d7191c"))+ 
    xlim(xlim+c(-2,3)) + ylim(ylim+c(-2,2))
  
  baseData <- map_data("state")
  p <- p + geom_path(data=baseData, aes(x=long, y=lat, group=group), 
                     colour="black")
  
  return(p)
}


p.causal.effect.pi <- c()
p.causal.effect.pi.sd <- c()
t.pi <- c(1850:1900)
t.star <- sum((t.pi-mean(t))/sd(t))
for ( k in 1:ncol(eta.mean.post.mat)){
  s <- 0
  for (l in t.pi){
    diff.intercept.pi <- eta.mean.post.mat[,k][3]-eta.mean.post.mat[,k][1]
    diff.trend.pi <- (l-mean(t))/sd(t) * (eta.mean.post.mat[,k][4]-eta.mean.post.mat[,k][2])
    s <- s + diff.intercept.pi + diff.trend.pi
    
  }
  p.causal.effect.pi <- append(p.causal.effect.pi,s)
}
p.causal.effect.pi <- p.causal.effect.pi/51

p.causal.effect.pi.1 <- matrix(NA,  nrow =10000 , ncol= 250)
for ( i in 1:nrow(eta.chain.mat)){
  eta.chain.mat.loc <- matrix(eta.chain.mat[i,], nrow = 7, ncol = 250)
  for ( k in 1:ncol(eta.chain.mat.loc)){
    diff.intercept.post.pi <- eta.chain.mat.loc[,k][3]-eta.chain.mat.loc[,k][1]
    diff.trend.post.pi <- t.star/51 * (eta.chain.mat.loc[,k][4]-eta.chain.mat.loc[,k][2])
    p.causal.effect.pi.1[i,k] <- diff.intercept.post.pi + diff.trend.post.pi
  }
}
for (l in 1:ncol(p.causal.effect.pi.1)){
  p.causal.effect.pi.sd[l] <- sd(p.causal.effect.pi.1[,l])
}


p1.1 <- map.heatmap(lat = annual_max_temp_1$lat, lon = annual_max_temp_1$lon,
                    data = p.causal.effect.pi)
p1.2 <- map.heatmap(lat = annual_max_temp_1$lat, lon = annual_max_temp_1$lon,
                    data = p.causal.effect.pi.sd)

p1 <- plot_grid(p1.1, p1.2, nrow = 1, ncol = 2)