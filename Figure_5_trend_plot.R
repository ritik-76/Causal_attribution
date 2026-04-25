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


eta.mean.post.mat <- matrix(eta.mean.post, nrow = 7, ncol = 250)

trend.diff <- c()

trend.factual <- c()
for ( k in 1:ncol(eta.mean.post.mat)){
  trend.factual[k] <- eta.mean.post.mat[,k][4]
  trend.diff[k] <- eta.mean.post.mat[,k][4]-eta.mean.post.mat[,k][2]
}
p.causal.trend <- matrix(NA,  nrow =10000 , ncol= 250)
p.causal.trend.diff <- matrix(NA,  nrow =10000 , ncol= 250)

for ( i in 1:nrow(eta.chain.mat)){
  eta.chain.mat.loc <- matrix(eta.chain.mat[i,], nrow = 7, ncol = 250)
  for ( k in 1:ncol(eta.chain.mat.loc)){
    p.causal.trend.diff[i,k] <- eta.chain.mat.loc[,k][4]-eta.chain.mat.loc[,k][2]
    p.causal.trend [i,k] <- eta.chain.mat.loc[,k][4]
  }
}

p.causal.trend.diff.sd <- c()
p.causal.trend.sd <- c()
for ( l in 1:ncol(p.causal.effect.1)){
  
  p.causal.trend.diff.sd[l] <- sd(p.causal.trend.diff[,l])
  p.causal.trend.sd[l] <- sd(p.causal.trend[,l])
}


p3.1 <- map.heatmap(lat = annual_max_temp_1$lat, lon = annual_max_temp_1$lon,
                    data = trend.factual)
p3.2 <- map.heatmap(lat = annual_max_temp_1$lat, lon = annual_max_temp_1$lon,
                    data = p.causal.trend.sd)
p3.3 <- map.heatmap(lat = annual_max_temp_1$lat, lon = annual_max_temp_1$lon,
                    data = trend.diff)
p3.4 <- map.heatmap(lat = annual_max_temp_1$lat, lon = annual_max_temp_1$lon,
                    data = p.causal.trend.diff.sd)
library(cowplot)

p3 <- plot_grid(p3.1, p3.2, p3.3, p3.4, nrow = 2, ncol = 2)