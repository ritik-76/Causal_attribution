
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


eta.mean.post.mat <- matrix(eta.mean.post, nrow = 7, ncol = 250)

p.causal.effect <- c()

for ( k in 1:ncol(eta.mean.post.mat)){
  diff.intercept <- eta.mean.post.mat[,k][3]-eta.mean.post.mat[,k][1]
  p.causal.effect[k] <- diff.intercept
  
}

p.causal.effect.1 <- matrix(NA,  nrow =10000 , ncol= 250)

for ( i in 1:nrow(eta.chain.mat)){
  eta.chain.mat.loc <- matrix(eta.chain.mat[i,], nrow = 7, ncol = 250)
  for ( k in 1:ncol(eta.chain.mat.loc)){
    diff.intercept <- eta.chain.mat.loc[,k][3]-eta.chain.mat.loc[,k][1]
    p.causal.effect.1[i,k] <- diff.intercept
   
  }
}

p.causal.sd <- c()
for ( l in 1:ncol(p.causal.effect.1)){
  p.causal.sd[l] <- sd(p.causal.effect.1[,l])
}


p0.1 <- map.heatmap(lat = annual_max_temp_1$lat, lon = annual_max_temp_1$lon,
                    data = p.causal.effect)
p0.2 <- map.heatmap(lat = annual_max_temp_1$lat, lon = annual_max_temp_1$lon,
                    data = p.causal.sd)

p0 <- plot_grid(p0.1, p0.2, nrow = 1, ncol = 2)