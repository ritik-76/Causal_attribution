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
  
  usabox <- make_bbox(lon = c(min(lon), max(lon)), lat = c(max(lat), min(lat)), f =.1)
  us <- get_map(location = usabox, zoom = 3, maptype = "terrain", source = "google")
  dataframe <- data.frame(x = lon, y = lat, data = data)
  
  p <- ggmap(us)
  p <- p + geom_tile(dataframe, mapping = aes(x = x, y = y, fill = data), 
                     alpha = 1) + # 
    labs(title=mainTitle, x="Longitude", y="Latitude") + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          plot.title = element_text(size=20, hjust = 0.5),
          legend.text = element_text(size = 14), 
          legend.title = element_text(size = 14), 
          legend.key.size = unit(0.3, "in")) + 
    scale_fill_manual(
      values = c(
        "selected" = "#1b9e77",   # teal
        "not.selected" = "#eeeeee"
      ),
      guide = "none"
    ) + 
    xlim(xlim+c(-2,3)) + ylim(ylim+c(-2,2))
  
  baseData <- map_data("state")
  p <- p + geom_path(data=baseData, aes(x=long, y=lat, group=group), 
                     colour="black")
  
  return(p)
}

p.causal.effect.1 <- matrix(NA,  nrow =10000 , ncol= 250)

for ( i in 1:nrow(eta.chain.mat)){
  eta.chain.mat.loc <- matrix(eta.chain.mat[i,], nrow = 7, ncol = 250)
}


for ( i in 1:nrow(eta.chain.mat)){
  eta.chain.mat.loc <- matrix(eta.chain.mat[i,], nrow = 7, ncol = 250)
  for ( k in 1:ncol(eta.chain.mat.loc)){
    diff.intercept <- eta.chain.mat.loc[,k][3]-eta.chain.mat.loc[,k][1]
    p.causal.effect.1[i,k] <- diff.intercept
    
  }
}

library(ExceedanceTools)

z.statistics.0 <- c()
z.statistics.1 <- c()

u <- 0.35
u1 <- 0.65
for ( i in 1:ncol(p.causal.effect.1)){
  z.statistics.0[i] <- mean(p.causal.effect.1[,i]-u)/
    (mean((p.causal.effect.1[,i]-u)^2))^0.5
  z.statistics.1[i] <- mean(p.causal.effect.1[,i]-u1)/
    (mean((p.causal.effect.1[,i]-u1)^2))^0.5
  
}

causal.credible.uc <- confreg( obj = t(p.causal.effect.1),
                               level = u,
                               statistic = z.statistics.0,
                               conf.level = 0.95,
                               direction = ">",
                               type = "o",
                               method = "test",
                               greedy = FALSE)
causal.credible.uc1 <- confreg( obj = t(p.causal.effect.1),
                                level = u1,
                                statistic = z.statistics.1,
                                conf.level = 0.95,
                                direction = ">",
                                type = "o",
                                method = "test",
                                greedy = FALSE)





uc_region <- c()
uc_region_1 <- c()

for ( i in 1:nrow(annual_max_temp_1)){
  uc_region[i] <- ifelse(i %in% causal.credible.uc$confidence, "selected", "not.selected")
  uc_region_1[i] <- ifelse(i %in% causal.credible.uc1$confidence, "selected", "not.selected")
  
}


p4.1 <- map.heatmap(lat = annual_max_temp_1$lat, lon = annual_max_temp_1$lon,
                    data = uc_region)
p4.2 <- map.heatmap(lat = annual_max_temp_1$lat, lon = annual_max_temp_1$lon,
                    data = uc_region_1)

library(cowplot)

p4 <- plot_grid(p4.1, p4.2, nrow = 1, ncol = 2)
