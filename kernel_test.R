# Author: Mathis Messager
# Creation date: June 2018
# Purpose: Visualize and create weighted kernels to create traffic-related pollution heatmaps using ArcMap Focal Statistics tool 
# Link to Arcmap tool: http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/focal-statistics.htm

library(raster)
library(plyr)
library(ggplot2)
library(data.table)

#### ------ Get example traffic map to get resolution from ----- #####
res <- 'C:/Mathis/ICSL/stormwater/results/bing'
reclassmlc <- file.path(res,"181204_02_00_class_mlc.tif")

#### ------ Function to produce kernel ----- #####
weightedkernel <- function(raster_file, maxdist, shape, p, outkernel) {
  #Read in raster and get x and y cell size
  r <- raster(raster_file)
  Xres <- xres(r)
  Yres <- yres(r)
  
  #Create matrix of distance from center pixel
  rad <- max(ceiling(maxdist/Xres),ceiling(maxdist/Yres)) #Compute radius of kernel
  matrix_row <- matrix(abs(rep(seq(0, rad * 2) - rad,
                               rad + 1)),
                       nrow = rad * 2 + 1, ncol = rad * 2 + 1) #Matrix of East-West distances from center axis
  
  matrix_col <- matrix(abs(rep(seq(0, rad * 2) - rad,
                               each = rad * 2 + 1)),
                       nrow = rad * 2 + 1, ncol = rad * 2 + 1) #Matrix of North-South distances from center axis
  
  matrix_dist <- sqrt((Yres*matrix_row)^2+(Xres*matrix_col)^2) #Matrix of total distance from center pixel
  
  #Compute relative pollution index
  if (shape=='linear')   matrix_weight <- 1-(matrix_dist/maxdist)
  if (shape=='power')  {
    #print('power')
    pow <- 1/p
    matrix_weight <- 1-(matrix_dist^pow / maxdist^pow)
  } 
  if (shape=='log') {
    #print('log')
    matrix_weight <- 1-(log(matrix_dist+1) / log(maxdist+1))
  }
  if (shape=='quartic_kernel') { #Check this https://gis.stackexchange.com/questions/32300/what-type-of-kernel-density-function-does-arcmap-use
    matrix_weight <- 3/pi*(1-(matrix_dist^2)/(maxdist^2))^2
  }
  
  #Make sure that all values are >= 0
  matrix_weight[matrix_dist>maxdist] <- 0
  #Create kernel header for ArcMap
  matrix_weight <- rbind(c(rad * 2 + 1, rad * 2 + 1, rep(NA, rad*2-1)),
                         matrix_weight)
  #Write out kernel to txt file
  write.table(matrix_weight, file=outkernel,col.names=FALSE,row.names=FALSE,quote=FALSE,na="")
}

#### ------ Run function for a range of maximum distances for pollution to travel and kernel shapes ----- #####
distlist <- c(10, 50, 100, 200, 300, 500)
powlist <- c(1,2,3) #powers to apply for kernel
for (i in distlist) {
  logout <- file.path(res, paste0('kernel_log', i, '.txt'))
  print(logout)
  weightedkernel(reclassmlc, i, 'log', outkernel=logout)
  for (j in powlist) {
    powout <- file.path(res, paste0('kernel_pow', i, '_', j, '.txt'))
    print(powout)
    weightedkernel(reclassmlc, i, 'power', p = j,  outkernel=powout)
  }
}


#### ------ Visualize all kernels in 2-D space ----- #####
#Create distance values for all kernel sizes and powers
x <- ldply(distlist, function(x) {
  seql <- seq(-x,x)
  if (length(seql)<2*max(distlist)){
    pad <- rep(0,(2*max(distlist)-(length(seql)-1))/2)
  } else {
    pad <- NULL
  }
  return(c(pad,seql, pad))
  })
x <- t(x)
colnames(x) <- distlist
xmelt <- melt(x, var.name = 'distlist')
setDT(xmelt)[xmelt$value == 0, value := NA]
xbind <- data.table(do.call("rbind", replicate(length(powlist), xmelt, simplify = FALSE)),
                    p = rep(powlist, each=nrow(xmelt)))

#Run 2-d kernel generation
xbind[,ypow:= 1-((abs(value)^(1/p))/(Var2^(1/p)))] #Power
xbind[,ylog2:= 1-(log2(abs(value)+1))/(log2(Var2+1))] #Logarithm
#Plot 'em
ggplot(xbind, aes(x=value, y=ypow, color = interaction(p, Var2))) +
  geom_line(size=1.2) +
  geom_line(aes(y= ylog2, group=factor(Var2)), color = 'black', size=1.2, alpha=1/2) +
  labs(x='Distance (meters)', y='Relative value') + theme_bw()