res <- 'F:/Levin_Lab/stormwater/results/bing'
reclassmlc <- file.path(res,"180617_14_30_reclass_mlc.tif")
# dist <- 250
# shape='log'
# p = 3
# outkernel = 'F:/Levin_Lab/stormwater/results/logkernel.txt'

weightedkernel <- function(raster_file, maxdist, shape, p, outkernel) {
  r <- raster(raster_file)
  Xres <- xres(r)
  Yres <- yres(r)
  rad <- max(ceiling(maxdist/Xres),ceiling(maxdist/Yres))
  matrix_row <- matrix(abs(rep(seq(0,rad*2)-rad,rad+1)),nrow=rad*2+1,ncol=rad*2+1)
  matrix_col <- matrix(abs(rep(seq(0,rad*2)-rad,each=rad*2+1)),nrow=rad*2+1,ncol=rad*2+1)
  matrix_dist <- sqrt((Yres*matrix_row)^2+(Xres*matrix_col)^2)
  
  if (shape=='linear')   matrix_weight <- 1-(matrix_dist/maxdist)
  if (shape=='power')  {
    print('power')
    pow <- 1/p
    matrix_weight <- 1-(matrix_dist^pow / maxdist^pow)
  } 
  if (shape=='log') {
    print('log')
    matrix_weight <- 1-(log(matrix_dist+1) / log(maxdist+1))
  }
  if (shape=='quartic_kernel') { #Check this https://gis.stackexchange.com/questions/32300/what-type-of-kernel-density-function-does-arcmap-use
    matrix_weight <- 3/pi*(1-(matrix_dist^2)/(maxdist^2))^2
  }
  matrix_weight[matrix_dist>dist] <- 0
  matrix_weight <- rbind(c(rad*2+1,rad*2+1,rep(NA, rad*2-1)),matrix_weight)
  write.table(matrix_weight, file=outkernel,col.names=FALSE,row.names=FALSE,quote=FALSE,na="")
}

weightedkernel(reclassmlc, 100, 'log', outkernel='F:/Levin_Lab/stormwater/results/logkernel100.txt')
weightedkernel(reclassmlc, 100, 'power',p=2, outkernel='F:/Levin_Lab/stormwater/results/powerkernel100.txt')



# #Curve illustrations
# dist=100
# p <- 2
# x <- seq(-dist,dist)
# #Power
# pow <- 1/p
# y <- 1-((abs(x)^pow)/(dist^pow))
# qplot(x,y)
# #Logarithm
# y <- 1-(log2(abs(x)+1))/(log2(dist+1))
# qplot(x,y)
# #

