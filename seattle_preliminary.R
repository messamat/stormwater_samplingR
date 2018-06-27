library(ggplot2)
library(data.table)
library(foreign)

setwd("C:/Mathis/ICSL/stormwater/results/")
BLMtrees <- read.dbf('BLMtrees_tab.dbf') #Traffic values for Big Leaf Maple
NMtrees <- read.dbf('NMtrees_tab.dbf') #Traffic values for Norway Maple
aadt <- read.dbf('heat_aadt.dbf') #Raster attribute table for Average Annual Daily Traffic heat map (AADT)
spdlm <- read.dbf('heat_spdlm.dbf') #Raster attribute table for speed limit heat map
streets <- read.dbf('Seattle_roads.dbf')

#Function to bin raster attribute table
bin_rastertab<- function(tab, binsize) {
  minbin = min(tab$Value)-min(tab$Value)%%binsize+binsize
  maxbin = max(tab$Value)-max(tab$Value)%%binsize+binsize
  bins <- data.frame(bin=c(with(tab,rep(seq(minbin,maxbin, binsize),each=binsize)),maxbin),
                     val=with(tab,seq(minbin-binsize,maxbin)))
  tab_hist <- merge(tab, bins, by.x='Value', by.y='val')
  tab_hist <- setDT(tab_hist)[,sum(Count), .(bin)]
}

########################################################################################################################
# Analyze distribution of traffic values for trees vs. whole city
########################################################################################################################
#-----------------BIG LEAF MAPLE --------------
#Speed limit
ggplot(BLMtrees, aes(x=heat_spdlm)) +
  geom_ribbon(data=spdlm, aes(x=Value, ymin=0, ymax=Count/10),stat='identity', fill='darkgrey') +
  geom_histogram(fill='red', bins=100, alpha=1/2)+
  geom_histogram(data=NMtrees, fill='blue', bins=100, alpha=1/5)+
  scale_y_sqrt() + 
  theme_classic()

#AADT
aadt_hist <- bin_rastertab(aadt, 200)
ggplot(BLMtrees, aes(x=heat_AADT)) +
  geom_ribbon(data=aadt_hist, aes(x=bin, ymin=0, ymax=V1/10),stat='identity', fill='darkgrey') +
  geom_histogram(fill='red', bins=100, alpha=1/3)+
  geom_histogram(data=NMtrees, fill='blue', bins=100, alpha=1/5)+
  scale_y_sqrt(expand=c(0,0)) + 
  scale_x_continuous(expand=c(0,0))+
  theme_classic()


#######################################################################################
# Analyze how functional class relates to AADT and Speed limit
#http://streetsillustrated.seattle.gov/street-types/street-classification/






