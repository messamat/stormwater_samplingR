library(ggplot2)
library(data.table)
library(foreign)

setwd("C:/Mathis/ICSL/stormwater/src/stormwater_samplingR")
source('internal_functions.R')

setwd("C:/Mathis/ICSL/stormwater/results/")
trees <- read.dbf('trees_tab.dbf') #Traffic values for Big Leaf Maple
aadt <- read.dbf('heat_aadt.dbf') #Raster attribute table for Average Annual Daily Traffic heat map (AADT)
spdlm <- read.dbf('heat_spdlm.dbf') #Raster attribute table for speed limit heat map
streets <- read.dbf('Seattle_roads.dbf')

########################################################################################################################
# Format data
########################################################################################################################
str(trees)
#Create industrial zoning boolean variable
trees[trees$CLASS_DESC=='Manufacturing/Industrial' & !is.na(trees$CLASS_DESC),'industrial'] <- 1
trees[trees$CLASS_DESC!='Manufacturing/Industrial' | is.na(trees$CLASS_DESC),'industrial'] <- 0
#Compute tree density and plot % white people vs. density
setDT(trees)[,treekm2:=length(UNITID)*1000000/ALAND10, .(GEOID10)]
trees$percwhite <- with(trees, DP0080003/DP0080001)
ggplot(trees, aes(x=percwhite, y=treekm2)) + geom_point()
#
treesel <- trees[grepl('^Acer platanoides.+', trees$SCIENTIFIC) | trees$SCIENTIFIC=='Acer macrophyllum',]
treesel$SCIENTIFIC_simple <- treesel$SCIENTIFIC
treesel[grepl('^Acer platanoides.+', treesel$SCIENTIFIC), 'SCIENTIFIC_simple'] <- 'Acer platanoides' 

########################################################################################################################
# Analyze distribution of traffic values for trees vs. whole city
########################################################################################################################
#Speed limit
spdlm$Value <- spdlm$Value+1
spdlm_hist <- bin_rastertab(spdlm, 50, tukey=T, rep=100, rastertab=T)
ggplot() +
  geom_rect(data=spdlm_hist, aes(xmin=binmin, xmax=binmax, ymin=0, ymax=count/100),fill='lightgrey')+
  geom_histogram(data=treesel,aes(x=heat_spdlm, fill=SCIENTIFIC_simple), bins=100, alpha=1/2, position='identity')+
  scale_y_sqrt(expand=c(0,0)) + 
  scale_x_log10(breaks=c(1,10,100,1000,10000)) + 
  theme_classic()

#AADT with logarithmic bins
aadt$Value <- aadt$Value+1
aadt_hist <- bin_rastertab(aadt, 50, tukey=T, rep=100, rastertab=T)
ggplot() +
  geom_rect(data=aadt_hist, aes(xmin=binmin, xmax=binmax, ymin=0, ymax=count/100),fill='lightgrey')+
  geom_histogram(data=treesel,aes(x=heat_AADT, fill=SCIENTIFIC_simple), bins=100, alpha=1/2, position='identity')+
  scale_y_sqrt(expand=c(0,0)) + 
  scale_x_log10() + 
  theme_classic()

#Congestion

#% white-only people
ggplot() +
  geom_histogram(data=treesel,aes(x=percwhite, fill=SCIENTIFIC_simple), bins=100, alpha=1/2, position='identity')+
  scale_y_sqrt(expand=c(0,0)) + 
  scale_x_log10() + 
  theme_classic()


#######################################################################################
#Create sampling design
#For speed limit
treesel$Value <- treesel$heat_spdlm
treespdlm_bins <- bin_rastertab(treesel, nbins=20, tukey=T, rep=100, rastertab=F)
ggplot() +
  geom_rect(data=treespdlm_bins, aes(xmin=binmin, xmax=binmax, ymin=0, ymax=count/100),fill='lightgrey')+
  scale_y_sqrt(expand=c(0,0)) + 
  scale_x_log10() + 
  theme_classic()
treesel$spdlm_bin <- .bincode(treesel$heat_spdlm, breaks=c(min(treesel$heat_spdlm),treespdlm_bins$binmax))

#For AADT
treesel$Value <- treesel$heat_AADT
treeaadt_bins <- bin_rastertab(treesel, nbins=20, tukey=T, rep=100, rastertab=F)
treesel$aadt_bin <- .bincode(treesel$heat_AADT, breaks=c(min(treesel$heat_AADT),treeaadt_bins$binmax))

ggplot(treesel, aes(x=aadt_bin, y=spdlm_bin)) + geom_point()

#For congestion

#For % white-only people
treesel$Value <-treesel$percwhite
treewhite_bins <- bin_rastertab(treesel, nbins=20, tukey=T, rep=100, rastertab=F)
treesel$white_bin <- .bincode(treesel$percwhite, breaks=c(min(treesel$percwhite),treewhite_bins$binmax))

#For zoning
length(which(treesel$industrial==1))



#######################################################################################
# Analyze how functional class relates to AADT and Speed limit
#http://streetsillustrated.seattle.gov/street-types/street-classification/

str(streets)
#Order factors
streets$ARTDESCRIP <- factor(streets$ARTDESCRIP, 
                             levels=c('Not Designated', 'Collector Arterial','Minor Arterial', 'Principal Arterial', 'State Route/Freeway','Interstate/Freeway'))

#Arterial classification only slightly discriminate among speed limits
ggplot(streets[!(streets$ARTDESCRIP %in% c('County Arterial', NA)),], aes(x=SPEEDLIMIT, fill=ARTDESCRIP)) + 
  geom_histogram()+
  facet_wrap(~ARTDESCRIP, scales='free_y', ncol=1) +
  theme_bw()

ggplot(streets[!(streets$ARTDESCRIP %in% c('County Arterial', NA)),], aes(x=ARTDESCRIP, y=SPEEDLIMIT, fill=ARTDESCRIP)) + 
  geom_boxplot() +
  theme_classic()

#Arterial classification  only slightly discriminate among speed limits
ggplot(streets[!(streets$ARTDESCRIP %in% c('County Arterial', NA)),], aes(x=ARTDESCRIP, y=AADT_inter, fill=ARTDESCRIP)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  theme_classic()

##############################################################################################
#Devise multivariate sampling procedure

