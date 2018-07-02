library(ggplot2)
library(data.table)
library(foreign)

setwd("C:/Mathis/ICSL/stormwater/src/stormwater_samplingR")
source('internal_functions.R')

setwd("C:/Mathis/ICSL/stormwater/results/")
trees <- read.dbf('trees_tab.dbf') #Traffic values for Seattle street trees
aadt <- read.dbf('heat_aadt.dbf') #Raster attribute table for Average Annual Daily Traffic heat map (AADT)
spdlm <- read.dbf('heat_spdlm.dbf') #Raster attribute table for speed limit heat map
bing <- read.dbf('heat_bing.dbf')
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
  scale_x_log10(name='Speed limit heat value',breaks=c(1,10,100,1000,10000)) + 
  guides(fill=guide_legend(title="Tree species")) + 
  theme_classic() + 
  theme(legend.position=c(0.20,0.8))

#AADT with logarithmic bins
aadt$Value <- aadt$Value+1
aadt_hist <- bin_rastertab(aadt, 50, tukey=T, rep=100, rastertab=T)
ggplot() +
  geom_rect(data=aadt_hist, aes(xmin=binmin, xmax=binmax, ymin=0, ymax=count/100),fill='lightgrey')+
  geom_histogram(data=treesel,aes(x=heat_AADT, fill=SCIENTIFIC_simple), bins=100, alpha=1/2, position='identity')+
  scale_y_sqrt(expand=c(0,0)) + 
  scale_x_log10(name='AADT heat value') + 
  guides(fill=guide_legend(title="Tree species")) + 
  theme_classic() + 
  theme(legend.position=c(0.20,0.8))

#Congestion
#AADT with logarithmic bins
bing$Value <- bing$Value+1
bing_hist <- bin_rastertab(bing, 100, tukey=T, rep=100, rastertab=T)
ggplot() +
  geom_rect(data=bing_hist, aes(xmin=binmin, xmax=binmax, ymin=0, ymax=count/10000),fill='lightgrey',size=3)+
  #geom_histogram(data=bing,aes(x=Value), bins=100, alpha=1/2, position='identity')+
  geom_histogram(data=treesel,aes(x=heat_bing_, fill=SCIENTIFIC_simple), bins=100, alpha=1/2, position='identity')+
  scale_y_sqrt(expand=c(0,0)) + 
  scale_x_sqrt(name='bing heat value') + 
  guides(fill=guide_legend(title="Tree species")) + 
  theme_classic() + 
  theme(legend.position=c(0.20,0.8))


#% white-only people
ggplot() +
  geom_histogram(data=treesel,aes(x=percwhite, fill=SCIENTIFIC_simple), bins=100, alpha=1/2, position='identity')+
  scale_y_sqrt(expand=c(0,0)) + 
  scale_x_log10() + 
  theme_classic()


#######################################################################################
AM <- trees[trees$SCIENTIFIC=='Acer macrophyllum',]
#Create sampling design
#For speed limit
AM$Value <- AM$heat_spdlm
treespdlm_bins <- bin_rastertab(AM, nbins=20, tukey=T, rep=100, rastertab=F)
ggplot() +
  geom_rect(data=treespdlm_bins, aes(xmin=binmin, xmax=binmax, ymin=0, ymax=count/100),fill='lightgrey')+
  scale_y_sqrt(expand=c(0,0)) + 
  scale_x_log10() + 
  theme_classic()
AM$spdlm_bin <- .bincode(AM$heat_spdlm, breaks=c(min(AM$heat_spdlm)-1,
                                                 treespdlm_bins$binmax[1:(nrow(treespdlm_bins)-1)],
                                                 max(treespdlm_bins$binmax)+1), right=T, include.lowest = T)

#For AADT
AM$Value <- AM$heat_AADT
treeaadt_bins <- bin_rastertab(AM, nbins=20, tukey=T, rep=100, rastertab=F)
AM$aadt_bin <- .bincode(AM$heat_AADT, breaks=c(min(AM$heat_AADT)-1,
                                                 treeaadt_bins$binmax[1:(nrow(treeaadt_bins)-1)],
                                                 max(treeaadt_bins$binmax)+1), right=F, include.lowest = T)

ggplot(AM, aes(x=aadt_bin, y=spdlm_bin, color=factor(industrial))) + geom_jitter() + geom_point(color='black',size=3) +
  labs(x='AADT bin', y='Speed limit bin')

#For congestion
AM$Value <- AM$heat_bing_
treebing_bins <- bin_rastertab(AM, nbins=20, tukey=T, rep=100, rastertab=F)
AM$bing_bin <- .bincode(AM$heat_bing_, breaks=c(min(AM$heat_bing_-1),
                                               treebing_bins$binmax[1:(nrow(treebing_bins)-1)],
                                                 max(treebing_bins$binmax)+1), right=F, include.lowest = T)

ggplot(AM, aes(x=bing_bin, y=aadt_bin, color=factor(industrial))) + geom_jitter() + geom_point(color='black',size=3) +
  labs(x='Congestion bin', y='AADT bin')

#For % white-only people
AM$Value <-AM$percwhite
treewhite_bins <- bin_rastertab(AM, nbins=5, tukey=T, rep=100, rastertab=F)
AM$white_bin <- .bincode(AM$percwhite, breaks=c(min(AM$percwhite)-1,
                                               treewhite_bins$binmax[1:(nrow(treewhite_bins)-1)],
                                               max(treewhite_bins$binmax)+1), right=T, include.lowest = T)
ggplot(AM, aes(x=aadt_bin, y=spdlm_bin, color=factor(white_bin))) + geom_jitter() + geom_point(color='black',size=3)

#For zoning
length(which(AM$industrial==1))


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

#Arterial classification moderately discriminates among AADT levels
ggplot(streets[!(streets$ARTDESCRIP %in% c('County Arterial', NA)),], aes(x=ARTDESCRIP, y=AADT_inter, fill=ARTDESCRIP)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  theme_classic()

##############################################################################################
#Devise multivariate sampling procedure

