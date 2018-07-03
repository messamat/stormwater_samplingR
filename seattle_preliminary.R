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
#str(trees)
#Create industrial zoning boolean variable
trees[trees$CLASS_DESC=='Manufacturing/Industrial' & !is.na(trees$CLASS_DESC),'industrial'] <- 1
trees[trees$CLASS_DESC!='Manufacturing/Industrial' | is.na(trees$CLASS_DESC),'industrial'] <- 0
#Compute tree density and plot % white people vs. density
setDT(trees)[,treekm2:=length(UNITID)*1000000/ALAND10, .(GEOID10)]
trees$percwhite <- with(trees, DP0080003/DP0080001)
ggplot(trees, aes(x=percwhite, y=treekm2)) + geom_point()
#
treesel <- as.data.frame(trees[grepl('^Acer platanoides.+', trees$SCIENTIFIC) | trees$SCIENTIFIC=='Acer macrophyllum',])
treesel$SCIENTIFIC_simple <- treesel$SCIENTIFIC
treesel[grepl('^Acer platanoides.+', treesel$SCIENTIFIC), 'SCIENTIFIC_simple'] <- 'Acer platanoides' 

########################################################################################################################
# Analyze distribution of traffic values for trees vs. whole city
########################################################################################################################
checkout_distrib=F
if (checkout_distrib==T) {
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
}

#######################################################################################
#AM <- trees[trees$SCIENTIFIC=='Acer macrophyllum',]
#Create sampling design
#For speed limit
treesel$Value <- treesel$heat_spdlm
treespdlm_bins <- bin_rastertab(treesel, nbins=10, tukey=T, rep=100, rastertab=F) #Create power-scaled bins based on speed limit heat value distribution of trees
ggplot() +
  geom_rect(data=treespdlm_bins, aes(xmin=binmin, xmax=binmax, ymin=0, ymax=count/100),fill='lightgrey')+
  scale_y_sqrt(expand=c(0,0)) + 
  scale_x_log10() + 
  theme_classic()
treesel$spdlm_bin <- .bincode(treesel$heat_spdlm, breaks=c(min(treesel$heat_spdlm)-1,
                                                 treespdlm_bins$binmax[1:(nrow(treespdlm_bins)-1)],
                                                 max(treespdlm_bins$binmax)+1), right=T, include.lowest = T) #Assign each tree to a bin 

#For AADT
treesel$Value <- treesel$heat_AADT 
treeaadt_bins <- bin_rastertab(treesel, nbins=10, tukey=T, rep=100, rastertab=F) #Create power-scaled bins based on AADT heat value distribution of trees
treesel$aadt_bin <- .bincode(treesel$heat_AADT, breaks=c(min(treesel$heat_AADT)-1,
                                                 treeaadt_bins$binmax[1:(nrow(treeaadt_bins)-1)],
                                                 max(treeaadt_bins$binmax)+1), right=F, include.lowest = T)

ggplot(treesel, aes(x=aadt_bin, y=spdlm_bin, color=factor(industrial))) + geom_jitter() + geom_point(color='black',size=3) +
  labs(x='AADT bin', y='Speed limit bin')

#For congestion
treesel$Value <- treesel$heat_bing_
treebing_bins <- bin_rastertab(treesel, nbins=10, tukey=T, rep=100, rastertab=F)
treesel$bing_bin <- .bincode(treesel$heat_bing_, breaks=c(min(treesel$heat_bing_-1),
                                               treebing_bins$binmax[1:(nrow(treebing_bins)-1)],
                                                 max(treebing_bins$binmax)+1), right=F, include.lowest = T)

ggplot(treesel, aes(x=bing_bin, y=aadt_bin, color=factor(spdlm_bin))) + geom_jitter() + geom_point(color='black',size=3) +
  labs(x='Congestion bin', y='AADT bin')

#For % white-only people
# treesel$Value <-treesel$percwhite
# treewhite_bins <- bin_rastertab(treesel, nbins=5, tukey=T, rep=100, rastertab=F)
# treesel$white_bin <- .bincode(treesel$percwhite, breaks=c(min(treesel$percwhite)-1,
#                                                treewhite_bins$binmax[1:(nrow(treewhite_bins)-1)],
#                                                max(treewhite_bins$binmax)+1), right=T, include.lowest = T)
# ggplot(treesel, aes(x=aadt_bin, y=spdlm_bin, color=factor(white_bin))) + geom_jitter() + geom_point(color='black',size=3)
# 
# #For zoning
# length(which(treesel$industrial==1))
##############################################################################################
#Devise multivariate sampling procedure
treeselcopy <- treesel

nrow(unique(treesel[,c('aadt_bin', 'spdlm_bin', 'bing_bin')])) #Number of unique variable bin combinations

#Get three trees with lowest and highest aadt, speed limit, and congestion heat value. 
ID='UNITID'
n=3
treesel[,ID] <- as.character(treesel[,ID])
treesel$tosample <- NA
#For each variable
for (var in c('heat_AADT','heat_spdlm','heat_bing_')) {
  print(var)
  mintreeID <- NULL
  i=n
  while (i > 0){
    mintrees <- treesel[(treesel[,var] %in% 
                           min(treesel[!(treesel[,ID] %in% mintreeID),var])),] #Get all trees with minimum variable value that is above that of the trees already sampled
    samplesize <- ifelse(nrow(mintrees)<=i,nrow(mintrees),i) #Determine sample size, as the minimum of [number of trees with desired value, remaining number of trees to sample]
    mintreeID[(n-i+1):(n-i+samplesize)] <- sample(x=as.character(mintrees[,ID]), size=1, replace=F) #Get ID of random sample of trees with desired value
    i=i-samplesize #Adjust remaining number of trees to sample
  }
  treesel[treesel[,ID] %in% mintreeID & !is.na(treesel$tosample),'tosample'] <-  paste0(treesel[treesel[,ID] %in% mintreeID & !is.na(treesel$tosample),'tosample'], 'min',var,'_')
  treesel[treesel[,ID] %in% mintreeID & is.na(treesel$tosample),'tosample'] <-  paste0('min',var,'_')
  
  #Same for max value
  maxtreeID <- NULL
  i=n
  while (i > 0){
    maxtrees <- treesel[(treesel[,var] %in% 
                           max(treesel[!(treesel[,ID] %in% maxtreeID),var])),] #Get all trees with maximum variable value that is above that of the trees already sampled
    samplesize <- ifelse(nrow(maxtrees)<=i,nrow(maxtrees),i) #Determine sample size, as the minimum of [number of trees with desired value, remaining number of trees to sample]
    maxtreeID[(n-i+1):(n-i+samplesize)] <- sample(x=as.character(maxtrees[,ID]), size=1, replace=F) #Get ID of random sample of trees with desired value
    i=i-samplesize #Adjust remaining number of trees to sample
  }
  treesel[treesel[,ID] %in% maxtreeID & !is.na(treesel$tosample),'tosample'] <-  paste0(treesel[treesel[,ID] %in% maxtreeID & !is.na(treesel$tosample),'tosample'], 'max',var,'_')
  treesel[treesel[,ID] %in% maxtreeID & is.na(treesel$tosample),'tosample'] <-  paste0('max',var,'_')
}


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
