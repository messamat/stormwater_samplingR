library(ggplot2)
library(data.table)
library(foreign)
library(sp)
library(plyr)
library(openxlsx)

rootdir = "C:/Mathis/ICSL/stormwater"
setwd(file.path(rootdir, "src/stormwater_samplingR"))
source('internal_functions.R')

setwd(file.path(rootdir,"results"))
trees <- read.dbf('trees_tab.dbf') #Traffic values for Seattle street trees
aadt <- read.dbf('heat_aadt.dbf') #Raster attribute table for Average Annual Daily Traffic heat map (AADT)
spdlm <- read.dbf('heat_spdlm.dbf') #Raster attribute table for speed limit heat map
bing <- read.dbf('heat_bing.dbf')
streets <- read.dbf('Seattle_roads.dbf')

#Sampled trees
fielddata <- read.xlsx(file.path(rootdir,"data/field_data/field_data_raw_20180730_edit.xlsx"), sheet=1,startRow=1, detectDates=T)
#Create Tree.ID.resampling to have the nearest tree when actually sampled tree was not in Seattle database
fielddatasel <- fielddata[fielddata$Include == 'Y' & !is.na(fielddata$Tree.ID.resampling),] #Ignore trees in parks

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


#Select either Norway or Big Leaf Maple, though no hybrid norway maple, and include all already-sampled trees
treesel <- as.data.frame(trees[(grepl('^Acer platanoides.+', trees$SCIENTIFIC)& !grepl('^.+x.+',  trees$SCIENTIFIC)) |
                                 trees$SCIENTIFIC == 'Acer macrophyllum' | 
                                 trees$UNITID %in% fielddatasel,])
treesel$SCIENTIFIC_simple <- treesel$SCIENTIFIC
treesel[treesel$UNITID %in% fielddatasel$Tree.ID.resampling,'sampled'] <- 'Y'

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
treeAADT_bins <- bin_rastertab(treesel, nbins=10, tukey=T, rep=100, rastertab=F) #Create power-scaled bins based on AADT heat value distribution of trees
treesel$AADT_bin <- .bincode(treesel$heat_AADT, breaks=c(min(treesel$heat_AADT)-1,
                                                 treeAADT_bins$binmax[1:(nrow(treeAADT_bins)-1)],
                                                 max(treeAADT_bins$binmax)+1), right=F, include.lowest = T)

ggplot(treesel, aes(x=AADT_bin, y=spdlm_bin, color=factor(industrial))) + geom_jitter() + geom_point(color='black',size=3) +
  labs(x='AADT bin', y='Speed limit bin')

#For congestion
treesel$Value <- treesel$heat_bing_
treebing_bins <- bin_rastertab(treesel, nbins=10, tukey=T, rep=100, rastertab=F)
treesel$bing_bin <- .bincode(treesel$heat_bing_, breaks=c(min(treesel$heat_bing_-1),
                                               treebing_bins$binmax[1:(nrow(treebing_bins)-1)],
                                                 max(treebing_bins$binmax)+1), right=F, include.lowest = T)

ggplot(treesel, aes(x=bing_bin, y=AADT_bin, color=factor(spdlm_bin))) + geom_jitter() + geom_point(color='black',size=3) +
  labs(x='Congestion bin', y='AADT bin')

#For % white-only people
# treesel$Value <-treesel$percwhite
# treewhite_bins <- bin_rastertab(treesel, nbins=5, tukey=T, rep=100, rastertab=F)
# treesel$white_bin <- .bincode(treesel$percwhite, breaks=c(min(treesel$percwhite)-1,
#                                                treewhite_bins$binmax[1:(nrow(treewhite_bins)-1)],
#                                                max(treewhite_bins$binmax)+1), right=T, include.lowest = T)
# ggplot(treesel, aes(x=AADT_bin, y=spdlm_bin, color=factor(white_bin))) + geom_jitter() + geom_point(color='black',size=3)
# 
# #For zoning
# length(which(treesel$industrial==1))


###################################################### Devise multivariate sampling procedure #################################################

#Make 'sampled colum'
#Make sure that no tree will be sampled in the same bin or within 100 m from that tree
treesel <- as.data.frame(treesel)
# Compute pairwise distances between all trees, and for each tree, find its closest neighboring tree
coordmatrix <- as.matrix(treesel[,c('SHAPE_LNG','SHAPE_LAT')])
dist <- as.data.frame(spDists(x = coordmatrix, y = coordmatrix, longlat = T, segments = FALSE, diagonal = FALSE)) #Compute pairwise distance among all trees
colnames(dist) <- as.character(treesel$UNITID)
dist$UNITID <- as.character(treesel$UNITID)
distmelt <-melt(setDT(dist), id.vars='UNITID', variable.name='UNITID.y', value.name = 'dist')
closest <- distmelt[UNITID!=UNITID.y, .SD[which.min(dist)], by = UNITID] #Find closest neighboring tree
treesel <- merge(treesel, closest, by='UNITID') #Merge with main tree dataset


#Create unique bin combination ID
ntotal = 40-length(unique(fielddatasel$`#`)) #Total number of trees to sample
binvars = c('AADT_bin', 'spdlm_bin', 'bing_bin')
treesel$bincomb <- adply(as.data.frame(treesel[,binvars]),1, function(x) paste(x, collapse=""))[,'V1'] #Compute ID for each unique combination of bins

############### Get three trees with lowest and highest aadt, speed limit, and congestion heat value. 
ID='UNITID'
n=3
maxdist = 0.1
set.seed(2)

distmelt_lim <- as.data.frame(distmelt[distmelt$dist<maxdist,])
treesel[,ID] <- as.character(treesel[,ID])
treesel$tosample <- NA
sampledID <- fielddata$Tree.ID.resampling
treesel$heat_bing <- treesel$heat_bing_
#For each variable
for (var in c('heat_AADT','heat_spdlm','heat_bing')) {
  print(var)
  #For min value
  mintreeID <- NULL
  i=n-nrow(treesel[!(is.na(treesel$tosample) & is.na(treesel$sampled)) & treesel[,paste(substr(var, 6, 15),'bin',sep='_')] %in% c(1,2),])
  while (i > 0){
    mintrees <- treesel[(treesel[,var] %in% 
                           min(treesel[!(treesel[,ID] %in% sampledID),var])),] #Get all trees with minimum variable value that is above that of the trees already sampled
    mintrees_sel <- mintrees[!(mintrees$UNITID %in% #Only keep trees that are not within the min distance of already selected trees
                                 as.character(distmelt_lim[distmelt_lim$UNITID %in% sampledID,'UNITID.y'])),] 
    if (nrow(mintrees_sel)==0) {
      sampledID <- c(sampledID, mintrees) #Trees already sampled
    } else {
      #samplesize <- ifelse(nrow(mintrees)<=i,nrow(mintrees),i) #Determine sample size, as the minimum of [number of trees with desired value, remaining number of trees to sample]
      samplesize <- 1
      mintreeID[n-i+samplesize] <- sample(x=as.character(mintrees_sel[,ID]), size=samplesize, replace=F) #Get ID of random sample of trees with desired value
      sampledID <- c(sampledID, mintreeID[(n-i+samplesize)])
      i=i-samplesize #Adjust remaining number of trees to sample
    }
  }
  treesel[treesel[,ID] %in% mintreeID & !is.na(treesel$tosample),'tosample'] <-  paste0(treesel[treesel[,ID] %in% mintreeID & !is.na(treesel$tosample),'tosample'], 'min',var,'_')
  treesel[treesel[,ID] %in% mintreeID & is.na(treesel$tosample),'tosample'] <-  paste0('min',var,'_')

  #Same for max value
  maxtreeID <- NULL
  i=n-nrow(treesel[!(is.na(treesel$tosample) & is.na(treesel$sampled)) & treesel[,paste(substr(var, 6, 15),'bin',sep='_')] == 10,])
  while (i > 0){
    maxtrees <- treesel[(treesel[,var] %in% 
                           max(treesel[!(treesel[,ID] %in% sampledID),var])),] #Get all trees with maximum variable value that is above that of the trees already sampled
    maxtrees_sel <- maxtrees[!(maxtrees$UNITID %in% #Only keep trees that are not within the max distance of already selected trees
                                 as.character(distmelt_lim[distmelt_lim$UNITID %in% sampledID,'UNITID.y'])),] 
    if (nrow(maxtrees_sel)==0) {
      sampledID <- c(sampledID, maxtrees) #Trees already sampled
    } else {
      #samplesize <- ifelse(nrow(maxtrees)<=i,nrow(maxtrees),i) #Determaxe sample size, as the maximum of [number of trees with desired value, remaining number of trees to sample]
      samplesize <- 1
      maxtreeID[n-i+samplesize] <- sample(x=as.character(maxtrees_sel[,ID]), size=samplesize, replace=F) #Get ID of random sample of trees with desired value
      sampledID <- c(sampledID, maxtreeID[(n-i+samplesize)])
      i=i-samplesize #Adjust remaining number of trees to sample
    }
  }
  treesel[treesel[,ID] %in% maxtreeID & !is.na(treesel$tosample),'tosample'] <-  paste0(treesel[treesel[,ID] %in% maxtreeID & !is.na(treesel$tosample),'tosample'], 'max',var,'_')
  treesel[treesel[,ID] %in% maxtreeID & is.na(treesel$tosample),'tosample'] <-  paste0('max',var,'_')
}

############### Randomly sample which bin to sample from and randomply sample within that bin ###############
nrest = 40 -(nrow(treesel[!is.na(treesel$tosample),])+length(unique(fielddatasel$`#`))) #Remaining number of trees to sample after picking the min and max ones
binunique <- unique(treesel$bincomb)[!(unique(treesel$bincomb) %in% unique(treesel[!(is.na(treesel$tosample) & is.na(treesel$sampled)),'bincomb']))] #Unique variable bin combinations not already sampled
binsample <- sample(binunique, nrest, replace=F) # Randomly draw which bins will be sampled

for (b in binsample) { #For each bin to sample
  print(b)
  sampledID <- NULL
  while (is.null(sampledID)) { 
    sampledID <- sample(treesel[treesel$bincomb==b & !(treesel$UNITID %in% fielddatasel$Tree.ID.resampling),'UNITID'], 1,replace=F) #Sample one tree
    treesu100 <- as.character(distmelt_lim[distmelt_lim$UNITID %in% treesel[!(is.na(treesel$tosample) & is.na(treesel$sampled)),'UNITID'], 'UNITID.y']) #List of trees within 100-m of already-sampled trees
    if (sampledID %in% treesu100) sampledID <- NULL  #Make sure that it's not within maxdist of any other sampled tree
    #Otherwise, re-sample within that bin
  }
  treesel[treesel$UNITID==sampledID,'tosample'] <- 'random'
}

#Check distribution of sampled trees
ggplot(treesel, aes(x=bing_bin, y=AADT_bin, color=factor(spdlm_bin))) + 
  geom_jitter()  +
  labs(x='Congestion bin', y='AADT bin') + 
  geom_point(data=treesel[!(is.na(treesel$tosample) & is.na(treesel$sampled)),], color='black',size=3)

save.image(paste0('random_sampling', Sys.Date()))
write.csv(treesel[!is.na(treesel$tosample),], 'preliminary_tree_sample_20180730.csv', row.names=F)