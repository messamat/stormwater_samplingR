library(data.table)
library(plyr)
library(stringr)
library(foreign)
library(rgdal)
library(raster)
library(ggplot2)
library(ggpmisc)
library(PeriodicTable)
data(periodicTable)

rootdir= "C:/Mathis/ICSL/stormwater" #UPDATE
setwd(file.path(rootdir, "results")) 
datadir = file.path(rootdir, "data")

#Import and format field data TO DO: make faster
fielddata <- read.csv(file.path(datadir,"field_data/field_data_raw_20180808_edit.csv"))
colnames(fielddata)[1] <- 'SiteID'
fielddata$SiteID <- as.character(fielddata$SiteID)
sel <- fielddata[!is.na(fielddata$XRFmin) & !is.na(fielddata$SiteID),]
fieldata_format <- data.frame()
for (row in seq(1,nrow(sel))) {
  extract <- sel[row,]
  for (xrf in seq(sel[row,'XRFmin'], sel[row,'XRFmax'])){
    #print(xrf)
    extract$XRFID <- xrf
    fieldata_format <- rbind(fieldata_format, extract)
  }
}

#Import and format ARTAX deconvolution results
artaxdir <- file.path(datadir, 'XRF20180808/PostBrukerCalibration/deconvolutionresults_XRF12_245_20180827')
tablist <- list.files(artaxdir)
for (i in 1:length(tablist)) {
  res <- read.csv(file.path(artaxdir,tablist[i]))
  res$XRFID <- str_extract(substr(tablist[i],41,43), "\\d+")
  if (i==1){
    deconvres <- res
  } else{
    deconvres <- rbind(deconvres, res)
  }
}
 <- dcast(deconvres, XRFID~Element, value.var='Net', fun.aggregate=sum)
#Normalize data by Rhodium photon count 
deconvrescastnorm <- cbind(XRFID=deconvrescast$XRFID, sweep(deconvrescast[,-1], MARGIN=1, FUN='/', STATS=deconvrescast$Rh))
#Merge datasets
df <- merge(fieldata_format, deconvrescastnorm, by='XRFID')
#Compute average
colnames(df)
artaxmean <- setDT(df)[,lapply(.SD,mean,na.rm=TRUE), by=c('SiteID','Pair'), .SDcols=29:51]
#Compute sd
artaxsd <- setDT(df)[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID','Pair'), .SDcols=29:51]
#Compute range
artaxrange <- setDT(df)[,lapply(.SD,function(x) max(x,na.rm=TRUE)-min(x,na.rm=TRUE)), by=c('SiteID','Pair'), .SDcols=29:51]
#Write out
write.dbf(fielddata_artaxmean, 'artaxmean_20180827.dbf')

#----------------- Assess variability of XRF measurements within trees
#Format data
artaxres <- merge(melt(artaxmean, id.vars=c('SiteID','Pair'), variable.name='Elem', value.name='mean'),
      melt(artaxsd, id.vars=c('SiteID','Pair'), variable.name='Elem', value.name='sd'),
      by=c('SiteID','Pair','Elem'))
artaxres$cv <- with(artaxres, sd/mean)
artaxres <- merge(artaxres, melt(artaxrange, id.vars=c('SiteID', 'Pair'), variable.name='Elem', value.name='range'),
      by=c('SiteID','Pair','Elem'))
artaxres <- merge(artaxres, periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

#Plot coefficient of variation distributions for every element
cvmean_label <- as.data.frame(artaxres[!(artaxres$Elem %in% c('Rh','Pd','Ar')),
                                       paste0('Mean CV: ',format(mean(cv),digits=2)),
                                       by=name])
ggplot(artaxres[!(artaxres$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_histogram()+
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~name, scales='free_x') + 
  geom_text(data=cvmean_label, 
            mapping=aes(x=Inf,y=Inf, label=V1, hjust=1, vjust=1))  + 
  theme_classic() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14))

#Plot relationship between elemental concentration and cv
ggplot(artaxres[!(artaxres$Elem %in% c('Rh','Pd','Ar')),], aes(x=mean, y=cv, color=name)) + 
  geom_point()+
  labs(x='Mean photon count (normalized)', y='Coefficient of variation') + 
  facet_wrap(~name, scales='free') +
  theme_classic() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14))

#----------------- Assess variability of XRF measurements among trees within a site







#------------------RELATE DATA TO POLLUTION PREDICTORS
#Import shapefile of sites
trees <- readOGR(dsn = file.path(datadir, 'field_data'), layer = 'sampling_sites_edit') 
colnames(trees@data)[1] <- 'SiteID'

#Join field data to site shapefile
treesxrf <- merge(trees, artaxmean, by=c('SiteID','Pair'))
treesxrf <- treesxrf[!is.na(treesxrf$Fe) | treesxrf$SiteID==1,]

#Define pollution predictor layers
heat_AADT <- raster(file.path(resdir, 'heataadt_int'))
heat_spdlm <- raster(file.path(resdir, 'heatspdlm_int'))
heat_bing <- raster(file.path(resdir, 'heatbing_int')) #Replace with heatbing_index for all PugetSound to be able to include Highway 2 sample
  
#Project site shapefile to same coordinate system as predictors
treesxrf <- spTransform(treesxrf, crs(heat_bing))
  
#Extract values of pollution predictors at sites
treesxrf$AADT <- extract(heat_AADT,treesxrf,method='simple')
treesxrf$spdlm <- extract(heat_spdlm,treesxrf,method='simple')
treesxrf$bing <- extract(heat_bing,treesxrf,method='simple')

treesxrf[is.na(treesxrf$AADT),'AADT'] <- 0
treesxrf[is.na(treesxrf$spdlm),'spdlm'] <- 0
treesxrf[is.na(treesxrf$bing),'bing'] <- 0

#Format data
treesxrf <- treesxrf[as.numeric(as.character(treesxrf$SiteID)) < 52,]
allelem <- colnames(treesxrf@data) %in% periodicTable$symb
metal_select <- c('Br','Co','Cr','Cu','Fe','Mn','Ni','Pb','Sr','Ti','V','Zn')
treesxrf_melt <- melt(treesxrf@data, id.vars=colnames(treesxrf@data)[!(colnames(treesxrf@data) %in% metal_select)])
treesxrf_melt <- merge(treesxrf_melt, periodicTable[,c('symb','name')], by.x='variable', by.y='symb')


#Plot Bing Congestion index vs. XRF data
ggplot(treesxrf_melt, aes(x=bing, y=value)) + 
  geom_point() + 
  labs(x='Bing congestion index', y='Photon count (normalized)') + 
  facet_wrap(~name, nrow=3, scales='free_y') +
  geom_smooth(method='lm', formula = y ~ x, se = T) + 
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic()

#Plot traffic volume index vs. XRF data
ggplot(treesxrf_melt, aes(x=AADT, y=value)) + 
  geom_point() + 
  labs(x='Traffic volume index', y='Photon count (normalized)') + 
  facet_wrap(~name, nrow=3, scales='free_y') +
  geom_smooth(method='lm', formula = y ~ x, se = T) + 
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic()

#Plot speed limit index vs. XRF data
ggplot(treesxrf_melt, aes(x=spdlm, y=value)) + 
  geom_point() + 
  labs(x='Speed limit index', y='Photon count (normalized)') + 
  facet_wrap(~name, nrow=3, scales='free_y') +
  geom_smooth(method='lm', formula = y ~ x, se = T) + 
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 3)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic()


#Multiple linear regressions
Znlm <- lm(Zn~bing+sqrt(AADT), data=treesxrf)
summary(Znlm)
qplot(x=predict(Znlm), y=treesxrf$Zn) + geom_smooth(method='lm') + 
  geom_text(aes(label=treesxrf$SiteID))

Felm <- lm(Fe~bing+sqrt(AADT)+spdlm, data=treesxrf)
summary(Felm)
qplot(x=predict(Felm), y=treesxrf$Fe) + geom_smooth(method='lm')

Culm <- lm(Cu~sqrt(AADT), data=treesxrf)
summary(Culm)
qplot(x=predict(Culm), y=treesxrf$Cu) + geom_smooth(method='lm')+ 
  geom_text(aes(label=treesxrf$SiteID))

Pblm <- lm(Pb~bing+sqrt(AADT)+spdlm, data=treesxrf)
summary(Pblm)
qplot(x=predict(Pblm), y=treesxrf$Pb) + geom_smooth(method='lm')+ 
  geom_text(aes(label=treesxrf$SiteID))

#Develop multilinear Kriging model
