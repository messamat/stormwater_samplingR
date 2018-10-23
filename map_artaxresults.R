#options(warn=-1)
library(data.table)
library(plyr)
library(stringr)
library(foreign)
library(rgdal)
library(raster)
library(ggplot2)
library(ggpmisc)
library(PeriodicTable)
library(kableExtra)
library(dplyr)
library(hexbin)
data(periodicTable)

rootdir <- "C:/Mathis/ICSL/stormwater" #UPDATE
resdir <- file.path(rootdir, "results")
datadir <- file.path(rootdir, "data")

vispal <- function(dat) {
  barplot(rep(1,length(dat)), col=dat)
}

#Element colors: http://jmol.sourceforge.net/jscolors/

#--------------------------------------------------------------------------------------
# Import data 

#Import and format field data 
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

#Import GIS data
trees <- readOGR(dsn = file.path(datadir, 'field_data'), layer = 'sampling_sites_edit') 
heat_AADT <- raster(file.path(resdir, 'heataadt_int'))
heat_spdlm <- raster(file.path(resdir, 'heatspdlm_int'))
heat_bing <- raster(file.path(resdir, 'heatbing_int')) #Replace with heatbing_index for all PugetSound to be able to include Highway 2 sample

#--------------------------------------------------------------------------------------
#Format XRF data

deconvrescast <- dcast(deconvres, XRFID~Element, value.var='Net', fun.aggregate=sum)
#Normalize data by Rhodium photon count 
deconvrescastnorm <- cbind(XRFID=deconvrescast$XRFID, sweep(deconvrescast[,-1], MARGIN=1, FUN='/', STATS=deconvrescast$Rh))
#Merge datasets
df <- merge(fieldata_format, deconvrescastnorm, by='XRFID')
#Compute average
artaxmean <- setDT(df)[,lapply(.SD,mean,na.rm=TRUE), by=c('SiteID','Pair'), .SDcols=29:51]
#Compute sd
artaxsd <- setDT(df)[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID','Pair'), .SDcols=29:51]
#Compute range
artaxrange <- setDT(df)[,lapply(.SD,function(x) max(x,na.rm=TRUE)-min(x,na.rm=TRUE)), by=c('SiteID','Pair'), .SDcols=29:51]
#Write out
write.dbf(artaxmean, 'artaxmean_20180827.dbf')

#--------------------------------------------------------------------------------------
#Assess within-tree variability of XRF measurements 

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
  geom_density()+
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~name) + 
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

#--------------------------------------------------------------------------------------
#Assess within-site variability of XRF measurements 

artaxmeansite <- artaxmean[,lapply(.SD,mean,na.rm=TRUE), by=c('SiteID'), .SDcols=3:25]
artaxsdsite <- artaxmean[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID'), .SDcols=3:25]
artaxressite <- merge(melt(artaxmeansite, id.vars='SiteID', variable.name='Elem', value.name='mean'),
                  melt(artaxsdsite, id.vars='SiteID', variable.name='Elem', value.name='sd'),
                  by=c('SiteID','Elem'))
artaxressite$cv <- with(artaxressite, sd/mean)
artaxressite <- merge(artaxressite, periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

ggplot(artaxressite[!(artaxressite$Elem %in% c('Rh','Pd','Ar')),], aes(x=cv, fill=name)) + 
  geom_density()+
  labs(x='Coefficient of variation', y='Count') + 
  facet_wrap(~name) + 
  theme_bw() + 
  theme(strip.background = element_rect(color = 'white', fill = 'lightgrey'),
        strip.text = element_text(size=14),
        panel.border = element_blank(),
        axis.line = element_line(color='black'))

#--------------------------------------------------------------------------------------
#Relate XRF data to pollution predictors 

#Import shapefile of sites
colnames(trees@data)[1] <- 'SiteID'

#Join field data to site shapefile
treesxrf <- merge(trees, artaxmean, by=c('SiteID','Pair'))
treesxrf <- treesxrf[!is.na(treesxrf$Fe) | treesxrf$SiteID==1,]

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
#metal_select <- c('Br','Co','Cr','Cu','Fe','Mn','Ni','Pb','Sr','Ti','V','Zn')
#metal_select <- c('Fe','Zn','Co','Cr') #for Phil 2018/10/02 'goodresults'
metal_select <- c('Ca','K','Sr','P') #for Phil 2018/10/02 'randomresults'

treesxrf_melt <- melt(treesxrf@data, id.vars=colnames(treesxrf@data)[!(colnames(treesxrf@data) %in% metal_select)])
treesxrf_melt <- merge(treesxrf_melt, periodicTable[,c('symb','name')], by.x='variable', by.y='symb')

#Check distribution of aadt, bing, and spdlm indices
qplot(sqrt(treesxrf$AADT))
qplot(sqrt(treesxrf$spdlm))
qplot(sqrt(treesxrf$bing))
setDT(treesxrf_melt)[,valuestand := (value-min(value, na.rm=T))/(max(value)-min(value)), by=variable]

rPal <- colorRampPalette(c('#fff5f0','#cb181d'))
bPal <- colorRampPalette(c('#deebf7','#2171b5'))
treesxrf_melt$valuefact <- factor(treesxrf_melt$valuestand, levels= unique(as.character(treesxrf_melt$valuestand[order(treesxrf_melt$valuestand,
                                                                                                                       treesxrf_melt$bing)])))

col1 <- t(col2rgb(bPal(100)[as.numeric(cut(treesxrf_melt$valuestand,breaks = 100))]))
colnames(col1) <- paste0(colnames(col1),'_1')
col2 = t(col2rgb(rPal(100)[as.numeric(cut(treesxrf_melt$bing,breaks = 100))]))
colnames(col2) <- paste0(colnames(col2),'_2')
cols <- data.table(cbind(col1, col2))

colsmix <- cols[,mixcolor(0.5, RGB(red_1, green_1,blue_1), RGB(red_2, green_2,blue_2))]
treesxrf_melt$colsmix <- rgb(colsmix@coords, maxColorValue=255)


#colsmix <- factor(colsmix, levels= unique(colsmix[order(treesxrf_melt$valuestand,
#                                                        treesxrf_melt$bing)]))

vispal(rgb(col1, maxColorValue=255))
vispal(rgb(col2, maxColorValue=255))
#vispal(as.character(colsmix))

#Plot Bing Congestion index vs. XRF data
bing_xrf <- ggplot(treesxrf_melt, aes(x=bing, y=as.numeric(as.character(valuefact)))) + 
  geom_point(aes(color=colsmix), size=3, alpha=1/4) +
  scale_color_identity() +
  labs(x='Traffic congestion (Bing index)', y='Concentration in moss (normalized photon count)') + 
  facet_wrap(~name, nrow=3, scales='free_y') +
  geom_smooth(method='lm', se = T, color='black') + 
  # stat_poly_eq(aes(label = paste(..rr.label..)), 
  #              label.x.npc = "right", label.y.npc = 0.02,
  #              formula = y ~ x, parse = TRUE, size = 5)+
  # stat_fit_glance(method = 'lm',
  #                 method.args = list(formula = y ~ x),
  #                 geom = 'text',
  #                 aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
  #                 label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        text = element_text(size=16),
        axis.title = element_blank())
bing_xrf


png(file.path(resdir, 'prelim_badXRFresults.png'), width=8, height=8, units = 'in', res=400)
#pdf(file.path(resdir, 'prelim_badXRFresults.pdf'), width=8, height=8)
bing_xrf
dev.off()

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
qplot(x=predict(Znlm), y=treesxrf$Zn) + geom_smooth(method='lm') 
#+ geom_text(aes(label=treesxrf$SiteID))

Felm <- lm(Fe~bing+sqrt(AADT)+spdlm, data=treesxrf)
summary(Felm)
qplot(x=predict(Felm), y=treesxrf$Fe) + geom_smooth(method='lm')

Culm <- lm(Cu~sqrt(AADT), data=treesxrf)
summary(Culm)
qplot(x=predict(Culm), y=treesxrf$Cu) + geom_smooth(method='lm') 

Pblm <- lm(Pb~bing+sqrt(AADT)+spdlm, data=treesxrf)
summary(Pblm)
qplot(x=predict(Pblm), y=treesxrf$Pb) + geom_smooth(method='lm')

#--------------------------------------------------------------------------------------
# Output overall variability table 
elemvar <- data.frame(Elem=colnames(deconvrescastnorm[,-1]),
                      mean=colMeans(deconvrescastnorm[,-1]),
                      min =t(setDT(deconvrescastnorm[,-1])[ , lapply(.SD, min)]),
                      max =t(setDT(deconvrescastnorm[,-1])[ , lapply(.SD, max)]),
                      sd= t(setDT(deconvrescastnorm[,-1])[ , lapply(.SD, sd)]))

#Generate table of statistics
df$SiteIDPair <- with(df, paste0(SiteID,Pair))
allelem <- colnames(df) %in% periodicTable$symb #Get all column which correspond to an element in the periodic table
#Apply to each element:
elemwise_comps <- sapply(colnames(df)[allelem], function(x) {
  #Analyze error within tree samples
  withintree <- summary(aov(as.formula(paste(x, "~SiteIDPair")), data=df))  
  #Analyze error among tree samples within sites
  withinsite <- summary(aov(as.formula(paste(x, "~SiteID")), data=df)) 
  #Regress pollutant predictors against XRF data for each element
  bing_lm <- summary(lm(as.formula(paste(x, "~bing")), data=treesxrf))
  aadt_lm <- summary(lm(as.formula(paste(x, "~sqrt(AADT)")), data=treesxrf))
  spdlm_lm <- summary(lm(as.formula(paste(x, "~sqrt(spdlm)")), data=treesxrf))
    
  return(c(tree_error=withintree[[1]]$`Sum Sq`[2]/sum(withintree[[1]]$`Sum Sq`),
           site_error=(withinsite[[1]]$`Sum Sq`[2]-withintree[[1]]$`Sum Sq`[2])/sum(withinsite[[1]]$`Sum Sq`),
           total_error = withinsite[[1]]$`Sum Sq`[2]/sum(withinsite[[1]]$`Sum Sq`),
           bing_r2 = bing_lm$adj.r.squared, 
           bing_sig = bing_lm$coefficients[8],
           aadt_r2 = aadt_lm$adj.r.squared,
           aadt_sig = aadt_lm$coefficients[8],
           spdlm_r2 = spdlm_lm$adj.r.squared,
           spdlm_sig = spdlm_lm$coefficients[8]
           ))
})
elemvar <- cbind(elemvar, t(elemwise_comps))
#Get rid of useless elements (signals from XRF Tracer detector and collimator, not in sample)
elemvar <- elemvar[!(elemvar$Elem %in% c('Ar','Rh','Pd')),]
#Get full element names
elemvar <- merge(elemvar, periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

#Format and output table
options(knitr.table.format = "html") 
elemvar[order(-elemvar$bing_r2),] %>%
  mutate(
    congestion_r2 = cell_spec(format(bing_r2, digits=3), bold = ifelse(bing_sig < 0.05, TRUE, FALSE)),
    volume_r2 = cell_spec(format(aadt_r2, digits=3), bold = ifelse(aadt_sig < 0.05, TRUE, FALSE)),
    speedlimit_r2 = cell_spec(format(spdlm_r2, digits=3), bold = ifelse(spdlm_sig < 0.05, TRUE, FALSE))
  ) %>%
  select(-c(9:14)) %>%
  kable(format = "html", escape = F) %>%
  kable_styling("striped", full_width = F) 

#%>%
#  save_kable(file.path(resdir, 'XRFrestab_20180828.doc'), self_contained=T)

#------------- Develop multilinear bayesian Kriging model  --------------------


