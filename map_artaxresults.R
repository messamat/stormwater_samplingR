#options(warn=-1)
library(data.table)
library(xlsx)
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
library(colorspace)
library(plotly)
library(listviewer)
data(periodicTable)

rootdir <- "C:/Mathis/ICSL/stormwater" #UPDATE
resdir <- file.path(rootdir, "results")
datadir <- file.path(rootdir, "data")

#Element colors: http://jmol.sourceforge.net/jscolors/

#Define functions
deconvolution_import <- function(artaxdir, idstart, idend) {
  #Collate all results from ARTAX deconvolution outputs
  tablist <- list.files(artaxdir)
  for (i in 1:length(tablist)) {
    res <- read.csv(file.path(artaxdir,tablist[i]))
    res$XRFID <- str_extract(substr(tablist[i], idstart, idend), "[0-9A-B_]+")
    if (i==1){
      deconvres <- res
    } else{
      deconvres <- rbind(deconvres, res)
    }
  }
  return(setDT(deconvres))
}

vispal <- function(dat) {
  #Create a color palette graph
  barplot(rep(1,length(dat)), col=dat)
}

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

#Import and format field XRF deconvolution results
fieldres <- deconvolution_import(file.path(datadir, 'XRF20181210/PostBrukerCalibration/deconvolutionresults_XRF12_245_20180827'),
                                 idstart = 41, idend= 43)

#Import raw XRF lab data
labxrf_list <- data.table(filename = grep('\\.csv$', list.files(file.path(datadir, 'XRF20181210/PelletMeasurements')) ,value=T))
for (i in labxrf_list$filename) {
  print(i)
  xrfrec <- read.csv(file.path(datadir, 'XRF20181210/PelletMeasurements', i))
  labxrf_list[filename == i, `:=`(cps = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Valid Count Last Packet',])),
                                  c_cum = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Valid Accumulated Counts',])),
                                  duration = as.numeric(as.character(xrfrec[rownames(xrfrec) == 'Duration Time',]))
  )]
}
setDT(labxrf_list)[, `:=`(c_cum_ps = c_cum/duration,
                          SiteID = regmatches(labxrf_list$filename, regexpr('[0-9A-Z]+[_][0-9]', labxrf_list$filename)))]

# problem_recs <- setDT(labxrf_list)[SiteID %in% c('3A_1', '3A_2','6B_2','7A_1', '7A_2','19A_1', '19A_2', '20A_1',
#                                                  '20B_1', '20B_2','22A_1', '22A_2', '28B_2', '32A_1', '33B_1', 
#                                                  '33B_2','34A_1','36A_2', '37A_1', '37A_2','40A_1', '42A_1', 
#                                                  '42B_1','42B_2','44A_1', '44A_2', '45A_1', '45A_2','46A_1','46A_2',
#                                                  '48A_1','49A_1', '49A_2','53A_1', '53A_2', '54B_2','55B_1', '55B_2',
#                                                  '57A_1', '57A_2','59A_1', '59A_2', '61A_1', '61A_2'),]
labxrf_list[c_cum_ps > 90000,][!(labxrf_list[c_cum_ps > 90000, SiteID] %in% problem_recs$SiteID)]

#Import and format lab XRF deconvolution results
labres <- deconvolution_import(file.path(datadir, 'XRF20181210/PelletMeasurements/deconvolutionresults_labXRF_20181212'),
                               idstart = 29, idend= 33)
#Import ICP-OES data
ICPdat <- read.xlsx(file.path(datadir, '/ICP_OES_20181206/mathis_ICP_120618-2_edit20181210.xlsx'), 1)
ICPthresholds <- read.xlsx(file.path(datadir, '/ICP_OES_20181206/mathis_ICP_120618-2_edit20181210.xlsx'), 2)

#Import GIS data
trees <- readOGR(dsn = file.path(datadir, 'field_data'), layer = 'sampling_sites_edit') 
heat_AADT <- raster(file.path(resdir, 'heataadt_int'))
heat_spdlm <- raster(file.path(resdir, 'heatspdlm_int'))
heat_bing <- raster(file.path(resdir, 'heatbing_int')) #Replace with heatbing_index for all PugetSound to be able to include Highway 2 sample

#--------------------------------------------------------------------------------------
#Format XRF data

#Cast while summing net photon counts across electron transitions
fieldrescast <- dcast(setDT(fieldres), XRFID~Element, value.var='Net', fun.aggregate=sum) 
fieldrescast[, XRFID := as.numeric(gsub('[_]', '', XRFID))]
labrescast <- dcast(setDT(labres), XRFID~Element, value.var='Net', fun.aggregate=sum)
#Average lab XRF results over multiple measurements for a given pellet
labrescast <- labrescast[, `:=`(SiteID = gsub('[A-B].*', '', XRFID),
                                Pair = gsub('[0-9_]+', '', XRFID),
                                XRFID = NULL)]
labrescast <- labrescast[, lapply(.SD, mean), by = .(SiteID, Pair)]

#Normalize data by Rhodium photon count
fieldrescastnorm <- fieldrescast[, lapply(.SD, function(x) {x/Rh}), by = XRFID]
labrescastnorm <- labrescast[, lapply(.SD, function(x) {x/Rh}), by = .(SiteID, Pair)]
  
#Merge datasets
fieldt <- setDT(fieldata_format)[fieldrescastnorm, on='XRFID']
labdt <- setDT(fieldata_format)[labrescastnorm, on =  .(SiteID, Pair)]

#Compute average
field_artaxstats <- fieldt[, sapply(.SD, function(x) list(mean=mean(x, na.rm=T),
                                                          sd=sd(x, na.rm=T),
                                                          range=max(x,na.rm=TRUE)-min(x,na.rm=TRUE))), by=c('SiteID','Pair'), .SDcols = 29:51]
setnames(field_artaxstats, c('SiteID', 'Pair', 
                             paste(rep(colnames(fieldt)[29:51],each=3),
                                   c('mean', 'sd', 'range'), sep='_')))

#Write out mean values
field_artaxmean <- field_artaxstats[, .SD, .SDcols = c(1,2, grep('mean', colnames(field_artaxstats)))] 
write.dbf(field_artaxmean, 'field_artaxmean_20180827.dbf')

#--------------------------------------------------------------------------------------
#Assess within-tree variability of XRF measurements 
artaxres <- melt(field_artaxstats, id.vars=c('SiteID','Pair'), variable.name='Elem_stats')
artaxres[, `:=`(Elem = sub('(.*)[_](.*)', '\\1', Elem_stats),
                stats = sub('(.*)[_](.*)', '\\2', Elem_stats))] #Replaces the whole match with the first group in the regex pattern
artaxres <- dcast(artaxres, SiteID+Pair+Elem~stats, value.var = 'value')
artaxres[, cv := sd/mean]
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
field_artaxmeansite <- field_artaxmean[,lapply(.SD,mean,na.rm=TRUE), by=c('SiteID'), .SDcols=3:25]
artaxsdsite <- field_artaxmean[,lapply(.SD,sd,na.rm=TRUE), by=c('SiteID'), .SDcols=3:25]
artaxressite <- merge(melt(field_artaxmeansite, id.vars='SiteID', variable.name='Elem', value.name='mean'),
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
# Check whether tree species mattered at all

#--------------------------------------------------------------------------------------
#Relate XRF field data to ICP results

#Format data 
ICPdat[ICPdat == 'TR'] <- 0
ICPdat[ICPdat == 'ND'] <- 0
ICPdat <- setDT(ICPdat)[!(SAMPLE.SET %in% c('BLK', 'QC-1', NA)),]
ICPdat <- ICPdat[, lapply(.SD, as.character), by=SAMPLE.SET]
ICPdat <- ICPdat[, lapply(.SD, as.numeric), by=SAMPLE.SET]

#Set aside and compute mean over duplicate measurements
ICPdup <- ICPdat[substr(ICPdat$SAMPLE.SET, 1,3) %in% substr(grep('DUP', ICPdat$SAMPLE.SET, value=T), 1,3),]
ICPdat[, SAMPLE.SET := substr(SAMPLE.SET, 1,3)] 
ICPmean <- ICPdat[, lapply(.SD, mean), by= SAMPLE.SET]

#Melt dataset to merge with XRF data
ICPmelt <- melt(setDT(ICPdat), id.vars = 'SAMPLE.SET', variable.name = 'Elem', value.name = 'ICP')
ICPmelt[, `:=`(SiteID = gsub('[A-Z]', '', SAMPLE.SET),
               Pair = gsub('[0-9]', '', SAMPLE.SET))]

#Merge with XRF data
ICPmerge <- ICPmelt[artaxres, on = .(SiteID, Pair, Elem)]
unique(ICPmelt$Elem)[!(unique(ICPmelt$Elem) %in% unique(artaxres$Elem))]
unique(artaxres$Elem)[!(unique(artaxres$Elem) %in% unique(ICPmelt$Elem))]

#Plot comparison
ICPfield_plot <- ggplot(ICPmerge[!is.na(ICP),], aes(x=mean, y=ICP, label=SiteID)) + 
  geom_point() + 
  geom_smooth(span=0.9) +
  facet_wrap(~Elem, scales='free') +
  theme_classic()

plotly_json(ICPfield_plot)
ggplotly(ICPfield_plot) %>%
  style(hoverinfo = "none", traces = 15:42)

ICPmerge[!is.na(ICP),round(summary(lm(ICP~mean))$adj.r.squared, 2), by=Elem]

#--------------------------------------------------------------------------------------
#Relate XRF field data to XRF lab data


#Bayesian deconvolution of 7A2 matches a lot better


#--------------------------------------------------------------------------------------
#Relate XRF data to pollution predictors 

#Import shapefile of sites
colnames(trees@data)[1] <- 'SiteID'

#Join field data to site shapefile
treesxrf <- merge(trees, field_artaxmean, by=c('SiteID','Pair'))
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
metal_select <- c('Br','Co','Cr','Cu','Fe','Mn','Ni','Pb','Sr','Ti','V','Zn')
#metal_select <- c('Fe','Zn','Co','Cr') #for Phil 2018/10/02 'goodresults'
#metal_select <- c('Ca','K','Sr','P') #for Phil 2018/10/02 'randomresults'

treesxrf_melt <- melt(treesxrf@data, id.vars=colnames(treesxrf@data)[!(colnames(treesxrf@data) %in% metal_select)], variable.name = 'Elem')
treesxrf_melt <- merge(treesxrf_melt, periodicTable[,c('symb','name')], by.x='Elem', by.y='symb')

#Check distribution of aadt, bing, and spdlm indices
qplot(sqrt(treesxrf$AADT))
qplot(sqrt(treesxrf$spdlm))
qplot(sqrt(treesxrf$bing))
setDT(treesxrf_melt)[,valuestand := (value-min(value, na.rm=T))/(max(value)-min(value)), by=Elem]

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
  stat_poly_eq(aes(label = paste(..rr.label..)),
               label.x.npc = "right", label.y.npc = 0.02,
               formula = y ~ x, parse = TRUE, size = 5)+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x.npc = 'right', label.y.npc = 0.12, size = 3) +
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
# Relate ICP data to pollution predictors 
treesICP <- setDT(treesxrf_melt)[ICPmelt, on=.(SiteID, Pair, Elem)]

#
ggplot(treesICP, aes(x=bing, y=ICP)) + 
  geom_point() + 
  labs(x='Speed limit index', y='Concentration (ug/g)') + 
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


#Speed limit
ggplot(treesICP, aes(x=spdlm, y=ICP)) + 
  geom_point() + 
  labs(x='Speed limit index', y='Concentration (ug/g)') + 
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

#Traffic volume
ggplot(treesICP, aes(x=AADT, y=ICP)) + 
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


#--------------------------------------------------------------------------------------
# Output overall variability table 
elemvar <- data.frame(Elem=colnames(fieldrescastnorm[,-1]),
                      mean=colMeans(fieldrescastnorm[,-1]),
                      min =t(setDT(fieldrescastnorm[,-1])[ , lapply(.SD, min)]),
                      max =t(setDT(fieldrescastnorm[,-1])[ , lapply(.SD, max)]),
                      sd= t(setDT(fieldrescastnorm[,-1])[ , lapply(.SD, sd)]))

#Generate table of statistics
dt$SiteIDPair <- with(dt, paste0(SiteID,Pair))
allelem <- colnames(dt) %in% periodicTable$symb #Get all column which correspond to an element in the periodic table
#Apply to each element:
elemwise_comps <- sapply(colnames(dt)[allelem], function(x) {
  #Analyze error within tree samples
  withintree <- summary(aov(as.formula(paste(x, "~SiteIDPair")), data=dt))  
  #Analyze error among tree samples within sites
  withinsite <- summary(aov(as.formula(paste(x, "~SiteID")), data=dt)) 
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


